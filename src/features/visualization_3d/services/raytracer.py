import numpy as np
import math
import time


# ---------------------------------------------------------------------------
# Module-level worker function (must be top-level for spawn-mode pickling)
# ---------------------------------------------------------------------------

def _raytrace_chunk(args):
    """
    Render a single tile of the image.

    Parameters (packed in a single tuple for compatibility with
    ``ParallelExecutor.map``):
        x0, y0, x1, y1       – pixel bounds of this chunk
        spheres_arr           – np.ndarray of shape (N, 7)  [x, y, z, r, R, G, B]
        triangles_arr         – np.ndarray of shape (M, 12) [v0..v2, R, G, B]
        L                     – np.ndarray(3,) normalised light direction
        bg_col                – np.ndarray(3,) background colour

    Returns:
        (x0, y0, x1, y1, rgba)  where *rgba* is uint8 array (h_c, w_c, 4).
    """
    (x0, y0, x1, y1, spheres_arr, triangles_arr, L, bg_col) = args

    w_c = x1 - x0
    h_c = y1 - y0

    # Setup rays – orthographic projection
    X, Y = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1))
    O_chunk = np.stack(
        [X.ravel(), Y.ravel(), np.full(w_c * h_c, -10000)], axis=-1
    ).astype(np.float32)

    num_rays = O_chunk.shape[0]

    t_min = np.full(num_rays, np.inf, dtype=np.float32)
    col_chunk = np.tile(bg_col, (num_rays, 1))
    normal_chunk = np.zeros((num_rays, 3), dtype=np.float32)
    hit_mask = np.zeros(num_rays, dtype=bool)

    # ---- Intersect Spheres ----
    if len(spheres_arr) > 0:
        for i in range(len(spheres_arr)):
            cx_s, cy_s, cz_s, r_s, r_col, g_col, b_col = spheres_arr[i]

            dx = O_chunk[:, 0] - cx_s
            dy = O_chunk[:, 1] - cy_s

            a = dx * dx + dy * dy
            disc = r_s * r_s - a

            hit = disc >= 0
            if np.any(hit):
                t_hit = cz_s - np.sqrt(disc[hit])
                closer = t_hit < t_min[hit]
                if np.any(closer):
                    valid_hits = np.zeros(num_rays, dtype=bool)
                    valid_hits[hit] = closer

                    t_min[valid_hits] = t_hit[closer]
                    col_chunk[valid_hits] = [r_col, g_col, b_col]
                    hit_mask[valid_hits] = True

                    px = O_chunk[valid_hits, 0]
                    py = O_chunk[valid_hits, 1]
                    pz = t_min[valid_hits]

                    normal_chunk[valid_hits, 0] = (px - cx_s) / r_s
                    normal_chunk[valid_hits, 1] = (py - cy_s) / r_s
                    normal_chunk[valid_hits, 2] = (pz - cz_s) / r_s

    # ---- Intersect Triangles (Möller–Trumbore) ----
    if len(triangles_arr) > 0:
        for i in range(len(triangles_arr)):
            v0 = triangles_arr[i, 0:3]
            v1 = triangles_arr[i, 3:6]
            v2 = triangles_arr[i, 6:9]
            t_col = triangles_arr[i, 9:12]

            e1 = v1 - v0
            e2 = v2 - v0
            pvec = np.array([-e2[1], e2[0], 0.0], dtype=np.float32)
            det = np.dot(e1, pvec)

            if abs(det) < 1e-8:
                continue
            inv_det = 1.0 / det

            tvec = O_chunk - v0
            u = np.sum(tvec * pvec, axis=1) * inv_det

            valid = (u >= 0) & (u <= 1)
            if not np.any(valid):
                continue

            qvec = np.zeros_like(tvec)
            qvec[:, 0] = tvec[:, 1] * e1[2] - tvec[:, 2] * e1[1]
            qvec[:, 1] = tvec[:, 2] * e1[0] - tvec[:, 0] * e1[2]
            qvec[:, 2] = tvec[:, 0] * e1[1] - tvec[:, 1] * e1[0]

            v = qvec[:, 2] * inv_det

            valid = valid & (v >= 0) & (u + v <= 1)
            if not np.any(valid):
                continue

            t = np.sum(e2 * qvec, axis=1) * inv_det
            valid = valid & (t > 1e-5) & (t < t_min)

            if np.any(valid):
                t_min[valid] = t[valid]
                col_chunk[valid] = t_col
                hit_mask[valid] = True

                n = np.cross(e1, e2)
                ni = n / np.linalg.norm(n)
                if ni[2] > 0:
                    ni = -ni
                normal_chunk[valid] = ni

    # ---- Shading ----
    if np.any(hit_mask):
        N = normal_chunk[hit_mask]
        C = col_chunk[hit_mask]

        diff = np.clip(np.sum(N * L, axis=1), 0.0, 1.0)

        amb = 0.3

        V = np.array([0, 0, -1], dtype=np.float32)
        H = L + V
        H = H / np.linalg.norm(H)
        spec = np.clip(np.sum(N * H, axis=1), 0.0, 1.0) ** 32

        final_c = C * (amb + diff[:, None] * 0.7) + spec[:, None] * 0.4
        col_chunk[hit_mask] = np.clip(final_c, 0.0, 1.0)

    # ---- Pack RGBA ----
    col_chunk_8 = (col_chunk * 255).astype(np.uint8)
    alpha_chunk = np.where(hit_mask, 255, 255).astype(np.uint8)

    rgba = np.hstack([col_chunk_8, alpha_chunk[:, None]])
    rgba = rgba.reshape((h_c, w_c, 4))

    return (x0, y0, x1, y1, rgba)


class NumpyRayTracer:
    """
    A pure-numpy offline ray-tracer.
    Supports rendering spheres (atoms) and triangles (protein cartoon meshes).
    """

    def __init__(self):
        pass

    def render(self, molecule, width, height, rot_x, rot_y, pan_x, pan_y, zoom, render_mode='ball_and_stick', progress_callback=None, executor=None):
        """
        Generates a ray-traced image.
        Returns an RGBA numpy array of shape (height, width, 4).
        """
        # --- 1. Gather Scene Data ---
        print("[RayTracer] Gathering scene data...")
        
        spheres = [] # list of [x, y, z, r, R, G, B]
        triangles = [] # list of [v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, R, G, B]
        
        # Camera transform matrix (rotation only, translation handled by ray origin)
        cos_x = math.cos(math.radians(rot_x))
        sin_x = math.sin(math.radians(rot_x))
        cos_y = math.cos(math.radians(rot_y))
        sin_y = math.sin(math.radians(rot_y))
        
        def transform(x, y, z):
            # Rotate Y
            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            # Rotate X
            y1 = y * cos_x - z1 * sin_x
            z2 = y * sin_x + z1 * cos_x
            return x1 * zoom, y1 * zoom, z2 * zoom

        cx = width / 2.0 + pan_x
        cy = height / 2.0 + pan_y

        # Gather atoms (Spheres)
        if render_mode not in ('cartoon', 'ribbon', 'backbone'):
            from src.features.visualization_3d.ui.mol_viewer_3d import DISPLAY_RADIUS
            for atom in molecule.atoms:
                if not atom.has_coords: continue
                # Transform to view space
                vx, vy, vz = transform(atom.x, atom.y, atom.z)
                # Shift offset to screen space centers (keep Z as depth)
                sx = cx + vx
                sy = cy - vy
                sz = vz
                
                radius = DISPLAY_RADIUS.get(atom.symbol, 0.35) * zoom * 0.6
                hex_c = atom.element.color.lstrip('#')
                r, g, b = tuple(int(hex_c[i:i+2], 16) for i in (0, 2, 4))
                
                spheres.append([sx, sy, sz, radius, r/255.0, g/255.0, b/255.0])
                
        elif render_mode == 'cartoon':
            # Gather Mesh
            from src.features.visualization_3d.services.protein_rendering import _cartoon_gen
            verts, tris, colors = _cartoon_gen.get_mesh(molecule)
            if verts is not None and len(verts) > 0:
                # Transform vertices
                t_verts = np.zeros_like(verts)
                # Apply same transform manually to all vertices
                vx = verts[:, 0] * zoom
                vy = verts[:, 1] * zoom
                vz = verts[:, 2] * zoom
                
                rx = vx * cos_y + vz * sin_y
                rz = -vx * sin_y + vz * cos_y
                ry = vy * cos_x - rz * sin_x
                rz2 = vy * sin_x + rz * cos_x
                
                t_verts[:, 0] = cx + rx
                t_verts[:, 1] = cy - ry
                t_verts[:, 2] = rz2
                
                v0 = t_verts[tris[:, 0]]
                v1 = t_verts[tris[:, 1]]
                v2 = t_verts[tris[:, 2]]
                
                # Combine
                tris_data = np.hstack([v0, v1, v2, colors])
                triangles = tris_data.tolist()

        spheres_arr = np.array(spheres, dtype=np.float32) if spheres else np.empty((0, 7), dtype=np.float32)
        triangles_arr = np.array(triangles, dtype=np.float32) if triangles else np.empty((0, 12), dtype=np.float32)
        
        # --- 2. Ray Tracing Setup ---
        img = np.zeros((height, width, 4), dtype=np.uint8)
        
        # Light ray (view space) -> coming from top-left-front
        lx, ly, lz = 0.5, -0.5, -1.0
        ll = math.sqrt(lx*lx + ly*ly + lz*lz)
        L = np.array([lx/ll, ly/ll, lz/ll], dtype=np.float32)
        
        # Process in chunks to save memory
        chunk_size = 64
        total_chunks = ((height + chunk_size - 1) // chunk_size) * ((width + chunk_size - 1) // chunk_size)

        print(f"[RayTracer] Rendering {width}x{height} image in {total_chunks} chunks...")

        bg_col = np.array([0.1, 0.1, 0.12], dtype=np.float32)

        # Build argument tuples for every chunk
        chunk_args = []
        for y0 in range(0, height, chunk_size):
            y1 = min(height, y0 + chunk_size)
            for x0 in range(0, width, chunk_size):
                x1 = min(width, x0 + chunk_size)
                chunk_args.append(
                    (x0, y0, x1, y1, spheres_arr, triangles_arr, L, bg_col)
                )

        # Dispatch: parallel via executor, or sequential fallback
        use_parallel = (
            executor is not None
            and total_chunks > 1
        )

        if use_parallel:
            print(f"[RayTracer] Using ParallelExecutor ({total_chunks} chunks)")
            try:
                results = executor.map(
                    _raytrace_chunk, chunk_args, timeout=300
                )
            except Exception as exc:
                print(f"[RayTracer] Parallel dispatch failed ({exc}), "
                      "falling back to sequential")
                results = [_raytrace_chunk(a) for a in chunk_args]
        else:
            results = [_raytrace_chunk(a) for a in chunk_args]

        # Assemble tiles into image
        for idx, (rx0, ry0, rx1, ry1, rgba) in enumerate(results):
            img[ry0:ry1, rx0:rx1] = rgba

            if progress_callback and (idx + 1) % 4 == 0:
                pct = ((idx + 1) / total_chunks) * 100
                if progress_callback(pct):
                    return None  # Cancelled

        return img
