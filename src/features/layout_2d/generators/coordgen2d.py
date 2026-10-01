"""
BKChem-style Professional 2D Coordinate Generator

Uses the deterministically robust OASA template-growth algorithm
to replace physics engines completely, solving topological singularities
and enforcing $120^circ$ constraints natively via regular polygons.
"""

from collections import deque, defaultdict
import math
import src.vendors.oasa.smiles as oasa_smiles
import src.vendors.oasa.coords_generator as oasa_cg
from src.vendors.oasa_bridge import domain_to_oasa_mol, oasa_prepare_for_layout

_DEBUG = False


class CoordinateGenerator2D:
    """Professional 2D coordinate generator wrapping OASA BKChem geometry."""
    
    BOND_LENGTH = 1.0
    
    def __init__(self, molecule, method='deterministic', force_regenerate=True):
        self.molecule = molecule
        self.coords = {}
        self.method = method
        self.force_regenerate = force_regenerate

    def generate(self):
        """Generate professional 2D coordinates natively utilizing OASA framework."""
        if not self.molecule or len(self.molecule.atoms) == 0:
            return {}
            
        # Check for existing coordinates - skip if force_regenerate is True
        if not self.force_regenerate:
            has_coords = all(hasattr(a, 'x2d') and a.x2d is not None for a in self.molecule.atoms)
            if has_coords:
                self.coords = {a.index: [a.x2d * self.BOND_LENGTH, a.y2d * self.BOND_LENGTH] for a in self.molecule.atoms}
                self._center_coords()
                return self.coords

        # Toggle: 'deterministic' uses the Golden Layout (SMILES Roundtrip)
        if self.method == 'deterministic':
            from src.features.layout_2d.methods.golden_layout import generate_golden_layout
            print("[DEBUG] Triggering Professional Golden Layout (SMILES Roundtrip)...") if _DEBUG else None
            self.coords = generate_golden_layout(self.molecule)
        else:
            # Traditional OASA calculation (older fallback)
            o_mol, atom_map = domain_to_oasa_mol(self.molecule)
            # Kekulize and clean up BEFORE layout — without this OASA
            # produces overlapping, non-deterministic layouts for
            # fused-ring systems because it can't orient aromatic rings.
            oasa_prepare_for_layout(o_mol)
            oasa_cg.calculate_coords(o_mol, bond_length=1.0, force=1)
            self.coords = {}
            for internal_idx, o_v in atom_map.items():
                x = float(o_v.x) if o_v.x is not None else None
                y = float(o_v.y) if o_v.y is not None else None
                if x is not None and y is not None:
                    self.coords[internal_idx] = [x, y]
            
        self._center_coords()
        
        # 1. Clean up overlapping chains / sterics from raw OASA output
        self._optimize_layout()

        # 2. Align the primary axis of the molecule horizontally using PCA
        self._align_pca()
        self._center_coords()
            
        # Cache internally back to the graph object
        for internal_idx in self.coords:
            x, y = self.coords[internal_idx]
            atom = self.molecule.get_atom(internal_idx)
            atom.x2d = x
            atom.y2d = y
            atom.z2d = 0.0

        return self.coords

    def _align_pca(self):
        """Align the longest molecular axis horizontally and canonicalize
        orientation so that re-running the upstream OASA layout (which is
        non-deterministic in its rigid rotation/reflection across calls)
        always produces the same final pose.

        The problem: OASA's vendored ``coords_generator`` uses Python sets
        keyed on vertex-object identity. Fresh atom object ids on each file
        read → different set iteration order → OASA picks rings in a
        different order → produces a layout that is the *same* (pairwise
        distances match to 1e-6, verified empirically) but under an
        arbitrary rigid motion. We fix this by fully canonicalizing the
        pose below using stable atom indices — not geometry-derived keys,
        which can tie in symmetric molecules.

        Canonical pose rule:
        1. Translate centroid to origin.
        2. PCA → rotate principal axis to +x.
        3. Flip x so that the atom with the SMALLEST INDEX among those with
           |x| > eps has x > 0. Atom indices are stable across reads.
        4. Flip y so that the atom with the SMALLEST INDEX among those with
           |y| > eps has y > 0.
        """
        import math
        if not self.coords or len(self.coords) < 2:
            return

        # --- 1. Centroid and PCA eigenvector ----------------------------
        cx = sum(c[0] for c in self.coords.values()) / len(self.coords)
        cy = sum(c[1] for c in self.coords.values()) / len(self.coords)

        cxx = cyy = cxy = 0.0
        for x, y in self.coords.values():
            dx, dy = x - cx, y - cy
            cxx += dx * dx
            cyy += dy * dy
            cxy += dx * dy

        b = -(cxx + cyy)
        c = cxx*cyy - cxy*cxy
        det = math.sqrt(max(0, b*b - 4*c))
        l1 = (-b + det) / 2.0

        vx_1, vy_1 = cxy, l1 - cxx
        vx_2, vy_2 = l1 - cyy, cxy
        if vx_1*vx_1 + vy_1*vy_1 > vx_2*vx_2 + vy_2*vy_2:
            vx, vy = vx_1, vy_1
        else:
            vx, vy = vx_2, vy_2

        norm = math.hypot(vx, vy)
        if norm < 1e-6:
            # Degenerate (single atom / colinear); nothing to align
            return
        vx /= norm
        vy /= norm

        # --- 2. Rotate principal axis to +x -----------------------------
        angle = -math.atan2(vy, vx)
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        rotated = {}
        for idx in self.coords:
            x, y = self.coords[idx]
            dx, dy = x - cx, y - cy
            rx = dx * cos_a - dy * sin_a
            ry = dx * sin_a + dy * cos_a
            rotated[idx] = [rx, ry]

        # --- 3. Canonicalize x-sign via smallest-index stable key -------
        EPS = 1e-4
        flip_x = False
        for idx in sorted(rotated.keys()):
            x = rotated[idx][0]
            if abs(x) > EPS:
                if x < 0:
                    flip_x = True
                break
        if flip_x:
            for idx in rotated:
                rotated[idx][0] = -rotated[idx][0]

        # --- 4. Canonicalize y-sign via smallest-index stable key -------
        flip_y = False
        for idx in sorted(rotated.keys()):
            y = rotated[idx][1]
            if abs(y) > EPS:
                if y < 0:
                    flip_y = True
                break
        if flip_y:
            for idx in rotated:
                rotated[idx][1] = -rotated[idx][1]

        # --- 5. Write back (centroid stays at the old centroid so the
        # molecule doesn't visually teleport; _center_coords runs after)
        for idx in self.coords:
            rx, ry = rotated[idx]
            self.coords[idx] = [rx + cx, ry + cy]

    def _optimize_layout(self, iterations=150):
        """
        Resolve chain and ring overlaps via a collision-aware spring system.

        Performance: Uses a spatial grid (cell-list) for steric repulsion
        instead of the naive O(n²) all-pairs check. Each iteration is ~O(n)
        for the repulsion phase.

        Convergence: 5 consecutive low-movement iterations (not just one)
        must occur before the loop exits early — the force field can jitter
        around a local minimum with alternating small/large steps, which the
        old single-step criterion accepted as "converged" prematurely.
        """

        # 1. Identify Rings
        rings = []
        if hasattr(self.molecule, '_rings') and self.molecule._rings:
            rings = self.molecule._rings
        else:
            rings = self.molecule.find_rings()

        ring_atoms = set()
        for r in rings:
            ring_atoms.update(r)

        # 2. Pre-calculate Ring Centers
        def get_ring_center(ring_indices):
            xs = [self.coords[idx][0] for idx in ring_indices if idx in self.coords]
            ys = [self.coords[idx][1] for idx in ring_indices if idx in self.coords]
            if not xs: return None
            return [sum(xs)/len(xs), sum(ys)/len(ys)]

        # Build set of ring-internal atom pairs to skip
        ring_pair_skip = set()
        for r in rings:
            for i in range(len(r)):
                for j in range(i + 1, len(r)):
                    pair = (min(r[i], r[j]), max(r[i], r[j]))
                    ring_pair_skip.add(pair)

        # Optimize iteratively. min_dist (1.05) is deliberately slightly
        # larger than BOND_LENGTH (1.0) so any pair of non-bonded atoms
        # is guaranteed to repel even when nominally "touching" — this is
        # the setting that prevents the visible overlap bug in mol2 files.
        damping = 0.3
        min_dist = 1.05
        grid_cell_size = min_dist * 1.5
        converge_streak = 0  # consecutive low-movement iterations

        for iteration in range(iterations):
            forces = {idx: [0.0, 0.0] for idx in self.coords}
            
            # A. Bond Spring Forces (Maintains connectivity)
            for bond in self.molecule.bonds:
                i, j = bond.begin_atom_idx, bond.end_atom_idx
                if i not in self.coords or j not in self.coords:
                    continue
                
                x1, y1 = self.coords[i]
                x2, y2 = self.coords[j]
                dx, dy = x2 - x1, y2 - y1
                dist = math.hypot(dx, dy)
                
                if dist > 0.1:
                    target = self.BOND_LENGTH
                    strength = 0.15
                    if i in ring_atoms and j in ring_atoms:
                        strength = 0.05
                    
                    mag = (dist - target) * strength
                    fx, fy = mag * dx / dist, mag * dy / dist
                    
                    forces[i][0] += fx
                    forces[i][1] += fy
                    forces[j][0] -= fx
                    forces[j][1] -= fy
                        
            # B. Rigid-Ring Repulsion
            ring_centers = [get_ring_center(r) for r in rings]
            for r1_idx in range(len(rings)):
                c1 = ring_centers[r1_idx]
                if not c1: continue
                for r2_idx in range(r1_idx + 1, len(rings)):
                    c2 = ring_centers[r2_idx]
                    if not c2: continue
                    dx, dy = c2[0] - c1[0], c2[1] - c1[1]
                    dist = math.hypot(dx, dy)
                    clash_threshold = 1.8 
                    if 0.05 < dist < clash_threshold:
                        mag = -0.3 * (clash_threshold - dist)
                        fx, fy = mag * dx / dist, mag * dy / dist
                        for atom_idx in rings[r1_idx]:
                            forces[atom_idx][0] += fx
                            forces[atom_idx][1] += fy
                        for atom_idx in rings[r2_idx]:
                            forces[atom_idx][0] -= fx
                            forces[atom_idx][1] -= fy

            # C. Spatial-grid steric repulsion (O(n) instead of O(n²))
            idxs = list(self.coords.keys())
            grid = defaultdict(list)
            for idx in idxs:
                gx = int(self.coords[idx][0] / grid_cell_size)
                gy = int(self.coords[idx][1] / grid_cell_size)
                grid[(gx, gy)].append(idx)

            for (gx, gy), cell_atoms in grid.items():
                neighbor_atoms = []
                for dx_c in (-1, 0, 1):
                    for dy_c in (-1, 0, 1):
                        neighbor_atoms.extend(grid.get((gx + dx_c, gy + dy_c), []))

                for i_pos in range(len(cell_atoms)):
                    idx1 = cell_atoms[i_pos]
                    for idx2 in neighbor_atoms:
                        if idx2 <= idx1:
                            continue
                        if self.molecule.get_bond_between(idx1, idx2):
                            continue
                        pair_key = (min(idx1, idx2), max(idx1, idx2))
                        if pair_key in ring_pair_skip:
                            continue

                        x1, y1 = self.coords[idx1]
                        x2, y2 = self.coords[idx2]
                        ddx, ddy = x2 - x1, y2 - y1
                        d = math.hypot(ddx, ddy)
                        
                        if 0.01 < d < min_dist:
                            mag = -0.5 * (min_dist - d) / d
                            fx, fy = mag * ddx / d, mag * ddy / d
                            forces[idx1][0] += fx
                            forces[idx1][1] += fy
                            forces[idx2][0] -= fx
                            forces[idx2][1] -= fy
            
            max_move = 0.0
            for idx in forces:
                fx, fy = forces[idx]
                move_x, move_y = fx * damping, fy * damping
                limit = 0.2
                move_x = max(-limit, min(limit, move_x))
                move_y = max(-limit, min(limit, move_y))
                
                self.coords[idx][0] += move_x
                self.coords[idx][1] += move_y
                max_move = max(max_move, math.hypot(move_x, move_y))
                
            if max_move < 0.002:
                converge_streak += 1
                if converge_streak >= 5:
                    break
            else:
                converge_streak = 0

        # 3. Hard de-overlap pass — guarantee that no two non-bonded
        # atoms remain within 0.95 Å of each other. Alternates
        # deoverlap sweeps with settling phases (which include
        # repulsion, so they don't reintroduce overlap) until the
        # worst-case pair is clearly above the bond length. The final
        # pass is a deoverlap with no settling after, which is safe
        # because the immediately-prior settle already balanced bond
        # forces around the separated atoms.
        for _ in range(3):
            self._hard_deoverlap_pass(ring_pair_skip, min_dist=0.95)
            self._settle_after_deoverlap(ring_atoms, ring_pair_skip, iterations=25)
        self._hard_deoverlap_pass(ring_pair_skip, min_dist=0.95)

    def _hard_deoverlap_pass(self, ring_pair_skip, min_dist=0.55):
        """Separate any non-bonded atom pair closer than ``min_dist`` by
        pushing them apart along the line joining them. Deterministic:
        pairs are processed in sorted (i, j) order so repeated runs do
        exactly the same thing."""
        idxs = sorted(self.coords.keys())
        n = len(idxs)
        bonded_pairs = set()
        for bond in self.molecule.bonds:
            a, b = bond.begin_atom_idx, bond.end_atom_idx
            bonded_pairs.add((min(a, b), max(a, b)))

        for _sweep in range(5):
            moved = False
            for ii in range(n):
                for jj in range(ii + 1, n):
                    i, j = idxs[ii], idxs[jj]
                    pair = (i, j)
                    if pair in bonded_pairs or pair in ring_pair_skip:
                        continue
                    x1, y1 = self.coords[i]
                    x2, y2 = self.coords[j]
                    dx, dy = x2 - x1, y2 - y1
                    d = math.hypot(dx, dy)
                    if d < min_dist:
                        # Push apart to min_dist + small margin
                        target = min_dist + 0.02
                        if d < 1e-6:
                            # Exactly overlapping — deterministic offset
                            # based on the pair indices (not RNG).
                            ang = (i * 37 + j * 17) % 360
                            dx = math.cos(math.radians(ang))
                            dy = math.sin(math.radians(ang))
                            d = 1.0
                        push = (target - d) * 0.5
                        ux, uy = dx / d, dy / d
                        self.coords[i] = [x1 - ux * push, y1 - uy * push]
                        self.coords[j] = [x2 + ux * push, y2 + uy * push]
                        moved = True
            if not moved:
                break

    def _settle_after_deoverlap(self, ring_atoms, ring_pair_skip, iterations=25):
        """Short spring-relaxation sweep to fix tensions introduced by the
        hard de-overlap pass. Uses bond springs AND steric repulsion,
        because bond-springs-only settling will pull previously-separated
        overlapping atoms back together via their bonded neighbors.
        """
        min_dist = 1.00  # repulsion kicks in just under bond length
        damping = 0.25
        # Precompute bonded-pairs set for fast membership
        bonded_pairs = set()
        for bond in self.molecule.bonds:
            a, b = bond.begin_atom_idx, bond.end_atom_idx
            bonded_pairs.add((min(a, b), max(a, b)))

        for _ in range(iterations):
            forces = {idx: [0.0, 0.0] for idx in self.coords}
            # Bond springs
            for bond in self.molecule.bonds:
                i, j = bond.begin_atom_idx, bond.end_atom_idx
                if i not in self.coords or j not in self.coords:
                    continue
                x1, y1 = self.coords[i]
                x2, y2 = self.coords[j]
                dx, dy = x2 - x1, y2 - y1
                dist = math.hypot(dx, dy)
                if dist > 0.1:
                    mag = (dist - self.BOND_LENGTH) * 0.25
                    fx, fy = mag * dx / dist, mag * dy / dist
                    forces[i][0] += fx; forces[i][1] += fy
                    forces[j][0] -= fx; forces[j][1] -= fy
            # Steric repulsion — O(n²) but N is small and this runs for
            # few iterations. Without this, bond springs drag previously-
            # deoverlapped atoms back on top of each other.
            idxs = list(self.coords.keys())
            for ii in range(len(idxs)):
                for jj in range(ii + 1, len(idxs)):
                    i, j = idxs[ii], idxs[jj]
                    pair = (min(i, j), max(i, j))
                    if pair in bonded_pairs or pair in ring_pair_skip:
                        continue
                    x1, y1 = self.coords[i]
                    x2, y2 = self.coords[j]
                    dx, dy = x2 - x1, y2 - y1
                    d = math.hypot(dx, dy)
                    if 0.01 < d < min_dist:
                        # Strong repulsion — scale to fully counter the
                        # bond-spring pull for this pair
                        mag = -0.8 * (min_dist - d) / d
                        fx, fy = mag * dx / d, mag * dy / d
                        forces[i][0] += fx; forces[i][1] += fy
                        forces[j][0] -= fx; forces[j][1] -= fy
            for idx, (fx, fy) in forces.items():
                mv_x = max(-0.15, min(0.15, fx * damping))
                mv_y = max(-0.15, min(0.15, fy * damping))
                self.coords[idx][0] += mv_x
                self.coords[idx][1] += mv_y

    def _place_missing_atoms(self, missing_atoms):
        """Place atoms that OASA failed to coordinate using their bonded neighbors."""
        import math
        placed_count = 0
        
        for atom_idx in missing_atoms:
            # Find neighbors that have coordinates
            neighbor_coords = []
            for bond in self.molecule.bonds:
                if bond.begin_atom_idx == atom_idx and bond.end_atom_idx in self.coords:
                    neighbor_coords.append(self.coords[bond.end_atom_idx])
                elif bond.end_atom_idx == atom_idx and bond.begin_atom_idx in self.coords:
                    neighbor_coords.append(self.coords[bond.begin_atom_idx])
            
            if neighbor_coords:
                # Place at average of neighbor positions plus small offset
                avg_x = sum(c[0] for c in neighbor_coords) / len(neighbor_coords)
                avg_y = sum(c[1] for c in neighbor_coords) / len(neighbor_coords)
                # Add small random offset to avoid exact overlap
                angle = (atom_idx % 6) * math.pi / 3  # Distribute in 6 directions
                offset = self.BOND_LENGTH * 0.8
                self.coords[atom_idx] = [avg_x + offset * math.cos(angle), 
                                         avg_y + offset * math.sin(angle)]
                placed_count += 1
            else:
                # Isolated atom - place in expanding circle
                angle = atom_idx * 0.5
                radius = 3.0 + atom_idx * 0.5
                self.coords[atom_idx] = [radius * math.cos(angle), radius * math.sin(angle)]
                placed_count += 1
        
        if placed_count > 0 and _DEBUG:
            print(f"[DEBUG] Placed {placed_count} missing atoms via fallback")

    def _center_coords(self):
        """Center the molecule coordinates around (0,0)."""
        if not self.coords:
            return
            
        min_x = min(c[0] for c in self.coords.values())
        max_x = max(c[0] for c in self.coords.values())
        min_y = min(c[1] for c in self.coords.values())
        max_y = max(c[1] for c in self.coords.values())
        
        cx = (min_x + max_x) / 2
        cy = (min_y + max_y) / 2
        
        for idx in self.coords:
            self.coords[idx][0] -= cx
            self.coords[idx][1] -= cy
