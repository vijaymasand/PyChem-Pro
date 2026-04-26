"""
Advanced Protein Visualization Module — PyMOL/Jmol-style Cartoon and Ribbon Representations

Provides professional-quality protein structure visualization.
Qt imports are done lazily inside functions to avoid import errors.

Author: Cascade AI Assistant
Date: March 24, 2026
"""

import math
import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from enum import Enum

from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom
from src.features.visualization_3d.services.cartoon_generator import CartoonGenerator

# Global generator instance (LOD auto-detected per molecule)
_cartoon_gen = CartoonGenerator()


class SecondaryStructure(Enum):
    """Secondary structure types following DSSP notation."""
    HELIX = "H"          # Alpha helix
    SHEET = "E"          # Extended strand (beta sheet)
    COIL = "C"           # Coil/loop
    TURN = "T"           # Hydrogen-bonded turn
    BEND = "S"           # Bend
    BRIDGE = "B"         # Beta bridge
    THREE_HELIX = "G"    # 3-10 helix
    PI_HELIX = "I"       # Pi helix


@dataclass
class Residue:
    """Represents a protein residue with backbone atoms."""
    index: int
    chain_id: str
    res_seq: int
    res_name: str
    ca_atom: Optional[Atom] = None
    c_atom: Optional[Atom] = None
    n_atom: Optional[Atom] = None
    o_atom: Optional[Atom] = None
    ss_type: SecondaryStructure = SecondaryStructure.COIL
    b_factor: float = 0.0


@dataclass
class Chain:
    """Represents a protein chain with its residues."""
    chain_id: str
    residues: List[Residue]


class ProteinStructure:
    """Parsed protein structure with secondary structure assignment."""
    
    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        self.chains: Dict[str, Chain] = {}
        self._parse_structure()
        self._detect_secondary_structure()
    
    def _parse_structure(self):
        """Parse molecule into chains and residues."""
        residue_data: Dict[Tuple[str, int], Dict] = {}
        
        for atom in self.molecule.atoms:
            if not hasattr(atom, 'chain_id') or not hasattr(atom, 'res_seq'):
                continue
            
            key = (atom.chain_id or 'A', atom.res_seq)
            if key not in residue_data:
                residue_data[key] = {
                    'atoms': [],
                    'res_name': getattr(atom, 'res_name', 'UNK'),
                    'b_factors': []
                }
            
            residue_data[key]['atoms'].append(atom)
            if hasattr(atom, 'b_factor'):
                residue_data[key]['b_factors'].append(atom.b_factor)
        
        for (chain_id, res_seq), data in sorted(residue_data.items()):
            residue = Residue(
                index=len(self.chains.get(chain_id, Chain(chain_id, [])).residues),
                chain_id=chain_id,
                res_seq=res_seq,
                res_name=data['res_name'],
                b_factor=np.mean(data['b_factors']) if data['b_factors'] else 0.0
            )
            
            for atom in data['atoms']:
                if hasattr(atom, 'pdb_name'):
                    name = atom.pdb_name.strip()
                    if name == 'CA':
                        residue.ca_atom = atom
                    elif name == 'C':
                        residue.c_atom = atom
                    elif name == 'N':
                        residue.n_atom = atom
                    elif name == 'O':
                        residue.o_atom = atom
            
            if residue.ca_atom is not None:
                if chain_id not in self.chains:
                    self.chains[chain_id] = Chain(chain_id, [])
                self.chains[chain_id].residues.append(residue)
    
    def _detect_secondary_structure(self):
        """Detect secondary structure.
        
        Strategy:
        1. First, try to use HELIX/SHEET records from PDB file (most accurate)
        2. If no PDB records, use DSSP hydrogen-bond energy algorithm
        """
        # Try PDB records first (like PyMOL does)
        pdb_applied = self._apply_pdb_ss_records()
        
        if pdb_applied:
            print("[DSSP] Using PDB HELIX/SHEET records for SS assignment")
            return
        
        # Fallback: DSSP hydrogen-bond energy algorithm
        print("[DSSP] No PDB SS records found, computing DSSP from coordinates...")
        self._detect_ss_dssp()
    
    def _apply_pdb_ss_records(self) -> bool:
        """Apply secondary structure from PDB HELIX/SHEET records.
        
        Returns True if PDB records were found and applied.
        """
        helix_ranges = self.molecule.properties.get('helix_ranges', [])
        sheet_ranges = self.molecule.properties.get('sheet_ranges', [])
        
        if not helix_ranges and not sheet_ranges:
            return False
        
        # Initialize all to coil
        for chain in self.chains.values():
            for r in chain.residues:
                r.ss_type = SecondaryStructure.COIL
        
        # Apply helix records: (chain_id, start_res, end_res)
        for chain_id, start_res, end_res in helix_ranges:
            if chain_id in self.chains:
                for r in self.chains[chain_id].residues:
                    if start_res <= r.res_seq <= end_res:
                        r.ss_type = SecondaryStructure.HELIX
        
        # Apply sheet records: (chain_id, start_res, end_res)
        for chain_id, start_res, end_res in sheet_ranges:
            if chain_id in self.chains:
                for r in self.chains[chain_id].residues:
                    if start_res <= r.res_seq <= end_res:
                        r.ss_type = SecondaryStructure.SHEET
        
        # Log assignment counts
        counts = {}
        for chain in self.chains.values():
            for r in chain.residues:
                ss = r.ss_type.value
                counts[ss] = counts.get(ss, 0) + 1
        print(f"[DSSP] PDB SS assignment: {counts}")
        
        return True
    
    def _detect_ss_dssp(self):
        """DSSP hydrogen-bond energy algorithm for SS detection.
        
        Uses the Kabsch-Sander electrostatic model:
        E = 0.084 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol
        A hydrogen bond exists when E < -0.5 kcal/mol.
        """
        for chain in self.chains.values():
            residues = chain.residues
            n = len(residues)
            
            if n < 4:
                continue
            
            # Initialize all to coil
            for r in residues:
                r.ss_type = SecondaryStructure.COIL
            
            # Step 1: Build H-bond energy matrix (sparse)
            hbond_energy = {}
            for i in range(n):
                for j in range(n):
                    if abs(i - j) < 2:
                        continue
                    e = self._dssp_hbond_energy(residues, i, j)
                    if e is not None and e < -0.5:
                        hbond_energy[(i, j)] = e
            
            # Step 2: Detect n-turns
            turns = {3: {}, 4: {}, 5: {}}
            for turn_n in [3, 4, 5]:
                for i in range(n - turn_n):
                    if (i, i + turn_n) in hbond_energy:
                        turns[turn_n][i] = True
            
            # Step 3: Assign helices from consecutive turns
            for i in range(n - 5):
                if turns[4].get(i) and turns[4].get(i + 1):
                    for j in range(i + 1, min(i + 5, n)):
                        if residues[j].ss_type != SecondaryStructure.SHEET:
                            residues[j].ss_type = SecondaryStructure.HELIX
            
            for i in range(n - 4):
                if turns[3].get(i) and turns[3].get(i + 1):
                    for j in range(i + 1, min(i + 4, n)):
                        if residues[j].ss_type == SecondaryStructure.COIL:
                            residues[j].ss_type = SecondaryStructure.THREE_HELIX
            
            for i in range(n - 6):
                if turns[5].get(i) and turns[5].get(i + 1):
                    for j in range(i + 1, min(i + 6, n)):
                        if residues[j].ss_type == SecondaryStructure.COIL:
                            residues[j].ss_type = SecondaryStructure.PI_HELIX
            
            # Step 4: Detect beta bridges
            bridge = [[False] * n for _ in range(n)]
            for i in range(1, n - 1):
                for j in range(i + 2, n - 1):
                    if (i - 1, j) in hbond_energy and (j, i + 1) in hbond_energy:
                        bridge[i][j] = bridge[j][i] = True
                    if (i, j) in hbond_energy and (j, i) in hbond_energy:
                        bridge[i][j] = bridge[j][i] = True
                    if j + 1 < n and (i - 1, j + 1) in hbond_energy and (j - 1, i + 1) in hbond_energy:
                        bridge[i][j] = bridge[j][i] = True
            
            for i in range(n):
                if any(bridge[i][j] for j in range(n)):
                    if residues[i].ss_type == SecondaryStructure.COIL:
                        residues[i].ss_type = SecondaryStructure.SHEET
            
            for i in range(1, n - 1):
                if (residues[i].ss_type == SecondaryStructure.COIL and
                    residues[i-1].ss_type == SecondaryStructure.SHEET and
                    residues[i+1].ss_type == SecondaryStructure.SHEET):
                    residues[i].ss_type = SecondaryStructure.SHEET
            
            # Step 5: Assign turns
            for turn_n in [4, 3, 5]:
                for i in turns[turn_n]:
                    for j in range(i, min(i + turn_n, n)):
                        if residues[j].ss_type == SecondaryStructure.COIL:
                            residues[j].ss_type = SecondaryStructure.TURN
            
            # Step 6: Clean up short segments
            self._clean_short_segments(residues, SecondaryStructure.HELIX, 3)
            self._clean_short_segments(residues, SecondaryStructure.SHEET, 2)
            
            counts = {}
            for r in residues:
                ss = r.ss_type.value
                counts[ss] = counts.get(ss, 0) + 1
            print(f"[DSSP] Chain {chain.chain_id} computed SS: {counts}")
    
    def _dssp_hbond_energy(self, residues: List['Residue'], 
                           donor_idx: int, acceptor_idx: int) -> Optional[float]:
        """Calculate DSSP hydrogen-bond energy using Kabsch-Sander formula.
        
        H-bond: CO(donor) → NH(acceptor)
        E = 0.084 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol
        
        H position estimated from the PREVIOUS residue's C atom (not the
        acceptor's own C).
        """
        donor_res = residues[donor_idx]
        acceptor_res = residues[acceptor_idx]
        
        c_atom = donor_res.c_atom
        o_atom = donor_res.o_atom
        n_atom = acceptor_res.n_atom
        
        if not all([c_atom, o_atom, n_atom]):
            return None
        if not all([c_atom.has_coords, o_atom.has_coords, n_atom.has_coords]):
            return None
        
        # Get the PREVIOUS residue's C atom for H-position estimation
        # This is the C atom from residue (acceptor_idx - 1) in the chain
        c_prev = None
        if acceptor_idx > 0:
            prev_res = residues[acceptor_idx - 1]
            if prev_res.c_atom and prev_res.c_atom.has_coords:
                c_prev = prev_res.c_atom
        
        if c_prev:
            dx = n_atom.x - c_prev.x
            dy = n_atom.y - c_prev.y
            dz = n_atom.z - c_prev.z
            d = math.sqrt(dx*dx + dy*dy + dz*dz)
            if d < 0.01:
                return None
            h_x = n_atom.x + dx / d * 1.0
            h_y = n_atom.y + dy / d * 1.0
            h_z = n_atom.z + dz / d * 1.0
        else:
            # Fallback for first residue
            if acceptor_res.ca_atom and acceptor_res.ca_atom.has_coords:
                dx = n_atom.x - acceptor_res.ca_atom.x
                dy = n_atom.y - acceptor_res.ca_atom.y
                dz = n_atom.z - acceptor_res.ca_atom.z
                d = math.sqrt(dx*dx + dy*dy + dz*dz)
                if d < 0.01:
                    return None
                h_x = n_atom.x + dx / d * 1.0
                h_y = n_atom.y + dy / d * 1.0
                h_z = n_atom.z + dz / d * 1.0
            else:
                return None
        
        r_ON = self._coord_distance(o_atom.x, o_atom.y, o_atom.z, n_atom.x, n_atom.y, n_atom.z)
        r_CH = self._coord_distance(c_atom.x, c_atom.y, c_atom.z, h_x, h_y, h_z)
        r_OH = self._coord_distance(o_atom.x, o_atom.y, o_atom.z, h_x, h_y, h_z)
        r_CN = self._coord_distance(c_atom.x, c_atom.y, c_atom.z, n_atom.x, n_atom.y, n_atom.z)
        
        if any(d < 0.01 for d in [r_ON, r_CH, r_OH, r_CN]):
            return None
        
        energy = 0.084 * 332.0 * (1.0/r_ON + 1.0/r_CH - 1.0/r_OH - 1.0/r_CN)
        return energy
    
    def _coord_distance(self, x1, y1, z1, x2, y2, z2) -> float:
        """Calculate distance between two coordinate sets."""
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    def _clean_short_segments(self, residues: List['Residue'], 
                              ss_type: SecondaryStructure, min_length: int):
        """Remove secondary structure segments shorter than min_length."""
        n = len(residues)
        i = 0
        while i < n:
            if residues[i].ss_type == ss_type:
                start = i
                while i < n and residues[i].ss_type == ss_type:
                    i += 1
                length = i - start
                if length < min_length:
                    for j in range(start, i):
                        residues[j].ss_type = SecondaryStructure.COIL
            else:
                i += 1
    
    def _atom_distance(self, atom1: Atom, atom2: Atom) -> Optional[float]:
        """Calculate distance between two atoms."""
        if not atom1 or not atom2:
            return None
        if not atom1.has_coords or not atom2.has_coords:
            return None
        
        dx = atom1.x - atom2.x
        dy = atom1.y - atom2.y
        dz = atom1.z - atom2.z
        return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    def _vector_between(self, atom1: Atom, atom2: Atom) -> Optional[Tuple[float, float, float]]:
        """Get vector from atom1 to atom2."""
        if not atom1 or not atom2:
            return None
        if not atom1.has_coords or not atom2.has_coords:
            return None
        return (atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z)
    
    def _angle_between(self, v1: Tuple[float, float, float], 
                      v2: Tuple[float, float, float]) -> Optional[float]:
        """Calculate angle between two vectors in degrees."""
        dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
        mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
        mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
        
        if mag1 == 0 or mag2 == 0:
            return None
        
        cos_angle = max(-1, min(1, dot / (mag1 * mag2)))
        return math.degrees(math.acos(cos_angle))


class SplineCalculator:
    """Catmull-Rom spline interpolation for smooth ribbon paths."""
    
    @staticmethod
    def catmull_rom_point(p0: Tuple[float, float, float],
                         p1: Tuple[float, float, float],
                         p2: Tuple[float, float, float],
                         p3: Tuple[float, float, float],
                         t: float,
                         alpha: float = 0.5) -> Tuple[float, float, float]:
        """Calculate Catmull-Rom spline point."""
        t2 = t * t
        t3 = t2 * t
        
        b0 = -alpha * t + 2 * alpha * t2 - alpha * t3
        b1 = 1 + (alpha - 3) * t2 + (2 - alpha) * t3
        b2 = alpha * t + (3 - 2 * alpha) * t2 + (alpha - 2) * t3
        b3 = -alpha * t2 + alpha * t3
        
        x = b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0]
        y = b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1]
        z = b0 * p0[2] + b1 * p1[2] + b2 * p2[2] + b3 * p3[2]
        
        return (x, y, z)
    
    @staticmethod
    def generate_spline(points: List[Tuple[float, float, float]], 
                       num_segments: int = 10) -> List[Tuple[float, float, float]]:
        """Generate smooth spline points through control points."""
        if len(points) < 2:
            return points
        
        result = []
        extended = [points[0]] + points + [points[-1]]
        
        for i in range(1, len(extended) - 2):
            p0 = extended[i - 1]
            p1 = extended[i]
            p2 = extended[i + 1]
            p3 = extended[i + 2]
            
            for j in range(num_segments):
                t = j / num_segments
                point = SplineCalculator.catmull_rom_point(p0, p1, p2, p3, t)
                result.append(point)
        
        result.append(points[-1])
        return result


def render_protein_cartoon(painter, molecule: Molecule,
                          width: int, height: int,
                          rot_x: float = 0, rot_y: float = 0, rot_z: float = 0,
                          pan_x: float = 0, pan_y: float = 0,
                          zoom: float = 1.0,
                          color_scheme: str = "secondary_structure",
                          use_ssao: bool = False,
                          use_gouraud: bool = False,
                          is_interacting: bool = False):
    """
    Render protein cartoon using the CPPCartoon-faithful mesh generator.
    
    Optimizations:
    - ALL shading is computed vectorized via NumPy before the draw loop
    - Back-face culling eliminates ~50% of triangles
    - Pre-computed QColor array minimizes per-triangle object creation
    - Thin same-color pen eliminates visible triangle edges (anti-aliasing)
    - Optional Gouraud normal smoothing blends normals across shared vertices
    - INTERACTIVE LOD: Switches to 4 subdivisions during rotation for speed.
    """
    from PySide6.QtGui import QColor, QPen, QBrush, QPolygonF, QPainter, QPainterPath, QImage
    from PySide6.QtCore import QPointF, Qt
    import time
    
    t0_total = time.time()
    
    # Antialiasing scale factor for smooth rendering
    scale_factor = 1 if is_interacting else 2  # 2x supersampling for static view
    
    # EARLY CACHE CHECK — skip ALL math if camera hasn't changed
    cache_key = (id(molecule), round(rot_x, 4), round(rot_y, 4), round(rot_z, 4),
                 round(pan_x, 2), round(pan_y, 2), round(zoom, 4),
                 width, height, is_interacting)
    
    if hasattr(render_protein_cartoon, '_img_cache') and \
       render_protein_cartoon._img_cache.get('key') == cache_key:
        painter.drawImage(0, 0, render_protein_cartoon._img_cache['img'])
        return
    
    # 1. Get the cached 3D mesh with higher quality for smooth rendering
    current_steps = 3 if is_interacting else 24  # Reduced for fast interaction
    current_profile = 4 if is_interacting else 16  # Reduced for fast interaction
    
    t0_mesh = time.time()
    vertices, triangles, colors = _cartoon_gen.get_mesh(molecule, spline_steps=current_steps, profile_detail=current_profile)
    
    if vertices is None or len(vertices) == 0:
        return
    
    num_verts = len(vertices)
    num_tri = len(triangles)
    
    # 2. Vectorized 3D -> 2D projection
    cos_x = math.cos(math.radians(rot_x))
    sin_x = math.sin(math.radians(rot_x))
    cos_y = math.cos(math.radians(rot_y))
    sin_y = math.sin(math.radians(rot_y))
    cos_z = math.cos(math.radians(rot_z))
    sin_z = math.sin(math.radians(rot_z))
    cx = width / 2.0 + pan_x
    cy = height / 2.0 + pan_y
    
    vx = vertices[:, 0] * zoom
    vy = vertices[:, 1] * zoom
    vz = vertices[:, 2] * zoom
    
    # Rotate Y then X
    rx_tmp = vx * cos_y + vz * sin_y
    rz = -vx * sin_y + vz * cos_y
    ry_tmp = vy * cos_x - rz * sin_x
    rz2 = vy * sin_x + rz * cos_x
    
    # Rotate Z
    rx = rx_tmp * cos_z - ry_tmp * sin_z
    ry = rx_tmp * sin_z + ry_tmp * cos_z
    
    sx = cx + rx       # screen x
    sy = cy - ry       # screen y (inverted)
    
    # 3. Compute face normals in rotated space (vectorized)
    t0 = triangles[:, 0]
    t1 = triangles[:, 1]
    t2 = triangles[:, 2]
    
    # Edge vectors
    e1x = rx[t1] - rx[t0]
    e1y = ry[t1] - ry[t0]
    e1z = rz2[t1] - rz2[t0]
    e2x = rx[t2] - rx[t0]
    e2y = ry[t2] - ry[t0]
    e2z = rz2[t2] - rz2[t0]
    
    # Cross product for face normals
    fnx = e1y * e2z - e1z * e2y
    fny = e1z * e2x - e1x * e2z
    fnz = e1x * e2y - e1y * e2x
    
    # Normalize face normals
    fnlen = np.sqrt(fnx*fnx + fny*fny + fnz*fnz)
    valid = fnlen > 1e-8
    fnx[valid] /= fnlen[valid]
    fny[valid] /= fnlen[valid]
    fnz[valid] /= fnlen[valid]
    
    # 4. Back-face culling: skip triangles facing away from viewer
    front_facing = fnz < 0.0
    
    # ── Gouraud Normal Smoothing ──────────────────────────────────────
    # Accumulate face normals into per-vertex normals, then shade per vertex
    # instead of per face.  This eliminates the flat-shading facet artefact
    # that becomes visible at very high zoom levels.
    if use_gouraud:
        # Accumulate face normals into vertex normals
        vnx = np.zeros(num_verts, dtype=np.float64)
        vny = np.zeros(num_verts, dtype=np.float64)
        vnz = np.zeros(num_verts, dtype=np.float64)
        
        # Only accumulate from front-facing triangles
        ff_mask = front_facing & valid
        np.add.at(vnx, t0[ff_mask], fnx[ff_mask])
        np.add.at(vnx, t1[ff_mask], fnx[ff_mask])
        np.add.at(vnx, t2[ff_mask], fnx[ff_mask])
        np.add.at(vny, t0[ff_mask], fny[ff_mask])
        np.add.at(vny, t1[ff_mask], fny[ff_mask])
        np.add.at(vny, t2[ff_mask], fny[ff_mask])
        np.add.at(vnz, t0[ff_mask], fnz[ff_mask])
        np.add.at(vnz, t1[ff_mask], fnz[ff_mask])
        np.add.at(vnz, t2[ff_mask], fnz[ff_mask])
        
        # Normalize vertex normals
        vnlen = np.sqrt(vnx*vnx + vny*vny + vnz*vnz)
        vn_valid = vnlen > 1e-8
        vnx[vn_valid] /= vnlen[vn_valid]
        vny[vn_valid] /= vnlen[vn_valid]
        vnz[vn_valid] /= vnlen[vn_valid]
        
        # Use averaged vertex normals per triangle (average of 3 vertices)
        nx = (vnx[t0] + vnx[t1] + vnx[t2]) / 3.0
        ny = (vny[t0] + vny[t1] + vny[t2]) / 3.0
        nz = (vnz[t0] + vnz[t1] + vnz[t2]) / 3.0
        # Re-normalize
        nlen2 = np.sqrt(nx*nx + ny*ny + nz*nz)
        nv = nlen2 > 1e-8
        nx[nv] /= nlen2[nv]
        ny[nv] /= nlen2[nv]
        nz[nv] /= nlen2[nv]
    else:
        # Flat shading — use face normals directly
        nx, ny, nz = fnx, fny, fnz
    
    # 5. Z-sort only front-facing triangles (back-to-front)
    tri_z = (rz2[t0] + rz2[t1] + rz2[t2]) / 3.0
    # Set back-facing triangles to -inf so they sort to the end
    tri_z[~front_facing] = -1e9
    sort_indices = np.argsort(tri_z)[::-1]
    
    # 6. Vectorized shading computation (Lambertian + Blinn-Phong specular)
    lx, ly, lz = 0.3, -0.6, -0.7
    ll = math.sqrt(lx*lx + ly*ly + lz*lz)
    lx /= ll; ly /= ll; lz /= ll
    
    # Diffuse: dot(normal, light)
    dot_nl = nx * lx + ny * ly + nz * lz
    diffuse = np.clip(dot_nl, 0.0, 1.0)
    
    # Specular (Blinn-Phong): halfway vector between light and view
    hx = lx
    hy = ly
    hz = lz - 1.0
    hl = math.sqrt(hx*hx + hy*hy + hz*hz)
    hx /= hl; hy /= hl; hz /= hl
    
    dot_nh = nx * hx + ny * hy + nz * hz
    specular = np.clip(dot_nh, 0.0, 1.0) ** 16 * 0.15
    
    ambient = 0.35
    
    if use_ssao and np.any(front_facing):
        z_min = np.min(tri_z[front_facing])
        z_max = np.max(tri_z[front_facing])
        z_range = max(1e-5, z_max - z_min)
        z_norm = (tri_z - z_min) / z_range
        
        ambient = 0.45 * (1.0 - z_norm * 0.8)
        rim = (1.0 - np.abs(nz)) ** 2.5
        shade = np.clip(ambient + diffuse * 0.65 + specular + rim * 0.35, 0.0, 1.2)
    else:
        shade = np.clip(ambient + diffuse * 0.65 + specular, 0.0, 1.2)
    
    # Per-triangle base color (from first vertex)
    tri_colors = colors[t0]  # shape (num_tri, 3)
    
    # Per-triangle final shaded color
    final_r = np.clip((tri_colors[:, 0] * shade + specular * 0.3) * 255, 0, 255).astype(np.uint8)
    final_g = np.clip((tri_colors[:, 1] * shade + specular * 0.3) * 255, 0, 255).astype(np.uint8)
    final_b = np.clip((tri_colors[:, 2] * shade + specular * 0.3) * 255, 0, 255).astype(np.uint8)
    
    t0_draw = time.time()
    
    # 7. CHUNKED Vectorized Rasterizer (memory-safe for 8GB RAM)
    active_indices = sort_indices[front_facing[sort_indices]]
    if len(active_indices) == 0:
        return

    N = len(active_indices)
    tri_act = triangles[active_indices]
    
    # Screen coords per triangle vertex (no scaling here)
    t_ax = sx[tri_act[:, 0]]; t_ay = sy[tri_act[:, 0]]
    t_bx = sx[tri_act[:, 1]]; t_by = sy[tri_act[:, 1]]
    t_cx = sx[tri_act[:, 2]]; t_cy = sy[tri_act[:, 2]]
    
    tr = final_r[active_indices]
    tg = final_g[active_indices]
    tb = final_b[active_indices]

    # RGBA pixel buffer with higher resolution for antialiasing
    buf_h = height * scale_factor
    buf_w = width * scale_factor
    buf = np.zeros((buf_h, buf_w, 4), dtype=np.uint8)
    
    MAX_BB = 20  # Increased for higher resolution
    CHUNK = 3000  # Process 3000 triangles at a time (~4MB per chunk)
    
    gy, gx = np.meshgrid(np.arange(MAX_BB, dtype=np.float32),
                         np.arange(MAX_BB, dtype=np.float32), indexing='ij')
    
    for c_start in range(0, N, CHUNK):
        c_end = min(c_start + CHUNK, N)
        sl = slice(c_start, c_end)
        Nc = c_end - c_start
        
        c_ax = t_ax[sl]; c_ay = t_ay[sl]
        c_bx = t_bx[sl]; c_by = t_by[sl]
        c_cx = t_cx[sl]; c_cy = t_cy[sl]
        
        # Bounding boxes (scaled for higher resolution)
        mn_x = np.floor(np.minimum(np.minimum(c_ax, c_bx), c_cx) * scale_factor).astype(np.int32)
        mx_x = np.ceil(np.maximum(np.maximum(c_ax, c_bx), c_cx) * scale_factor).astype(np.int32)
        mn_y = np.floor(np.minimum(np.minimum(c_ay, c_by), c_cy) * scale_factor).astype(np.int32)
        mx_y = np.ceil(np.maximum(np.maximum(c_ay, c_by), c_cy) * scale_factor).astype(np.int32)
        
        bb_w = mx_x - mn_x + 1
        bb_h = mx_y - mn_y + 1
        small = (bb_w <= MAX_BB) & (bb_h <= MAX_BB)
        si = np.where(small)[0]
        
        if len(si) > 0:
            PX = mn_x[si, np.newaxis, np.newaxis] + gx[np.newaxis]
            PY = mn_y[si, np.newaxis, np.newaxis] + gy[np.newaxis]
            
            AX = c_ax[si, np.newaxis, np.newaxis] * scale_factor
            AY = c_ay[si, np.newaxis, np.newaxis] * scale_factor
            BX = c_bx[si, np.newaxis, np.newaxis] * scale_factor
            BY = c_by[si, np.newaxis, np.newaxis] * scale_factor
            CX = c_cx[si, np.newaxis, np.newaxis] * scale_factor
            CY = c_cy[si, np.newaxis, np.newaxis] * scale_factor
            
            e0 = (BX-AX)*(PY-AY) - (BY-AY)*(PX-AX)
            e1 = (CX-BX)*(PY-BY) - (CY-BY)*(PX-BX)
            e2 = (AX-CX)*(PY-CY) - (AY-CY)*(PX-CX)
            
            inside = ((e0>=0)&(e1>=0)&(e2>=0))|((e0<=0)&(e1<=0)&(e2<=0))
            inb = (PX>=0)&(PX<buf_w)&(PY>=0)&(PY<buf_h)
            mask = inside & inb
            
            fpx = PX[mask].astype(np.intp)
            fpy = PY[mask].astype(np.intp)
            tidx = np.broadcast_to(np.arange(len(si))[:,np.newaxis,np.newaxis],
                                   (len(si),MAX_BB,MAX_BB))[mask]
            
            cr = tr[sl][si]; cg = tg[sl][si]; cb = tb[sl][si]
            buf[fpy, fpx, 0] = cr[tidx]
            buf[fpy, fpx, 1] = cg[tidx]
            buf[fpy, fpx, 2] = cb[tidx]
            buf[fpy, fpx, 3] = 255
        
        # Large triangles fallback (scaled)
        for li in np.where(~small)[0]:
            _mn_x = max(0, int(mn_x[li])); _mx_x = min(buf_w-1, int(mx_x[li]))
            _mn_y = max(0, int(mn_y[li])); _mx_y = min(buf_h-1, int(mx_y[li]))
            if _mn_x > _mx_x or _mn_y > _mx_y: continue
            px_r = np.arange(_mn_x, _mx_x+1, dtype=np.float32)
            py_r = np.arange(_mn_y, _mx_y+1, dtype=np.float32)
            LPX, LPY = np.meshgrid(px_r, py_r)
            le0=( float(c_bx[li])*scale_factor-float(c_ax[li])*scale_factor)*(LPY-float(c_ay[li])*scale_factor)-(float(c_by[li])*scale_factor-float(c_ay[li])*scale_factor)*(LPX-float(c_ax[li])*scale_factor)
            le1=( float(c_cx[li])*scale_factor-float(c_bx[li])*scale_factor)*(LPY-float(c_by[li])*scale_factor)-(float(c_cy[li])*scale_factor-float(c_by[li])*scale_factor)*(LPX-float(c_bx[li])*scale_factor)
            le2=( float(c_ax[li])*scale_factor-float(c_cx[li])*scale_factor)*(LPY-float(c_cy[li])*scale_factor)-(float(c_ay[li])*scale_factor-float(c_cy[li])*scale_factor)*(LPX-float(c_cx[li])*scale_factor)
            lm=((le0>=0)&(le1>=0)&(le2>=0))|((le0<=0)&(le1<=0)&(le2<=0))
            ys,xs=np.where(lm)
            buf[ys+_mn_y,xs+_mn_x,0]=tr[c_start+li]
            buf[ys+_mn_y,xs+_mn_x,1]=tg[c_start+li]
            buf[ys+_mn_y,xs+_mn_x,2]=tb[c_start+li]
            buf[ys+_mn_y,xs+_mn_x,3]=255

    # Downsample for antialiasing (2x supersampling)
    if scale_factor > 1:
        # Simple averaging downsampling for antialiasing
        buf_small = buf.reshape(height, scale_factor, width, scale_factor, 4).mean(axis=(1, 3)).astype(np.uint8)
        img = QImage(buf_small.data, width, height, width*4, QImage.Format.Format_RGBA8888)
    else:
        img = QImage(buf.data, buf_w, buf_h, buf_w*4, QImage.Format.Format_RGBA8888)
    img = img.copy()
    
    if not hasattr(render_protein_cartoon, '_img_cache'):
        render_protein_cartoon._img_cache = {}
    render_protein_cartoon._img_cache = {'key': cache_key, 'img': img}
    
    painter.drawImage(0, 0, img)
    
    draw_ms = (time.time()-t0_draw)*1000
    total_ms = (time.time()-t0_total)*1000




def render_protein_ribbon(painter, molecule: Molecule,
                         width: int, height: int,
                         rot_x: float = 0, rot_y: float = 0, rot_z: float = 0,
                         pan_x: float = 0, pan_y: float = 0,
                         zoom: float = 1.0,
                         color_scheme: str = "secondary_structure"):
    """
    Render smooth ribbon representation.
    """
    from PySide6.QtGui import QColor, QPen
    from PySide6.QtCore import Qt
    from src.shared.ui.theme import COLORS
    
    # Color definitions - use theme colors if available, fall back to defaults
    def get_ss_color(ss_type, default_rgb):
        """Get color for secondary structure type from theme or use default."""
        color_map = {
            SecondaryStructure.HELIX: 'ss_helix',
            SecondaryStructure.SHEET: 'ss_sheet',
            SecondaryStructure.COIL: 'ss_coil',
            SecondaryStructure.TURN: 'ss_turn',
        }
        
        theme_key = color_map.get(ss_type)
        if theme_key and theme_key in COLORS:
            hex_color = COLORS[theme_key]
            hex_color = hex_color.lstrip('#')
            r = int(hex_color[0:2], 16)
            g = int(hex_color[2:4], 16)
            b = int(hex_color[4:6], 16)
            return QColor(r, g, b)
        else:
            return QColor(*default_rgb)
    
    SS_COLORS = {
        SecondaryStructure.HELIX: get_ss_color(SecondaryStructure.HELIX, (220, 50, 50)),
        SecondaryStructure.SHEET: get_ss_color(SecondaryStructure.SHEET, (50, 150, 220)),
        SecondaryStructure.COIL: get_ss_color(SecondaryStructure.COIL, (180, 180, 180)),
        SecondaryStructure.TURN: get_ss_color(SecondaryStructure.TURN, (0, 212, 170)),
        SecondaryStructure.BEND: QColor(200, 200, 100),
        SecondaryStructure.BRIDGE: QColor(128, 0, 128),
        SecondaryStructure.THREE_HELIX: QColor(255, 100, 100),
        SecondaryStructure.PI_HELIX: QColor(200, 50, 50),
    }
    
    def rainbow_color(position, total):
        ratio = position / total if total > 0 else 0
        hue = int(240 * (1 - ratio))
        return QColor.fromHsv(hue, 255, 255)
    
    def bfactor_color(b_factor, min_bf=0, max_bf=100):
        if max_bf <= min_bf:
            return QColor(128, 128, 128)
        t = (b_factor - min_bf) / (max_bf - min_bf)
        t = max(0, min(1, t))
        r = int(255 * t)
        b = int(255 * (1 - t))
        return QColor(r, 0, b)
    
    protein = ProteinStructure(molecule)
    
    cos_x, sin_x = math.cos(math.radians(rot_x)), math.sin(math.radians(rot_x))
    cos_y, sin_y = math.cos(math.radians(rot_y)), math.sin(math.radians(rot_y))
    cos_z, sin_z = math.cos(math.radians(rot_z)), math.sin(math.radians(rot_z))
    cx, cy = width / 2 + pan_x, height / 2 + pan_y
    
    def project_point(x, y, z):
        x *= zoom
        y *= zoom
        z *= zoom
        x1 = x * cos_y + z * sin_y
        z1 = -x * sin_y + z * cos_y
        y1 = y * cos_x - z1 * sin_x
        z2 = y * sin_x + z1 * cos_x
        
        x2 = x1 * cos_z - y1 * sin_z
        y2 = x1 * sin_z + y1 * cos_z
        return (cx + x2, cy - y2, z2)
    
    CHAIN_COLORS = [
        QColor(0, 100, 255), QColor(255, 0, 0), QColor(0, 200, 0),
        QColor(255, 165, 0), QColor(128, 0, 128), QColor(0, 200, 200),
        QColor(255, 255, 0), QColor(255, 0, 255),
    ]
    
    chain_idx = 0
    for chain_id, chain in protein.chains.items():
        # Extract CA points and colors
        ca_points = []
        colors = []
        
        for residue in chain.residues:
            if residue.ca_atom and residue.ca_atom.has_coords:
                ca_points.append((residue.ca_atom.x, residue.ca_atom.y, residue.ca_atom.z))
                
                if color_scheme == "secondary_structure":
                    colors.append(SS_COLORS.get(residue.ss_type, QColor(180, 180, 180)))
                elif color_scheme == "rainbow":
                    colors.append(rainbow_color(residue.index, len(chain.residues)))
                elif color_scheme == "bfactor":
                    colors.append(bfactor_color(residue.b_factor))
                else:
                    colors.append(CHAIN_COLORS[chain_idx % len(CHAIN_COLORS)])
        
        if len(ca_points) < 2:
            continue
        
        # Generate spline
        smooth_points = SplineCalculator.generate_spline(ca_points, num_segments=20)
        projected = [project_point(*p) for p in smooth_points]
        
        # Draw ribbon with color transitions
        ribbon_width = 6.0
        segment_length = len(projected) // len(colors) if colors else 1
        
        for i in range(len(projected) - 1):
            color_idx = min(i // segment_length, len(colors) - 1)
            color = colors[color_idx] if colors else QColor(180, 180, 180)
            
            pen = QPen(color, ribbon_width)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            
            painter.setPen(pen)
            painter.drawLine(projected[i][0], projected[i][1], 
                           projected[i + 1][0], projected[i + 1][1])
        
        chain_idx += 1
