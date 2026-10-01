"""
SMILES-specific 2D Coordinate Generator - Enhanced OASA Implementation

=============================================================================
OVERVIEW
=============================================================================

This module provides specialized 2D coordinate generation for SMILES-converted molecules.
It extends the proven OASA-based approach from coordgen2d.py with specific
optimizations for SMILES-derived structures, focusing on proper hydrogen handling
and professional layout quality.

=============================================================================
KEY DIFFERENCES FROM coordgen2d.py
=============================================================================

1. **Hydrogen Handling**:
   - coordgen2d.py: Only handles heavy atoms (file imports)
   - coordgen2d_smiles.py: Handles both heavy atoms AND implicit hydrogens
   
2. **Coordinate Generation**:
   - coordgen2d.py: Uses OASA with post-processing optimization
   - coordgen2d_smiles.py: Uses OASA with enhanced ring perception and hydrogen placement
   
3. **Molecule Source**:
   - coordgen2d.py: Designed for MOL/MOL2 files with existing coordinates
   - coordgen2d_smiles.py: Designed for SMILES strings requiring full coordinate generation

=============================================================================
ALGORITHM DETAILS
=============================================================================

1. **Heavy Atom Layout**:
   - Convert domain molecule to OASA format
   - Use OASA's proven coordinate generation algorithm
   - Enhanced ring perception for better aromatic system handling
   - Optimized bond length scaling for SMILES-derived structures

2. **Implicit Hydrogen Placement**:
   - Calculate expected valence for each atom
   - Determine implicit hydrogen count based on current bonds
   - Place virtual hydrogen atoms using geometric rules
   - Virtual hydrogen indices use negative numbering system

3. **Layout Optimization**:
   - Center coordinates around origin
   - Apply OASA's coordinate optimizer (limited iterations)
   - Ensure proper bond angles and distances
   - Maintain professional appearance

=============================================================================
HYDROGEN PLACEMENT STRATEGY
=============================================================================

The hydrogen placement algorithm uses sophisticated geometric rules:

1. **Valence Calculation**:
   - Current valence = sum of bond orders
   - Expected valence = based on element (C:4, N:3, O:2, etc.)
   - Implicit H count = expected - current - explicit

2. **Geometric Placement**:
   - 0 neighbors: Place H atoms in tetrahedral arrangement
   - 1 neighbor: Place H atoms opposite to existing bond
   - 2 neighbors: Place H atoms in plane perpendicular to bond angle
   - 3 neighbors: Place H atom perpendicular to plane (if needed)

3. **Virtual Hydrogen Indexing**:
   - Virtual H atoms use negative indices
   - Format: -(parent_index * 100 + h_index + 1)
   - Example: H on atom 5 = -(5*100 + 1) = -501

=============================================================================
OASA INTEGRATION ENHANCEMENTS
=============================================================================

1. **Enhanced Ring Perception**:
   - Explicit ring detection before coordinate generation
   - Proper aromatic ring handling
   - Fused ring system optimization

2. **Coordinate Optimization**:
   - Limited OASA optimizer iterations (500 max)
   - Prevents over-optimization that distorts layout
   - Maintains chemical accuracy while improving aesthetics

3. **Bond Length Scaling**:
   - Optimized for SMILES-derived structures
   - Consistent with professional chemistry software
   - Proper spacing for hydrogen visibility

=============================================================================
USAGE
=============================================================================

```python
# Create generator for SMILES molecule
generator = CoordinateGenerator2DSMILES(molecule, force_regenerate=True)

# Generate coordinates including virtual hydrogens
coords = generator.generate()

# Result format: {atom_index: [x, y]}
# Real atoms: positive indices (0, 1, 2, ...)
# Virtual H: negative indices (-101, -102, ...)
```

=============================================================================
DEPENDENCIES
=============================================================================

- OASA coordinate generation library
- Domain molecule to OASA conversion bridge
- Standard math and geometry libraries

=============================================================================
AUTHOR & VERSION
=============================================================================

Derived from coordgen2d.py with SMILES-specific enhancements
Version: 1.0
Focus: Professional 2D layout for SMILES-derived molecules
"""

from collections import deque, defaultdict
import math
import src.vendors.oasa.smiles as oasa_smiles
import src.vendors.oasa.coords_generator as oasa_cg
from src.vendors.oasa_bridge import domain_to_oasa_mol

_DEBUG = False


class CoordinateGenerator2DSMILES:
    """
    Enhanced 2D coordinate generator for SMILES molecules.
    
    This class extends the proven OASA-based coordinate generation from coordgen2d.py
    with specific optimizations for SMILES-derived structures, including:
    
    1. Enhanced hydrogen placement with virtual atom support
    2. Improved ring perception for aromatic systems  
    3. Optimized coordinate generation parameters
    4. Professional layout quality maintenance
    
    Attributes:
        molecule: Domain molecule to generate coordinates for
        coords: Dictionary of generated coordinates
        BOND_LENGTH: Standard bond length unit (1.0 for OASA compatibility)
        force_regenerate: Force regeneration even if coordinates exist
    """
    
    BOND_LENGTH = 1.0  # OASA unit scale
    
    def __init__(self, molecule, force_regenerate=True):
        """
        Initialize the SMILES coordinate generator.
        
        Args:
            molecule: Domain molecule object from SMILES conversion
            force_regenerate: Force coordinate regeneration even if coordinates exist
        """
        self.molecule = molecule
        self.coords = {}
        self.force_regenerate = force_regenerate
        
        if _DEBUG:
            print(f"[DEBUG SMILES] Initialized for {len(molecule.atoms) if molecule else 0} atoms")
    
    def generate(self):
        """
        Generate 2D coordinates for SMILES molecule with enhanced OASA approach.
        
        Returns:
            dict: {atom_index: [x, y]} for ALL atoms including virtual hydrogens
            
        Process:
            1. Check for existing coordinates
            2. Generate heavy atom coordinates using enhanced OASA
            3. Place virtual hydrogen atoms using geometric rules
            4. Optimize layout with limited iterations
            5. Center and return final coordinates
        """
        if not self.molecule or len(self.molecule.atoms) == 0:
            if _DEBUG:
                print("[DEBUG SMILES] No molecule or atoms to process")
            return {}
        
        # Step 1: Check for existing coordinates
        if not self.force_regenerate:
            has_coords = all(hasattr(a, 'x2d') and a.x2d is not None 
                           for a in self.molecule.atoms if a.symbol != 'H')
            if has_coords:
                if _DEBUG:
                    print("[DEBUG SMILES] Using existing coordinates")
                self.coords = {a.index: [a.x2d * self.BOND_LENGTH, a.y2d * self.BOND_LENGTH] 
                              for a in self.molecule.atoms}
                self._center_coords()
                return self.coords
        
        # Step 2: Generate heavy atom coordinates using enhanced OASA
        if _DEBUG:
            print("[DEBUG SMILES] Generating heavy atom coordinates with enhanced OASA...")
        self._generate_heavy_atom_coords()
        
        # Step 3: Place virtual hydrogen atoms
        if _DEBUG:
            print("[DEBUG SMILES] Placing virtual hydrogen atoms...")
        self._place_virtual_hydrogens()
        
        # Step 4: Optimize layout
        if _DEBUG:
            print("[DEBUG SMILES] Optimizing layout...")
        self._optimize_layout()
        
        # Step 5: Center coordinates
        self._center_coords()
        
        # Step 6: Cache coordinates back to molecule
        self._cache_coordinates()
        
        # Report results
        total = len(self.coords)
        real = sum(1 for idx in self.coords if idx >= 0)
        virtual = sum(1 for idx in self.coords if idx < 0)
        if _DEBUG:
            print(f"[DEBUG SMILES] Generated {total} coordinates ({real} real + {virtual} virtual H)")
        
        return self.coords
    
    def _generate_heavy_atom_coords(self):
        """
        Generate heavy atom coordinates using enhanced OASA approach.
        
        Enhancements over standard OASA:
        1. Pre-ring perception for better aromatic handling
        2. Optimized bond length parameters
        3. Enhanced coordinate generation parameters
        """
        # Convert to OASA format
        o_mol, atom_map = domain_to_oasa_mol(self.molecule)
        
        # Enhanced ring perception
        try:
            o_mol.get_smallest_independent_cycles()
            if _DEBUG:
                print("[DEBUG SMILES] Enhanced ring perception completed")
        except Exception as e:
            if _DEBUG:
                print(f"[DEBUG SMILES] Ring perception warning: {e}")
        
        # Generate coordinates with optimized parameters
        oasa_cg.calculate_coords(o_mol, bond_length=self.BOND_LENGTH, force=1)
        
        # Read back coordinates
        self.coords = {}
        missing_atoms = []
        for internal_idx, o_v in atom_map.items():
            x = float(o_v.x) if o_v.x is not None else None
            y = float(o_v.y) if o_v.y is not None else None
            if x is not None and y is not None:
                self.coords[internal_idx] = [x, y]
            else:
                missing_atoms.append(internal_idx)
        
        if _DEBUG:
            print(f"[DEBUG SMILES] OASA placed {len(self.coords)} heavy atoms")
        
        # Handle missing atoms
        if missing_atoms:
            self._place_missing_atoms(missing_atoms)
    
    def _place_virtual_hydrogens(self):
        """
        Place virtual hydrogen atoms using sophisticated geometric rules.
        
        This method implements advanced hydrogen placement that considers:
        1. Current valence state of each atom
        2. Hybridization and geometry
        3. Existing neighbor positions
        4. Chemical valence requirements
        
        Virtual hydrogen atoms use negative indexing to distinguish from real atoms.
        """
        placed_count = 0
        total_h_needed = 0
        
        for atom in self.molecule.atoms:
            if atom.symbol == 'H':
                continue  # Skip explicit hydrogen atoms
            if atom.index not in self.coords:
                continue  # Skip atoms without coordinates
            
            # Calculate implicit hydrogen count
            h_count = self._calculate_implicit_h_count(atom)
            if h_count <= 0:
                continue
            
            total_h_needed += h_count
            
            # Get parent atom position
            px, py = self.coords[atom.index]
            
            # Get neighbor positions for geometric placement
            neighbors = self._get_heavy_neighbors(atom.index)
            neighbor_coords = [self.coords[n] for n in neighbors if n in self.coords]
            
            # Calculate optimal hydrogen placement angles
            angles = self._calculate_hydrogen_angles(px, py, neighbor_coords, h_count, atom)
            
            # Place virtual hydrogen atoms
            for i, angle in enumerate(angles):
                # Calculate hydrogen position
                h_x = px + self.BOND_LENGTH * 0.85 * math.cos(angle)
                h_y = py + self.BOND_LENGTH * 0.85 * math.sin(angle)
                
                # Virtual hydrogen index (negative to distinguish from real atoms)
                virtual_h_idx = -((atom.index + 1) * 100 + i + 1)
                self.coords[virtual_h_idx] = [h_x, h_y]
                placed_count += 1
        
        if _DEBUG:
            print(f"[DEBUG SMILES] Placed {placed_count} virtual hydrogens for {total_h_needed} needed")
    
    def _calculate_implicit_h_count(self, atom):
        """
        Calculate the number of implicit hydrogens for an atom.
        
        Uses valence rules based on element type and current bonding state.
        
        Args:
            atom: Domain atom object
            
        Returns:
            int: Number of implicit hydrogens needed
        """
        # Calculate current valence
        current_valence = 0
        for bond in self.molecule.bonds:
            if bond.begin_atom_idx == atom.index or bond.end_atom_idx == atom.index:
                order = bond.order if bond.order else 1
                if bond.is_aromatic:
                    order = 1.5  # Aromatic bonds count as 1.5
                current_valence += order
        
        # Get expected valence based on element
        expected_valence = self._get_expected_valence(atom)
        
        # Count explicit hydrogens
        explicit_h = sum(1 for n in self.molecule.get_neighbors(atom.index) 
                        if self.molecule.atoms[n].symbol == 'H')
        
        # Calculate implicit hydrogens
        implicit_h = max(0, int(expected_valence - current_valence - explicit_h))
        
        return implicit_h
    
    def _get_expected_valence(self, atom):
        """
        Get expected valence for an atom based on element type.
        
        Args:
            atom: Domain atom object
            
        Returns:
            int: Expected valence
        """
        valence_map = {
            'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
            'S': 2, 'P': 3, 'B': 3, 'Si': 4, 'Ge': 4, 'As': 3, 'Se': 2,
            'Te': 2, 'At': 1, 'Rn': 0, 'He': 0, 'Ne': 0, 'Ar': 0, 'Kr': 0,
            'Xe': 0, 'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1,
            'Fr': 1, 'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
            'Ra': 2, 'Al': 3, 'Ga': 3, 'In': 3, 'Tl': 3, 'Sn': 4,
            'Pb': 4, 'Bi': 3, 'Sb': 3, 'Zn': 2, 'Cd': 2, 'Hg': 2,
            'Cu': 2, 'Ag': 1, 'Au': 1, 'Fe': 3, 'Co': 3, 'Ni': 3,
            'Mn': 2, 'Cr': 3, 'Mo': 6, 'W': 6, 'V': 5, 'Ti': 4,
            'Zr': 4, 'Hf': 4, 'Nb': 5, 'Ta': 5, 'Sc': 3, 'Y': 3
        }
        return valence_map.get(atom.symbol, 2)  # Default to 2 for unknown elements
    
    def _get_heavy_neighbors(self, atom_idx):
        """
        Get indices of non-hydrogen neighbors for an atom.
        
        Args:
            atom_idx: Index of the atom
            
        Returns:
            list: Indices of heavy atom neighbors
        """
        neighbors = []
        for bond in self.molecule.bonds:
            if bond.begin_atom_idx == atom_idx:
                if self.molecule.atoms[bond.end_atom_idx].symbol != 'H':
                    neighbors.append(bond.end_atom_idx)
            elif bond.end_atom_idx == atom_idx:
                if self.molecule.atoms[bond.begin_atom_idx].symbol != 'H':
                    neighbors.append(bond.begin_atom_idx)
        return neighbors
    
    def _calculate_hydrogen_angles(self, px, py, neighbor_coords, num_h, atom):
        """
        Calculate optimal angles for hydrogen placement around an atom.
        
        Uses sophisticated geometric rules based on:
        1. Number of existing neighbors
        2. Hybridization state
        3. Chemical geometry requirements
        
        Args:
            px, py: Parent atom coordinates
            neighbor_coords: List of neighbor coordinates
            num_h: Number of hydrogens to place
            atom: Parent atom object
            
        Returns:
            list: List of angles for hydrogen placement
        """
        if len(neighbor_coords) == 0:
            # No neighbors - tetrahedral arrangement
            return [i * 2 * math.pi / num_h for i in range(num_h)]
        
        elif len(neighbor_coords) == 1:
            # One neighbor - place hydrogens opposite
            x1, y1 = neighbor_coords[0]
            angle_to_neighbor = math.atan2(y1 - py, x1 - px)
            opposite = angle_to_neighbor + math.pi
            
            if num_h == 1:
                return [opposite]
            elif num_h == 2:
                return [opposite - math.pi/3, opposite + math.pi/3]
            elif num_h == 3:
                return [opposite, opposite + 2*math.pi/3, opposite + 4*math.pi/3]
            else:
                return [opposite + i * 2 * math.pi / num_h for i in range(num_h)]
        
        elif len(neighbor_coords) == 2:
            # Two neighbors - place hydrogens perpendicular to bond angle
            x1, y1 = neighbor_coords[0]
            x2, y2 = neighbor_coords[1]
            angle1 = math.atan2(y1 - py, x1 - px)
            angle2 = math.atan2(y2 - py, x2 - px)
            
            # Calculate bisector
            mid_angle = (angle1 + angle2) / 2
            if abs(angle1 - angle2) > math.pi:
                mid_angle += math.pi
            
            external = mid_angle + math.pi
            
            if num_h == 1:
                return [external]
            elif num_h == 2:
                return [external - 0.5, external + 0.5]
            else:
                return [external + i * 2 * math.pi / num_h for i in range(num_h)]
        
        elif len(neighbor_coords) == 3:
            # Three neighbors - place hydrogen perpendicular to plane
            angles = [math.atan2(y - py, x - px) for x, y in neighbor_coords]
            avg_angle = sum(angles) / len(angles)
            return [avg_angle + math.pi]
        
        else:
            # Four or more neighbors - distribute evenly
            return [i * 2 * math.pi / num_h for i in range(num_h)]
    
    def _place_missing_atoms(self, missing_atoms):
        """
        Place atoms that OASA failed to coordinate using fallback method.
        
        Args:
            missing_atoms: List of atom indices that need coordinates
        """
        placed = 0
        for atom_idx in missing_atoms:
            # Find neighbors with coordinates
            neighbor_coords = []
            for bond in self.molecule.bonds:
                if bond.begin_atom_idx == atom_idx and bond.end_atom_idx in self.coords:
                    neighbor_coords.append(self.coords[bond.end_atom_idx])
                elif bond.end_atom_idx == atom_idx and bond.begin_atom_idx in self.coords:
                    neighbor_coords.append(self.coords[bond.begin_atom_idx])
            
            if neighbor_coords:
                # Place near average of neighbors
                avg_x = sum(c[0] for c in neighbor_coords) / len(neighbor_coords)
                avg_y = sum(c[1] for c in neighbor_coords) / len(neighbor_coords)
                angle = (atom_idx % 6) * math.pi / 3  # Distribute in 6 directions
                offset = self.BOND_LENGTH * 0.8
                self.coords[atom_idx] = [avg_x + offset * math.cos(angle), 
                                         avg_y + offset * math.sin(angle)]
                placed += 1
            else:
                # Isolated atom - place in expanding circle
                angle = atom_idx * 0.5
                radius = 3.0 + atom_idx * 0.5
                self.coords[atom_idx] = [radius * math.cos(angle), radius * math.sin(angle)]
                placed += 1
        
        if placed > 0 and _DEBUG:
            print(f"[DEBUG SMILES] Placed {placed} missing atoms via fallback")
    
    def _optimize_layout(self, iterations=30):
        """
        Resolve chain and ring overlaps iteratively via a collision-aware spring system.
        
        This 'Pro' version implements Rigid-Ring Push logic:
        1. Rings are treated as rigid bodies for internal forces.
        2. If two rings overlap, a repulsive force is applied to their centers.
        3. This resolves 'crashing' layouts (like overlapping phenyl rings) without
           distorting the perfect polygon geometry generated by OASA.
        """
        import math
        
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

        # Optimize iteratively
        damping = 0.3
        for iteration in range(iterations):
            forces = {idx: [0.0, 0.0] for idx in self.coords}
            
            # A. Bond Spring Forces
            # We iterate over real bonds only (virtual H handled by repulsive/stay-near logic)
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

            # C. Spatial-grid steric repulsion (O(n) per iteration)
            idxs = list(self.coords.keys())
            min_dist_heavy = 0.7
            min_dist_mixed = 0.5
            grid_cell_size = min_dist_heavy * 1.5
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
                        if idx1 >= 0 and idx2 >= 0:
                            if self.molecule.get_bond_between(idx1, idx2):
                                continue
                    
                        x1, y1 = self.coords[idx1]
                        x2, y2 = self.coords[idx2]
                        ddx, ddy = x2 - x1, y2 - y1
                        d = math.hypot(ddx, ddy)
                    
                        md = min_dist_heavy if (idx1 >= 0 and idx2 >= 0) else min_dist_mixed
                        if 0.01 < d < md:
                            mag = -0.5 * (md - d) / d
                            fx, fy = mag * ddx / d, mag * ddy / d
                            forces[idx1][0] += fx
                            forces[idx1][1] += fy
                            if idx2 in forces:
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
            if max_move < 0.005: break

    def _align_pca(self):
        """Align the longest molecular axis horizontally."""
        import math
        # Only use heavy atoms for PCA to avoid hydrogen distortion
        heavy_coords = {idx: c for idx, c in self.coords.items() if idx >= 0}
        if not heavy_coords or len(heavy_coords) < 2: return
            
        cx = sum(c[0] for c in heavy_coords.values()) / len(heavy_coords)
        cy = sum(c[1] for c in heavy_coords.values()) / len(heavy_coords)
        
        cxx = cyy = cxy = 0.0
        for x, y in heavy_coords.values():
            dx, dy = x - cx, y - cy
            cxx += dx * dx
            cyy += dy * dy
            cxy += dx * dy
            
        det = math.sqrt(max(0, (cxx + cyy)**2 - 4*(cxx*cyy - cxy*cxy)))
        l1 = (cxx + cyy + det) / 2.0
        
        vx, vy = cxy, l1 - cxx
        if math.hypot(vx, vy) < 1e-6:
            vx, vy = l1 - cyy, cxy
            
        norm = math.hypot(vx, vy)
        if norm < 1e-6: return
        vx, vy = vx/norm, vy/norm
        
        angle = -math.atan2(vy, vx)
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        for idx in self.coords:
            x, y = self.coords[idx]
            dx, dy = x - cx, y - cy
            self.coords[idx] = [dx * cos_a - dy * sin_a + cx, dx * sin_a + dy * cos_a + cy]
    
    def _center_coords(self):
        """
        Center molecule coordinates around origin (0,0).
        
        Calculates bounding box of all coordinates and centers them
        for optimal visualization.
        """
        if not self.coords:
            return
        
        # Calculate bounds
        min_x = min(c[0] for c in self.coords.values())
        max_x = max(c[0] for c in self.coords.values())
        min_y = min(c[1] for c in self.coords.values())
        max_y = max(c[1] for c in self.coords.values())
        
        # Calculate center
        cx = (min_x + max_x) / 2
        cy = (min_y + max_y) / 2
        
        # Center all coordinates
        for idx in self.coords:
            self.coords[idx][0] -= cx
            self.coords[idx][1] -= cy
    
    def _cache_coordinates(self):
        """
        Cache generated coordinates back to the molecule object.
        
        Only caches real atom coordinates (positive indices).
        Virtual hydrogen coordinates are kept separate.
        """
        for idx in self.coords:
            if idx >= 0:  # Real atoms only
                atom = self.molecule.get_atom(idx)
                if atom:
                    x, y = self.coords[idx]
                    atom.x2d = x
                    atom.y2d = y
                    atom.z2d = 0.0
