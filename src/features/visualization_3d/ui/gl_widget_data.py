import numpy as np
import time
import math
from src.shared.qt_compat import QOpenGLBuffer
from src.features.visualization_3d.ui.gl_widget_helpers import _element_color_float, _display_radius

class GLDataMixin:
    def set_atom_colors(self, colors_dict):
        """Set custom colors for specific atoms and rebuild GPU buffers."""
        if colors_dict is None:
            self.custom_atom_colors = {}
        else:
            self.custom_atom_colors.update(colors_dict)
            
        if self.molecule and self._colors is not None:
            old_to_new = getattr(self.molecule, 'properties', {}).get('_old_to_new_idx', {})
            for atom in self.molecule.atoms:
                idx = atom.index
                new_idx = old_to_new.get(idx, idx)
                if new_idx >= len(self._colors):
                    continue
                if idx in self.custom_atom_colors:
                    c = self.custom_atom_colors[idx]
                    if isinstance(c, tuple) and len(c) == 3:
                        if isinstance(c[0], int):
                            c = [x / 255.0 for x in c]
                    elif hasattr(c, 'redF'):
                        c = [c.redF(), c.greenF(), c.blueF()]
                    self._colors[new_idx] = c
                else:
                    self._colors[new_idx] = _element_color_float(atom.symbol)
            
            # Re-upload bond colors
            self._pack_bonds(self.molecule, self.molecule.atoms, self._colors) # The atoms array here is ignored for unpacking since we remap it inside _pack_bonds correctly now
            
            if self.gl_available:
                self._update_gl_buffers()
        self.update()

    def set_molecule(self, mol):
        """Load a molecule for rendering.

        Packs atom positions, colours, and radii into flat numpy arrays.
        """
        self.molecule = mol

        if mol is None or not hasattr(mol, 'atoms') or len(mol.atoms) == 0:
            self._positions = None
            self._colors = None
            self._radii = None
            self._symbols = None
            self._bond_starts = None
            self._bond_ends = None
            self._bond_start_colors = None
            self._bond_end_colors = None
            self.update()
            return

        # Split atoms into protein vs ligand to render ligand bonds/atoms over/during mesh rendering
        protein_atoms = []
        ligand_atoms = []
        for a in mol.atoms:
            if a.symbol == 'H' and not self.show_hydrogens:
                continue
            if getattr(a, 'is_hetatm', False):
                ligand_atoms.append(a)
            else:
                protein_atoms.append(a)
                
        # Remap atom indices for bonds
        old_to_new = {}
        ordered_atoms = protein_atoms + ligand_atoms
        for i, a in enumerate(ordered_atoms):
            old_to_new[a.index] = i
        mol.properties['_old_to_new_idx'] = old_to_new
        
        self._ligand_start = len(protein_atoms)
        n = len(ordered_atoms)

        positions = np.zeros((n, 3), dtype=np.float32)
        colors = np.zeros((n, 3), dtype=np.float32)
        radii = np.zeros(n, dtype=np.float32)
        symbols = []
        
        self._atom_is_backbone = np.zeros(n, dtype=bool)
        self._atom_is_sidechain = np.zeros(n, dtype=bool)
        self._atom_is_water = np.zeros(n, dtype=bool)
        self._atom_res_seqs = np.zeros(n, dtype=int)

        _AMINO_ACIDS = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
        _BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}
        _WATER_RES = {'HOH', 'WAT', 'SOL', 'DOD'}

        for i, atom in enumerate(ordered_atoms):
            if atom.has_coords:
                positions[i] = [atom.x, atom.y, atom.z if atom.z is not None else 0.0]
            c = _element_color_float(atom.symbol)
            colors[i] = c
            radii[i] = _display_radius(atom.symbol)
            symbols.append(atom.symbol)

            if (getattr(atom, 'res_name', '') or '').upper() in _WATER_RES:
                self._atom_is_water[i] = True

            rs = getattr(atom, 'res_seq', -1)
            if rs is not None:
                self._atom_res_seqs[i] = rs
                
            if i < self._ligand_start:
                # Protein atom (incl. modified residues like PTR / HID): classify
                # by atom name so its backbone joins the cartoon and its side
                # chain is hidden with the rest — never drawn as an overlapping
                # ball-and-stick "ligand".
                if (atom.pdb_name or '').strip() in _BACKBONE_ATOMS:
                    self._atom_is_backbone[i] = True
                else:
                    self._atom_is_sidechain[i] = True
            
        # Center to origin
        if n > 0:
            self._centroid = np.mean(positions, axis=0)
            positions -= self._centroid
        else:
            self._centroid = np.zeros(3, dtype=np.float32)

        self._positions = positions
        self._colors = colors
        self._radii = radii
        self._symbols = symbols

        # Pack bond data
        self._pack_bonds(mol, ordered_atoms, colors)

        # Auto-fit camera
        self._auto_fit()

        if self.gl_available:
            self._update_gl_buffers()

        self.update()

    def _update_gl_buffers(self):
        """Upload molecule data to GPU buffers."""
        if not self.gl_available: return
        import time
        t_upd = time.time()
        try:
            self.makeCurrent()
            gl = self.context().functions()
            
            # Calculate Atom Visibility Masks
            active_radii = self._radii * self.sphere_scale
            is_mesh_active = self.molecule and self.molecule.properties.get('is_protein')
            show_all_sidechains = getattr(self, 'show_sidechains', False)
            sidechain_res_vis = getattr(self, 'sidechain_res_vis', {})
            visible_sidechains = getattr(self, 'visible_sidechains', set())
            
            if is_mesh_active and hasattr(self, '_atom_is_backbone'):
                # Mask out all backbone atoms
                active_radii[self._atom_is_backbone] = 0.0

                # Mask out sidechains if not shown
                if not show_all_sidechains:
                    for i in range(len(active_radii)):
                        if self._atom_is_sidechain[i]:
                            rs = self._atom_res_seqs[i]
                            if not sidechain_res_vis.get(rs, False) and rs not in visible_sidechains:
                                active_radii[i] = 0.0

            # Hide crystallographic waters by default — a docking-prep structure
            # can carry hundreds of lone O atoms that render as scattered red dots
            # obscuring the protein.  (Set `show_waters = True` to reveal them.)
            if (not getattr(self, 'show_waters', False)
                    and hasattr(self, '_atom_is_water')):
                active_radii[self._atom_is_water] = 0.0

            # 1. Atoms Buffer (Center, Color, Radius, Offset)
            if self._vbo_atoms: self._vbo_atoms.destroy()
            n = len(self._positions)
            if n > 0:
                import numpy as np
                data = np.zeros(n * 6 * 9, dtype=np.float32)
                for i in range(6):
                    data[i*9+0::54] = self._positions[:, 0]
                    data[i*9+1::54] = self._positions[:, 1]
                    data[i*9+2::54] = self._positions[:, 2]
                    data[i*9+3::54] = self._colors[:, 0]
                    data[i*9+4::54] = self._colors[:, 1]
                    data[i*9+5::54] = self._colors[:, 2]
                    data[i*9+6::54] = active_radii
                offsets = np.array([[-1, -1], [1, -1], [1, 1], [-1, -1], [1, 1], [-1, 1]], dtype=np.float32)
                for i in range(6):
                    data[i*9+7::54] = offsets[i, 0]
                    data[i*9+8::54] = offsets[i, 1]
                from src.shared.qt_compat import QOpenGLBuffer
                self._vbo_atoms = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
                self._vbo_atoms.create()
                self._vbo_atoms.bind()
                self._vbo_atoms.allocate(data.tobytes(), len(data.tobytes()))
            
            # 2. Protein Mesh Buffer
            if self.molecule and self.molecule.properties.get('is_protein'):
                from src.features.visualization_3d.services.cartoon_generator import generate_cartoon_mesh
                from src.shared.ui.theme import COLORS
                # Compute secondary structure (DSSP) and propagate to atoms BEFORE
                # the mesh generator reads atom.ss_type.  Structures loaded without
                # HELIX/SHEET header records (PDBQT, many PDBs) arrive all-coil, so
                # without this the GL cartoon is a featureless worm instead of a
                # proper helix/sheet ribbon.  (The software renderer already does
                # this; the GL path was missing it.)  Guarded per-molecule so DSSP
                # runs once, not every frame.
                if not self.molecule.properties.get('_ss_propagated', False):
                    try:
                        from src.features.visualization_3d.services.protein_rendering import ProteinStructure
                        ProteinStructure(self.molecule)
                    except Exception:
                        pass
                    self.molecule.properties['_ss_propagated'] = True
                t_mesh = time.time()
                # Always use high quality mesh for GPU renderer (it can easily handle it)
                v, t, c = generate_cartoon_mesh(self.molecule, spline_steps=24, profile_detail=16, theme_colors=COLORS)
                if v is not None:
                    if self._vbo_mesh: self._vbo_mesh.destroy()
                    indices = t.flatten()
                    v_flat = v[indices] - self._centroid
                    c_flat = c[indices]
                    self._num_mesh_vertices = len(v_flat)
                    mdata = np.zeros(self._num_mesh_vertices * 9, dtype=np.float32)
                    mdata[0::9] = v_flat[:, 0]
                    mdata[1::9] = v_flat[:, 1]
                    mdata[2::9] = v_flat[:, 2]
                    mdata[3::9] = c_flat[:, 0]
                    mdata[4::9] = c_flat[:, 1]
                    mdata[5::9] = c_flat[:, 2]
                    # Vectorized normal computation (replaces slow per-triangle Python loop)
                    n_tris = self._num_mesh_vertices // 3
                    tri_verts = v_flat[:n_tris * 3].reshape(n_tris, 3, 3)  # (T, 3_verts, 3_xyz)
                    edge1 = tri_verts[:, 1, :] - tri_verts[:, 0, :]
                    edge2 = tri_verts[:, 2, :] - tri_verts[:, 0, :]
                    normals = np.cross(edge1, edge2)  # (T, 3)
                    mags = np.linalg.norm(normals, axis=1, keepdims=True)
                    mags[mags < 1e-6] = 1.0
                    normals /= mags
                    # Broadcast same normal to all 3 vertices of each triangle
                    normals_expanded = np.repeat(normals, 3, axis=0)  # (T*3, 3)
                    mdata[6::9] = normals_expanded[:, 0]
                    mdata[7::9] = normals_expanded[:, 1]
                    mdata[8::9] = normals_expanded[:, 2]
                    from src.shared.qt_compat import QOpenGLBuffer
                    self._vbo_mesh = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
                    self._vbo_mesh.create()
                    self._vbo_mesh.bind()
                    self._vbo_mesh.allocate(mdata.tobytes(), len(mdata.tobytes()))
                    self._ibo_mesh = None
                    # print(f"[GL] Mesh buffer updated in {time.time()-t_mesh:.3f}s")
            else:
                self._vbo_mesh = None
                
            # Calculate Bond Visibility Mask
            if self._bond_starts is not None and len(self._bond_starts) > 0:
                N_bonds = len(self._bond_starts)
                self._num_lines = N_bonds * 6
                
                line_data = np.empty((N_bonds, 6, 14), dtype=np.float32)
                
                starts = self._bond_starts.copy()
                ends = self._bond_ends.copy()
                
                if is_mesh_active and hasattr(self, '_bond_is_backbone'):
                    # Create boolean mask
                    bond_mask = np.ones(N_bonds, dtype=bool)
                    # Hide backbone bonds
                    if len(self._bond_is_backbone) == N_bonds:
                        bond_mask[self._bond_is_backbone] = False
                    
                    if not show_all_sidechains:
                        if not sidechain_res_vis and not visible_sidechains:
                            bond_mask[:self._ligand_bond_start // 6] = False
                        else:
                            # Use vectorized masking for explicitly visible sidechains
                            rs_array = self._bond_res_seq
                            mask_to_hide = np.ones(N_bonds, dtype=bool)
                            
                            for rs in visible_sidechains:
                                mask_to_hide[rs_array == rs] = False
                            for rs, vis in sidechain_res_vis.items():
                                if vis:
                                    mask_to_hide[rs_array == rs] = False
                            
                            # Hide sidechain bonds whose residue is NOT visible
                            # (we only apply this to protein bonds, i.e., before _ligand_bond_start)
                            is_protein_bond = np.arange(N_bonds) < (self._ligand_bond_start // 6)
                            bond_mask[mask_to_hide & (rs_array != -1) & is_protein_bond] = False
                    
                    starts[~bond_mask] = 0.0
                    ends[~bond_mask] = 0.0
                
                colors1 = self._bond_start_colors
                colors2 = self._bond_end_colors
                
                line_data[:, :, 0:3] = starts[:, np.newaxis, :]
                line_data[:, :, 3:6] = ends[:, np.newaxis, :]
                line_data[:, :, 6:9] = colors1[:, np.newaxis, :]
                line_data[:, :, 9:12] = colors2[:, np.newaxis, :]
                
                corners = np.array([
                    [-1.0, 0.0], [ 1.0, 0.0], [ 1.0, 1.0],
                    [-1.0, 0.0], [ 1.0, 1.0], [-1.0, 1.0]
                ], dtype=np.float32)
                line_data[:, :, 12:14] = corners
                
                line_data = line_data.reshape(-1, 14)
                
                from src.shared.qt_compat import QOpenGLBuffer
                self._vbo_lines = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
                self._vbo_lines.create()
                self._vbo_lines.bind()
                self._vbo_lines.allocate(line_data.tobytes(), len(line_data.tobytes()))
            else:
                self._vbo_lines = None
                self._num_lines = 0

            # print(f"[Performance] Total GL Buffer update took {time.time()-t_upd:.3f}s")

        except Exception as e:
            pass # print(f"[GL] Buffer update error: {e}")

    def _pack_bonds(self, mol, atoms, atom_colors):
        """Pack bond endpoint data into flat numpy arrays."""
        bonds = getattr(mol, 'bonds', None)
        if not bonds:
            self._bond_starts = None
            self._bond_ends = None
            self._bond_start_colors = None
            self._bond_end_colors = None
            self._ligand_bond_start = 0
            return

        n_atoms = len(atoms)
        
        prot_starts, prot_ends, prot_start_colors, prot_end_colors = [], [], [], []
        lig_starts, lig_ends, lig_start_colors, lig_end_colors = [], [], [], []

        old_to_new = getattr(mol, 'properties', {}).get('_old_to_new_idx', {})
        
        # Precompute rings for bond offsets (Skip for large structures to avoid 25s freeze)
        rings = []
        if hasattr(mol, 'find_rings') and len(atoms) < 1500:
            rings = mol.find_rings()
        bond_to_ring = {}
        for ring in rings:
            for i in range(len(ring)):
                a, b = ring[i], ring[(i+1) % len(ring)]
                bond_to_ring[frozenset([a, b])] = ring
                
        self._bond_is_backbone = []
        self._bond_res_seq = []
        
        for bond in bonds:
            bi = bond.begin_atom_idx
            ei = bond.end_atom_idx

            nbi = old_to_new.get(bi)
            nei = old_to_new.get(ei)

            if nbi is None or nei is None:
                continue
            
            # Map original indices to ordered indices
            new_bi = old_to_new.get(bi, bi)
            new_ei = old_to_new.get(ei, ei)
            
            if new_bi >= n_atoms or new_ei >= n_atoms:
                continue

            # Skip bonds touching a hidden water (see the radius mask above).
            if (not getattr(self, 'show_waters', False) and hasattr(self, '_atom_is_water')
                    and (self._atom_is_water[new_bi] or self._atom_is_water[new_ei])):
                continue

            # Fetch original atoms from the molecule for property checks
            a1 = mol.atoms[bi]
            a2 = mol.atoms[ei]
            
            if not a1.has_coords or not a2.has_coords:
                continue
            # Skip hydrogen bonds if hydrogens hidden
            if not self.show_hydrogens and (a1.symbol == 'H' or a2.symbol == 'H'):
                continue
            
            is_backbone_bond = self._atom_is_backbone[new_bi] and self._atom_is_backbone[new_ei]

            p1 = self._positions[new_bi]
            p2 = self._positions[new_ei]
            
            # Check if this bond belongs to the ligand
            is_ligand = getattr(a1, 'is_hetatm', False) or getattr(a2, 'is_hetatm', False)
            
            bond_lines = []
            if is_ligand and (bond.is_double or bond.is_aromatic or bond.order >= 1.5):
                # Calculate a 3D offset perpendicular to the bond
                d = p2 - p1
                length = np.linalg.norm(d)
                if length > 1e-4:
                    d = d / length
                    offset_dir = None
                    # Try to use ring center for offset direction
                    ring = bond_to_ring.get(frozenset([bi, ei]))
                    if ring:
                        center = np.zeros(3, dtype=np.float32)
                        valid_atoms = 0
                        for idx in ring:
                            r_new_idx = old_to_new.get(idx)
                            if r_new_idx is not None and r_new_idx < n_atoms:
                                center += self._positions[r_new_idx]
                                valid_atoms += 1
                        if valid_atoms > 0:
                            center /= valid_atoms
                            mid = (p1 + p2) / 2.0
                            to_center = center - mid
                            to_center -= np.dot(to_center, d) * d
                            norm = np.linalg.norm(to_center)
                            if norm > 1e-4:
                                offset_dir = to_center / norm
                    
                    if offset_dir is None:
                        # Fallback orthogonal vector
                        up = np.array([0.0, 1.0, 0.0], dtype=np.float32)
                        if abs(np.dot(d, up)) > 0.9:
                            up = np.array([1.0, 0.0, 0.0], dtype=np.float32)
                        offset_dir = np.cross(d, up)
                        offset_dir /= np.linalg.norm(offset_dir)
                    
                    offset = offset_dir * 0.15 * self.stick_scale
                    
                    if bond.is_aromatic or bond.order == 1.5:
                        bond_lines.append((p1, p2))
                        bond_lines.append((p1 + offset, p2 + offset))
                    else:
                        bond_lines.append((p1 + offset*0.7, p2 + offset*0.7))
                        bond_lines.append((p1 - offset*0.7, p2 - offset*0.7))
                else:
                    bond_lines.append((p1, p2))
            else:
                bond_lines.append((p1, p2))
            
            for start_p, end_p in bond_lines:
                if is_ligand:
                    lig_starts.append(start_p)
                    lig_ends.append(end_p)
                    lig_start_colors.append(atom_colors[new_bi])
                    lig_end_colors.append(atom_colors[new_ei])
                else:
                    prot_starts.append(start_p)
                    prot_ends.append(end_p)
                    prot_start_colors.append(atom_colors[new_bi])
                    prot_end_colors.append(atom_colors[new_ei])
                    self._bond_is_backbone.append(is_backbone_bond)
                    
                    rs = getattr(a1, 'res_seq', -1)
                    if rs is None or rs == -1:
                        rs = getattr(a2, 'res_seq', -1)
                    if rs is None:
                        rs = -1   # ligand/non-protein atoms carry no res_seq
                    self._bond_res_seq.append(rs)

        self._ligand_bond_start = len(prot_starts) * 6
        
        self._bond_is_backbone.extend([False] * len(lig_starts))
        self._bond_is_backbone = np.array(self._bond_is_backbone, dtype=bool)
        
        self._bond_res_seq.extend([-1] * len(lig_starts))
        self._bond_res_seq = np.array(self._bond_res_seq, dtype=int)

        starts = prot_starts + lig_starts
        ends = prot_ends + lig_ends
        start_colors = prot_start_colors + lig_start_colors
        end_colors = prot_end_colors + lig_end_colors

        if starts:
            self._bond_starts = np.array(starts, dtype=np.float32)
            self._bond_ends = np.array(ends, dtype=np.float32)
            self._bond_start_colors = np.array(start_colors, dtype=np.float32)
            self._bond_end_colors = np.array(end_colors, dtype=np.float32)
        else:
            self._bond_starts = None
            self._bond_ends = None
            self._bond_start_colors = None
            self._bond_end_colors = None

