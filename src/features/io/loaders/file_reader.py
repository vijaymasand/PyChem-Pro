"""
File Readers — Import molecules from MOL, SDF, and MOL2 file formats.

Reads MDL V2000 MOL/SDF and Tripos MOL2 formats into Molecule objects.
"""

import re
from src.core.domain.models.atom import Atom
from src.core.domain.models.bond import Bond, BondType
from src.core.domain.models.molecule import Molecule
from src.features.smiles_parser.rules.aromaticity import perceive_aromaticity


def _post_load_processing(mol: Molecule):
    """
    Perform essential perception and typing after loading from file.
    """
    # 1. Ring perception (essential for aromaticity and descriptors)
    mol.find_rings()
    
    # 2. Propagate explicit aromaticity from bonds (MOL2/SDF type 4)
    mol.propagate_aromaticity()
    
    # 3. Perceive aromaticity from Kekule structures (SDF/MOL alternating bonds)
    # We assign hybridization first as it's a prerequisite for perception
    mol.assign_hybridization()
    perceive_aromaticity(mol)
    
    # 4. Final typing
    mol.assign_sybyl_types()


def read_mol(filepath_or_string, is_string=False):
    """
    Read a molecule from MDL V2000 MOL file.

    Args:
        filepath_or_string: File path or MOL content string
        is_string: If True, treat first arg as content string

    Returns:
        Molecule object with atoms, bonds, and coordinates
    """
    if is_string:
        lines = filepath_or_string.split('\n')
    else:
        with open(filepath_or_string, 'r') as f:
            lines = f.read().split('\n')

    if len(lines) < 4:
        raise ValueError("MOL file too short")

    # Header: line 0 = name, line 1 = program/timestamp, line 2 = comment
    mol_name = lines[0].strip()
    mol = Molecule(name=mol_name)
    mol.properties['source'] = 'file'

    # Counts line (line 3)
    counts_line = lines[3]
    num_atoms = int(counts_line[0:3].strip())
    num_bonds = int(counts_line[3:6].strip())

    # Atom block: lines 4 to 4+num_atoms-1
    for i in range(num_atoms):
        line = lines[4 + i]
        x = float(line[0:10].strip())
        y = float(line[10:20].strip())
        z = float(line[20:30].strip())
        symbol = line[31:34].strip()

        # Charge code
        charge_code = 0
        if len(line) >= 39:
            try:
                charge_code = int(line[36:39].strip())
            except ValueError:
                pass

        formal_charge = _mdl_charge_decode(charge_code)

        atom = Atom(symbol, formal_charge=formal_charge)
        atom.x = x
        atom.y = y
        atom.z = z
        mol.add_atom(atom)

    # Bond block: lines 4+num_atoms to 4+num_atoms+num_bonds-1
    bond_start = 4 + num_atoms
    for i in range(num_bonds):
        line = lines[bond_start + i]
        a1 = int(line[0:3].strip()) - 1  # Convert to 0-indexed
        a2 = int(line[3:6].strip()) - 1
        bond_type_code = int(line[6:9].strip())

        stereo = 0
        if len(line) >= 12:
            try:
                stereo = int(line[9:12].strip())
            except ValueError:
                pass

        bt = _mdl_bond_type(bond_type_code)
        mol.add_bond(a1, a2, bt, stereo)

    # Properties block — look for M  CHG and M  ISO
    for i in range(bond_start + num_bonds, len(lines)):
        line = lines[i]
        if line.startswith('M  END'):
            break
        if line.startswith('M  CHG'):
            num_entries = int(line[6:9].strip())
            for j in range(num_entries):
                offset = 9 + j * 8
                atom_idx = int(line[offset:offset+4].strip()) - 1
                charge = int(line[offset+4:offset+8].strip())
                if 0 <= atom_idx < len(mol.atoms):
                    mol.atoms[atom_idx].formal_charge = charge
        elif line.startswith('M  ISO'):
            num_entries = int(line[6:9].strip())
            for j in range(num_entries):
                offset = 9 + j * 8
                atom_idx = int(line[offset:offset+4].strip()) - 1
                isotope = int(line[offset+4:offset+8].strip())
                if 0 <= atom_idx < len(mol.atoms):
                    mol.atoms[atom_idx].isotope = isotope

    _post_load_processing(mol)
    return mol


def read_sdf(filepath_or_string, is_string=False):
    """
    Read the first molecule from an SDF file.

    Args:
        filepath_or_string: File path or SDF content string
        is_string: If True, treat first arg as content string

    Returns:
        Molecule object
    """
    if is_string:
        content = filepath_or_string
    else:
        with open(filepath_or_string, 'r') as f:
            content = f.read()

    # Split by $$$$ and take first block
    blocks = content.split('$$$$')
    if not blocks or not blocks[0].strip():
        raise ValueError("No molecule found in SDF file")

    mol_block = blocks[0].strip()

    # Read as MOL
    mol = read_mol(mol_block, is_string=True)

    # Parse data fields
    data_section = mol_block.split('M  END')
    if len(data_section) > 1:
        data_text = data_section[1]
        field_pattern = re.compile(r'>\s*<([^>]+)>\s*\n(.+?)(?=\n>|\n\s*$)', re.DOTALL)
        for match in field_pattern.finditer(data_text):
            key = match.group(1).strip()
            value = match.group(2).strip()
            mol.properties[key] = value

    _post_load_processing(mol)
    return mol


def read_mol2(filepath_or_string, is_string=False):
    """
    Read a molecule from Tripos MOL2 file.

    Args:
        filepath_or_string: File path or MOL2 content string
        is_string: If True, treat first arg as content string

    Returns:
        Molecule object with atoms, bonds, coordinates, and partial charges
    """
    if is_string:
        content = filepath_or_string
    else:
        with open(filepath_or_string, 'r') as f:
            content = f.read()

    lines = content.split('\n')

    # Find section markers
    sections = {}
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('@<TRIPOS>'):
            section_name = stripped.replace('@<TRIPOS>', '')
            sections[section_name] = i

    mol_name = ""
    if 'MOLECULE' in sections:
        mol_start = sections['MOLECULE']
        if mol_start + 1 < len(lines):
            mol_name = lines[mol_start + 1].strip()

    mol = Molecule(name=mol_name)
    mol.properties['source'] = 'mol2'

    # Parse ATOM section
    if 'ATOM' not in sections:
        raise ValueError("No ATOM section found in MOL2 file")

    atom_start = sections['ATOM'] + 1
    atom_end = len(lines)
    # Find end of atom section (next section marker or EOF)
    for i in range(atom_start, len(lines)):
        if lines[i].strip().startswith('@<TRIPOS>'):
            atom_end = i
            break

    for i in range(atom_start, atom_end):
        line = lines[i].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue

        # atom_id atom_name x y z atom_type [subst_id subst_name charge]
        atom_name = parts[1]
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        atom_type = parts[5]

        # Extract element symbol from SYBYL type (e.g., "C.3" -> "C", "N.am" -> "N")
        symbol = atom_type.split('.')[0]
        # Handle special cases
        if symbol.upper() in ('CL', 'BR'):
            symbol = symbol[0].upper() + symbol[1].lower()
        elif len(symbol) == 1:
            symbol = symbol.upper()
        else:
            symbol = symbol[0].upper() + symbol[1:].lower()

        atom = Atom(symbol)
        atom.x = x
        atom.y = y
        atom.z = z
        atom.sybyl_type = atom_type

        # Parse partial charge
        if len(parts) >= 9:
            try:
                atom.partial_charge = float(parts[8])
            except ValueError:
                pass

        mol.add_atom(atom)

    # Parse BOND section
    if 'BOND' in sections:
        bond_start = sections['BOND'] + 1
        bond_end = len(lines)
        for i in range(bond_start, len(lines)):
            if lines[i].strip().startswith('@<TRIPOS>'):
                bond_end = i
                break

        for i in range(bond_start, bond_end):
            line = lines[i].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue

            # bond_id origin_atom_id target_atom_id bond_type
            a1 = int(parts[1]) - 1  # 0-indexed
            a2 = int(parts[2]) - 1
            bond_type_str = parts[3].lower()

            if bond_type_str == '2':
                bt = BondType.DOUBLE
            elif bond_type_str == '3':
                bt = BondType.TRIPLE
            elif bond_type_str == 'ar':
                bt = BondType.AROMATIC
            elif bond_type_str == 'am':
                bt = BondType.AMIDE
            else:
                bt = BondType.SINGLE

            if 0 <= a1 < len(mol.atoms) and 0 <= a2 < len(mol.atoms):
                mol.add_bond(a1, a2, bt)

    _post_load_processing(mol)
    return mol


def _mdl_charge_decode(charge_code):
    """Convert MDL V2000 charge code to formal charge."""
    code_map = {
        0: 0,
        1: 3,
        2: 2,
        3: 1,
        5: -1,
        6: -2,
        7: -3,
    }
    return code_map.get(charge_code, 0)


def _mdl_bond_type(code):
    """Convert MDL bond type code to BondType."""
    if code == 1:
        return BondType.SINGLE
    elif code == 2:
        return BondType.DOUBLE
    elif code == 3:
        return BondType.TRIPLE
    elif code == 4:
        return BondType.AROMATIC
    else:
        return BondType.SINGLE


def read_pdb(filepath):
    """
    Read a molecule/protein from PDB file format.
    
    Highly optimized for large proteins using NumPy for coordinate storage
    and batch processing. Handles 100k+ atoms in seconds.
    """
    import numpy as np
    import os
    import time
    
    t0 = time.time()
    mol = Molecule(name=os.path.basename(filepath))
    conect_records = []
    helix_ranges = []
    sheet_ranges = []
    serial_to_idx = {}
    
    atoms_data = [] # List of (serial, element, x, y, z, name, res, chain, seq, b_fact, is_het)
    
    # Use fast I/O
    with open(filepath, 'r', buffering=1024*1024) as f:
        for line in f:
            if len(line) < 6: continue
            record = line[0:6]
            
            if record.startswith('ATOM') or record.startswith('HETATM'):
                try:
                    serial = int(line[6:11])
                    name = line[12:16].strip()
                    res = line[17:20].strip()
                    chain = line[21]
                    seq = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    element = line[76:78].strip() if len(line) >= 78 else ''
                    if not element:
                        nc = name.lstrip('0123456789')
                        element = nc[:2].upper() if nc[:2].upper() in ('CL','BR','FE','ZN','MG','CA','NA','MN','CU','CO','NI') else nc[0].upper()
                    
                    b_fact = float(line[60:66]) if len(line) >= 66 else 0.0
                    is_het = record.startswith('HETATM')
                    
                    atoms_data.append((serial, element, x, y, z, name, res, chain, seq, b_fact, is_het))
                except: continue
                
            elif record.startswith('HELIX'):
                try: helix_ranges.append((line[19], int(line[21:25]), int(line[33:37])))
                except: pass
            elif record.startswith('SHEET'):
                try: sheet_ranges.append((line[21], int(line[22:26]), int(line[33:37])))
                except: pass
            elif record.startswith('CONECT'):
                try:
                    p = line[6:].split()
                    if len(p) >= 2:
                        o = int(p[0])
                        for t in p[1:]: conect_records.append((o, int(t)))
                except: pass

    mol.begin_bulk_load()
    for d in atoms_data:
        a = Atom(d[1])
        a.x, a.y, a.z = d[2], d[3], d[4]
        a.pdb_name, a.res_name, a.chain_id, a.res_seq = d[5], d[6], d[7], d[8]
        a.b_factor, a.is_hetatm = d[9], d[10]
        idx = mol.add_atom(a)
        serial_to_idx[d[0]] = idx

    # Fast SS assignment
    if helix_ranges or sheet_ranges:
        h_dict = {}
        for c, s, e in helix_ranges: h_dict.setdefault(c, []).append((s, e))
        s_dict = {}
        for c, s, e in sheet_ranges: s_dict.setdefault(c, []).append((s, e))
        
        for a in mol.atoms:
            cid, rseq = a.chain_id, a.res_seq
            a.ss_type = 'C'
            if cid in h_dict:
                for s, e in h_dict[cid]:
                    if s <= rseq <= e: a.ss_type = 'H'; break
            if a.ss_type == 'C' and cid in s_dict:
                for s, e in s_dict[cid]:
                    if s <= rseq <= e: a.ss_type = 'E'; break
    else:
        for a in mol.atoms: a.ss_type = 'C'

    # Add CONECT bonds
    done = set()
    for s1, s2 in conect_records:
        if s1 in serial_to_idx and s2 in serial_to_idx:
            i1, i2 = serial_to_idx[s1], serial_to_idx[s2]
            key = tuple(sorted((i1, i2)))
            if key not in done:
                mol.add_bond(i1, i2, BondType.SINGLE)
                done.add(key)

    if len(mol.atoms) > 0 and len(mol.bonds) < len(mol.atoms) * 0.8:
        _auto_bond_pdb(mol)
    
    mol.end_bulk_load()
    mol.properties.update({'helix_ranges': helix_ranges, 'sheet_ranges': sheet_ranges,
                           'is_protein': any(a.res_name in _AMINO_ACIDS for a in mol.atoms)})
    print(f"[Performance] read_pdb took {time.time()-t0:.3f}s for {len(mol.atoms)} atoms")
    return mol

    # Enable bulk-load mode to suppress ring-cache invalidation per atom/bond
    mol.begin_bulk_load()

    # Batch add atoms for better performance
    for serial, atom in atoms_to_add:
        idx = mol.add_atom(atom)
        serial_to_idx[serial] = idx

    # Efficient secondary structure assignment using range dictionaries
    if helix_ranges or sheet_ranges:
        # Convert ranges to dictionaries for O(1) lookup
        helix_dict = {}
        for chain, start, end in helix_ranges:
            if chain not in helix_dict:
                helix_dict[chain] = []
            helix_dict[chain].append((start, end))
            
        sheet_dict = {}
        for chain, start, end in sheet_ranges:
            if chain not in sheet_dict:
                sheet_dict[chain] = []
            sheet_dict[chain].append((start, end))
        
        # Assign secondary structure efficiently
        for atom in mol.atoms:
            if not hasattr(atom, 'chain_id') or not atom.chain_id:
                atom.ss_type = 'C'
                continue
                
            chain_id = atom.chain_id
            res_seq = atom.res_seq
            
            # Check helix
            if chain_id in helix_dict:
                for start, end in helix_dict[chain_id]:
                    if start <= res_seq <= end:
                        atom.ss_type = 'H'
                        break
                else:
                    # Check sheet
                    if chain_id in sheet_dict:
                        for start, end in sheet_dict[chain_id]:
                            if start <= res_seq <= end:
                                atom.ss_type = 'E'
                                break
                        else:
                            atom.ss_type = 'C'
                    else:
                        atom.ss_type = 'C'
            else:
                atom.ss_type = 'C'
    else:
        # No secondary structure data - set all to coil
        for atom in mol.atoms:
            atom.ss_type = 'C'

    # Add bonds from CONECT records (optimized)
    added_bonds = set()
    for serial1, serial2 in conect_records:
        if serial1 in serial_to_idx and serial2 in serial_to_idx:
            idx1 = serial_to_idx[serial1]
            idx2 = serial_to_idx[serial2]
            bond_key = tuple(sorted((idx1, idx2)))
            if bond_key not in added_bonds:
                mol.add_bond(idx1, idx2, BondType.SINGLE)
                added_bonds.add(bond_key)

    # Auto-detect bonds if CONECT records are missing or incomplete (e.g. only for ligands)
    if len(mol.atoms) > 0 and len(mol.bonds) < len(mol.atoms) * 0.8:
        _auto_bond_pdb(mol)

    # End bulk-load mode — invalidate ring cache once
    mol.end_bulk_load()

    # Store secondary structure ranges on the molecule
    mol.properties['helix_ranges'] = helix_ranges
    mol.properties['sheet_ranges'] = sheet_ranges
    mol.properties['is_protein'] = any(
        hasattr(a, 'res_name') and a.res_name in _AMINO_ACIDS
        for a in mol.atoms)

    return mol


_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEC', 'PYL',
}

# Covalent bond distance cutoffs by element pair (Angstroms)
_COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'P': 1.07,
    'F': 0.57, 'Cl': 0.99, 'Br': 1.14, 'I': 1.33, 'Se': 1.20,
    'Fe': 1.32, 'Zn': 1.22, 'Mg': 1.41, 'Ca': 1.76, 'Na': 1.66,
    'Mn': 1.39, 'Cu': 1.32, 'Co': 1.26, 'Ni': 1.24,
}


def _auto_bond_pdb(mol):
    """
    Auto-detect bonds in PDB with NumPy vectorization.
    Handles 50,000+ atoms in under a second.
    """
    import numpy as np
    import time
    from collections import defaultdict

    t_start = time.time()
    n = len(mol.atoms)
    if n == 0: return

    # 1. Residue-topology bonding (Backbone)
    residue_atoms = defaultdict(dict)
    for i, atom in enumerate(mol.atoms):
        c, s, n_ = getattr(atom, 'chain_id',''), getattr(atom, 'res_seq',None), getattr(atom, 'pdb_name',None)
        if s is not None and n_:
            name = n_.strip()
            if name: residue_atoms[(c, s)][name] = i

    added = set()
    def _add(i, j):
        if i == j: return
        key = tuple(sorted((i, j)))
        if key not in added:
            mol.add_bond(i, j, BondType.SINGLE)
            added.add(key)

    for (c, s), d in residue_atoms.items():
        for n1, n2 in [('N','CA'), ('CA','C'), ('C','O'), ('C','OXT')]:
            if n1 in d and n2 in d: _add(d[n1], d[n2])
        # Peptide bond C(i) -> N(i+1)
        if 'C' in d:
            nxt = residue_atoms.get((c, s+1))
            if nxt and 'N' in nxt:
                a1, a2 = mol.atoms[d['C']], mol.atoms[nxt['N']]
                if ((a1.x-a2.x)**2 + (a1.y-a2.y)**2 + (a1.z-a2.z)**2)**0.5 < 2.0:
                    _add(d['C'], nxt['N'])

    # 2. Spatial Grid with NumPy Vectorization
    if n > 100000: return # Safety cap
    
    pos = np.array([[a.x, a.y, a.z] for a in mol.atoms], dtype=np.float32)
    radii = np.array([_COVALENT_RADII.get(a.symbol, 1.5) for a in mol.atoms], dtype=np.float32)
    
    cell_size = 4.0
    grid = defaultdict(list)
    indices = np.floor(pos / cell_size).astype(np.int32)
    for i in range(n): grid[tuple(indices[i])].append(i)

    offsets = np.array(np.meshgrid([-1,0,1],[-1,0,1],[-1,0,1])).T.reshape(-1, 3)
    
    for cell, cell_atoms in grid.items():
        cell_atoms = np.array(cell_atoms)
        # Check current cell and neighbors
        for off in offsets:
            nb_cell = tuple(np.array(cell) + off)
            if nb_cell not in grid: continue
            nb_atoms = np.array(grid[nb_cell])
            
            # Vectorized distance check
            for i in cell_atoms:
                # To avoid N^2 in cell, only check j > i
                targets = nb_atoms[nb_atoms > i]
                if len(targets) == 0: continue
                
                d2 = np.sum((pos[targets] - pos[i])**2, axis=1)
                max_d = (radii[targets] + radii[i] + 0.45)**2
                
                mask = (d2 > 0.16) & (d2 < max_d)
                for j in targets[mask]:
                    _add(i, int(j))
    print(f"[Performance] _auto_bond_pdb took {time.time()-t_start:.3f}s")
