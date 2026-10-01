"""
Molecular Fingerprints — Morgan (ECFP-like) and Topological pattern fingerprints.
"""

import hashlib
from collections import defaultdict


def _hash_ints(ints):
    """Hash a list of integers into a 32-bit integer."""
    # Convert integers to a packed string or simply use hashlib
    s = ",".join(str(i) for i in ints).encode('utf-8')
    h = hashlib.md5(s).digest()
    # Take first 4 bytes as a 32-bit int
    return int.from_bytes(h[:4], 'little')


def get_morgan_fingerprint(mol, radius=2, n_bits=2048):
    """
    Generate a Morgan (ECFP-like) fingerprint bit vector.

    Args:
        mol: Molecule object
        radius: Number of iterations (Radius 2 = ECFP4)
        n_bits: Size of the exact bit vector
    
    Returns:
        List of integers (0 or 1) of length n_bits
    """
    if not mol.atoms:
        return [0] * n_bits

    bonds = mol.bonds
    ring_atoms = set()
    rings = mol.find_rings()
    for ring in rings:
        ring_atoms.update(ring)

    # 1. Initialize invariants (Daylight-like)
    # [num_heavy_neighbors, valence, atomic_num, isotope, formal_charge, is_ring]
    identifiers = {}
    for atom in mol.atoms:
        deg = atom.heavy_degree(bonds)
        val = atom.total_valence(bonds) - atom.total_h
        invar = [
            deg,
            val,
            atom.atomic_number,
            atom.isotope,
            atom.formal_charge,
            1 if atom.index in ring_atoms else 0
        ]
        identifiers[atom.index] = _hash_ints(invar)
        
    # Store all unique identifiers generated across all iterations
    all_identifiers = set(identifiers.values())

    # 2. Iterate up to radius
    for _ in range(radius):
        new_identifiers = {}
        for atom in mol.atoms:
            # Collect neighbor identifiers and bond orders
            neighbors = []
            for nb_idx, bond_idx in mol._adjacency.get(atom.index, []):
                bond = mol.get_bond(bond_idx)
                neighbors.append((bond.order, identifiers[nb_idx]))
                
            # Sort neighbors to ensure canonical order
            neighbors.sort()
            
            # Create feature tuple
            feature = [identifiers[atom.index]]
            for bo, nb_id in neighbors:
                feature.extend([int(bo), nb_id])
                
            new_id = _hash_ints(feature)
            new_identifiers[atom.index] = new_id
            all_identifiers.add(new_id)
            
        identifiers = new_identifiers

    # 3. Fold into bit vector
    bit_vector = [0] * n_bits
    for ident in all_identifiers:
        bit_idx = ident % n_bits
        bit_vector[bit_idx] = 1

    return bit_vector


def get_topological_fingerprint(mol, min_path=1, max_path=7, n_bits=2048):
    """
    Generate a Daylight-like topological (path-based) fingerprint bit vector.
    
    Args:
        mol: Molecule object
        min_path: Minimum path length in bonds
        max_path: Maximum path length in bonds
        n_bits: Size of bit vector
        
    Returns:
        List of integers (0 or 1) of length n_bits
    """
    if not mol.atoms:
        return [0] * n_bits

    paths = set()

    def dfs(current_node, current_path, current_bond_path):
        length = len(current_bond_path)
        if length >= min_path:
            # Canonicalize path to ensure forward/reverse are identical
            # Node features: atomic number
            # Edge features: bond order
            
            # Build forward sequence
            seq_fwd = []
            for i in range(len(current_path)):
                atom = mol.get_atom(current_path[i])
                seq_fwd.append(atom.atomic_number)
                if i < len(current_bond_path):
                    bond = mol.get_bond(current_bond_path[i])
                    seq_fwd.append(int(bond.order) + (10 if bond.is_aromatic else 0))
                    
            # Build reverse sequence
            seq_rev = seq_fwd[::-1]
            
            # Take lexicographically smaller one
            canonical_seq = min(seq_fwd, seq_rev)
            paths.add(tuple(canonical_seq))

        if length < max_path:
            for nb_idx, bond_idx in mol._adjacency.get(current_node, []):
                if nb_idx not in current_path:  # Avoid cycles
                    dfs(nb_idx, current_path + [nb_idx], current_bond_path + [bond_idx])

    for atom in mol.atoms:
        dfs(atom.index, [atom.index], [])

    # Fold into bit vector
    bit_vector = [0] * n_bits
    for path_seq in paths:
        h = _hash_ints(path_seq)
        bit_idx = h % n_bits
        bit_vector[bit_idx] = 1

    return bit_vector


def get_maccs_keys(mol):
    """
    Generate MACCS 166 structural keys.
    (Note: Full MACCS 166 implementation requires exactly 166 SMARTS/substructure definitions.
    This provides a highly-simplified subset mapping just to demonstrate the logic.
    In a real-world scenario, you would evaluate full SMARTS matches.)
    
    Returns:
        List of integers (0 or 1) of length 167 (1-indexed, index 0 is 0)
    """
    keys = [0] * 167
    if not mol.atoms:
        return keys

    # Simple structural counts as proxy for MACCS keys
    symbols = {a.symbol for a in mol.atoms}
    
    if 'F' in symbols or 'Cl' in symbols or 'Br' in symbols or 'I' in symbols:
        keys[114] = 1  # Halogen
    
    if 'O' in symbols:
        keys[139] = 1  # OH
        
    if 'N' in symbols:
        keys[161] = 1  # Nitrogen
        
    if 'S' in symbols:
        keys[137] = 1  # Sulfur
        
    if 'P' in symbols:
        keys[142] = 1  # Phosphorus

    rings = mol.find_rings()
    if len(rings) > 0:
        keys[162] = 1  # Ring
    
    if any(a.is_aromatic for a in mol.atoms):
        keys[165] = 1  # Aromatic ring
        
    # Just a placeholder subset!
    return keys


def tanimoto_similarity(fp1, fp2):
    """
    Calculate Tanimoto similarity coefficient between two bit vectors.
    Sim(A, B) = c / (a + b - c)
    where a is bits in A, b is bits in B, and c is bits in A intersection B.
    """
    if len(fp1) != len(fp2):
        raise ValueError("Fingerprints must be the same length")
        
    c = sum(1 for i, j in zip(fp1, fp2) if i and j)
    a = sum(1 for i in fp1 if i)
    b = sum(1 for i in fp2 if i)
    
    denominator = (a + b - c)
    if denominator == 0:
        return 0.0
    return c / denominator
