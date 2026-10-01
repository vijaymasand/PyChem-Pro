"""
Substructure Matcher Service — Uses exactly subgraph isomorphism (DFS backtracking) 
to match a query SMILES/SMARTS string against a target Molecule.
"""

from src.features.smiles_parser.services.parser import parse_smiles


def find_substructure_matches(target_mol, smiles_query) -> list:
    """
    Find instances of smiles_query within target_mol.
    
    Returns a list of lists of atom indices from the target_mol that match the query subgraph.
    """
    if not target_mol or not target_mol.atoms:
        return []

    try:
        # Parse SMILES query to a Molecule
        query_mol = parse_smiles(smiles_query, name="query", use_bkchem_tokenizer=True)
    except Exception as e:
        # Invalid query
        return []
        
    if not query_mol or not query_mol.atoms:
        return []

    # Pre-calculate query atom degrees and neighbors for fast lookup
    query_neighbors = []
    for q_atom in query_mol.atoms:
        query_neighbors.append([n for n, b in query_mol.get_neighbor_bonds(q_atom.index)])
        
    target_neighbors = []
    for t_atom in target_mol.atoms:
        target_neighbors.append([n for n, b in target_mol.get_neighbor_bonds(t_atom.index)])

    # Ensure ring perception is up-to-date in target for relaxed matching
    target_mol.find_rings()

    def is_atom_compatible(q_atom, t_atom):
        if q_atom.symbol != '*':
            # Handle aromatic and aliphatic matches
            q_sym = q_atom.symbol.lower()
            t_sym = t_atom.symbol.lower()
            if q_sym != t_sym:
                return False
                
        # Relax aromatic matching: structural files don't flag is_aromatic natively.
        # Fallback to is_in_ring logic.
        if q_atom.is_aromatic and not (t_atom.is_aromatic or target_mol.is_in_ring(t_atom.index)):
            return False
            
        return True

    def is_bond_compatible(q_bond, t_bond, q1, q2, t1, t2):
        if q_bond is None or t_bond is None:
            return False
        if q_bond.bond_type == t_bond.bond_type:
            return True
            
        is_q_aromatic = q_bond.is_aromatic or (query_mol.atoms[q1].is_aromatic and query_mol.atoms[q2].is_aromatic)
        is_t_aromatic = t_bond.is_aromatic or t_bond.bond_type == 4 or (target_mol.atoms[t1].is_aromatic and target_mol.atoms[t2].is_aromatic)
        
        if is_q_aromatic:
            # Relaxed matching: query aromatic bond matches target bond if target is marked aromatic
            # OR if the target bond is simply in a ring (fallback for structural files).
            if is_t_aromatic or t_bond.is_in_ring:
                return True
            
        # Allow SMARTS single/double bonds to loosely match aromatic bonds in target
        if is_t_aromatic and (q_bond.bond_type in (1, 2)):
            return True
            
        return False

    all_matches = []
    
    def backtrack(current_match):
        """
        DFS backtracking search.
        current_match: dict {q_idx: t_idx}
        """
        if len(current_match) == len(query_mol.atoms):
            all_matches.append(list(current_match.values()))
            return
            
        # Select next query atom to match (heuristic: connected to already matched to aggressively prune)
        next_q = None
        for matched_q in current_match:
            for neighbor in query_neighbors[matched_q]:
                if neighbor not in current_match:
                    next_q = neighbor
                    break
            if next_q is not None:
                break
                
        if next_q is None:
            # First atom or new component
            for i in range(len(query_mol.atoms)):
                if i not in current_match:
                    next_q = i
                    break
                    
        q_atom = query_mol.atoms[next_q]
        
        # Determine candidate target atoms
        candidates = []
        if any(neighbor in current_match for neighbor in query_neighbors[next_q]):
            for neighbor in query_neighbors[next_q]:
                if neighbor in current_match:
                    matched_t = current_match[neighbor]
                    for t_neighbor in target_neighbors[matched_t]:
                        if t_neighbor not in current_match.values():
                            candidates.append(t_neighbor)
            # Remove duplicates
            candidates = list(set(candidates))
        else:
            candidates = [i for i in range(len(target_mol.atoms)) if i not in current_match.values()]
            
        for t_idx in candidates:
            t_atom = target_mol.atoms[t_idx]
            
            if not is_atom_compatible(q_atom, t_atom):
                continue
                
            # Check bond compatibilities with ALL currently valid matched neighbors
            bonds_compatible = True
            for neighbor in query_neighbors[next_q]:
                if neighbor in current_match:
                    matched_t = current_match[neighbor]
                    q_bond = query_mol.get_bond_between(next_q, neighbor)
                    t_bond = target_mol.get_bond_between(t_idx, matched_t)
                    if not is_bond_compatible(q_bond, t_bond, next_q, neighbor, t_idx, matched_t):
                        bonds_compatible = False
                        break
            
            if bonds_compatible:
                current_match[next_q] = t_idx
                backtrack(current_match)
                del current_match[next_q]

    backtrack({})
    
    # Remove duplicate combinations (e.g. matching symmetry in benzene returns 12 maps, but 1 unique subset)
    unique_matches = []
    seen = set()
    for m in all_matches:
        frozen = frozenset(m)
        if frozen not in seen:
            seen.add(frozen)
            unique_matches.append(m)
            
    return unique_matches
