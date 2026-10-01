"""
Professional Golden Layout (SMILES Roundtrip) Implementation.

This module provides a 'Gold Standard' 2D layout by converting any molecule
to a SMILES string and back. This leverages OASA's highly optimized
template-growth engine (used during SMILES parsing) which is far more 
robust than standalone coordinate generation.
"""

from src.vendors.oasa_bridge import domain_to_oasa_mol, oasa_prepare_for_layout
import src.vendors.oasa.smiles as oasa_smiles
import src.vendors.oasa.coords_generator as oasa_cg

def generate_golden_layout(molecule):
    """
    Generate professional 2D coordinates using the SMILES roundtrip method.

    Args:
        molecule: Domain molecule object.

    Returns:
        dict: {atom_index: [x, y]} for the input molecule.
    """
    if not molecule or not molecule.atoms:
        return {}

    # 1. Convert Domain -> OASA (O1) and kekulize before SMILES/layout
    o1_mol, atom_map = domain_to_oasa_mol(molecule)
    oasa_prepare_for_layout(o1_mol)

    # 2. Derive SMILES from O1
    # We use a dedicated smiles instance to capture the atom processing order
    sm_engine = oasa_smiles.smiles()

    # We must ensure the molecule is connected for standard SMILES
    if not o1_mol.is_connected():
        # Fallback to standard OASA calculation if disconnected (salts, etc.)
        oasa_cg.calculate_coords(o1_mol, bond_length=1.0, force=1)
        coords = {}
        for idx, v in atom_map.items():
            coords[idx] = [float(v.x), float(v.y)]
        return coords
        
    try:
        # Generate SMILES. This internally fills sm_engine._processed_atoms
        smiles_str = sm_engine.get_smiles(o1_mol)
        
        # Mapping: SMILES_POSITION (heavy atoms only) -> ORIGINAL_DOMAIN_INDEX
        oasa_to_domain = {v: k for k, v in atom_map.items()}
        smiles_to_domain_map = []
        for v in sm_engine._processed_atoms:
            # We ONLY map atoms that were in our initial heavy atom map
            if v in oasa_to_domain:
                smiles_to_domain_map.append(oasa_to_domain[v])
            else:
                # This could be an OASA-generated hydrogen atom - skip for mapping
                pass
        
        # 3. Parse SMILES back to OASA (O2) with high-quality coordination
        # text_to_mol(calc_coords=1) uses the robust template-growth engine
        o2_mol = oasa_smiles.text_to_mol(smiles_str, calc_coords=1)
        
        if not o2_mol or not o2_mol.vertices:
            raise ValueError("OASA failed to re-parse generated SMILES")
            
        # 4. Map coordinates back
        # OASA's SMILES parser creates heavy atom vertices in the EXACT same sequence
        # they appeared in the SMILES (and same sequence we found them in _processed_atoms).
        final_coords = {}
        # Filter O2 vertices for non-H if any leaked in (usually they don't from text_to_mol)
        o2_heavy_vertices = [v for v in o2_mol.vertices if v.symbol != 'H']
        
        for i, o2_v in enumerate(o2_heavy_vertices):
            if i < len(smiles_to_domain_map):
                domain_idx = smiles_to_domain_map[i]
                final_coords[domain_idx] = [float(o2_v.x), float(o2_v.y)]
        
        # Verify we got everything
        if len(final_coords) < len(atom_map):
            print(f"[DEBUG] Golden Layout incomplete: {len(final_coords)}/{len(atom_map)} - falling back.")
            raise ValueError("Mapping loss")

        return final_coords
        
    except Exception as e:
        print(f"[DEBUG] Golden Layout failed: {e}. Falling back to standard OASA.")
        # Final fallback - standard coordinate generation. The prep helper
        # was already applied above before the SMILES generation attempt,
        # so OASA here also sees a kekulized skeleton and produces a
        # clean, deterministic layout.
        oasa_cg.calculate_coords(o1_mol, bond_length=1.0, force=1)
        coords = {}
        for idx, v in atom_map.items():
            coords[idx] = [float(v.x), float(v.y)]
        return coords
