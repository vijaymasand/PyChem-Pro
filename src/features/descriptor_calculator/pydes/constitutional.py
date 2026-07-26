"""
PyDes Constitutional Descriptors

This module generates constitutional descriptors (elemental counts, ring counts, etc.).
It dynamically creates variants to provide a comprehensive set of descriptors.
"""
from typing import Dict, Any

def calculate_constitutional(molecule) -> Dict[str, Any]:
    """Calculate all constitutional descriptors for the given molecule."""
    results = {}
    
    # Pre-compute basic properties
    atoms = molecule.atoms
    bonds = molecule.bonds
    
    num_atoms = len(atoms)
    num_bonds = len(bonds)
    
    results["nAtoms"] = num_atoms
    results["nBonds"] = num_bonds
    
    # Element counts
    elements = ['H', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B', 'Si']
    counts = {el: 0 for el in elements}
    
    heavy_atoms = 0
    for atom in atoms:
        sym = atom.symbol
        if sym in counts:
            counts[sym] += 1
        if sym != 'H':
            heavy_atoms += 1
            
    results["nHeavyAtoms"] = heavy_atoms
    
    for el in elements:
        results[f"n{el}"] = counts[el]
        
    # Element ratios
    if counts['C'] > 0:
        for el in elements:
            if el != 'C':
                results[f"Ratio_{el}C"] = counts[el] / counts['C']
    else:
        for el in elements:
            if el != 'C':
                results[f"Ratio_{el}C"] = 0.0
                
    # Proportions
    if num_atoms > 0:
        for el in elements:
            results[f"Prop_{el}"] = counts[el] / num_atoms
            
    # Degree counts
    degrees = [0] * num_atoms
    for bond in bonds:
        degrees[bond.begin_atom_idx] += 1
        degrees[bond.end_atom_idx] += 1
        
    for i in range(1, 7): # Degrees 1 to 6
        results[f"nDegree{i}"] = sum(1 for d in degrees if d == i)
        
    # Bond types
    single = 0
    double = 0
    triple = 0
    aromatic = 0
    
    for bond in bonds:
        if bond.is_single: single += 1
        elif bond.is_double: double += 1
        elif bond.is_triple: triple += 1
        
        if bond.is_aromatic: aromatic += 1
        
    results["nSingleBonds"] = single
    results["nDoubleBonds"] = double
    results["nTripleBonds"] = triple
    results["nAromaticBonds"] = aromatic
    
    # always emit the proportions, a molecule without bonds simply has zeros;
    # skipping the keys would punch NaN holes into a batch CSV
    denominator = num_bonds if num_bonds > 0 else 1
    results["Prop_SingleBonds"] = single / denominator
    results["Prop_DoubleBonds"] = double / denominator
    results["Prop_TripleBonds"] = triple / denominator
    results["Prop_AromaticBonds"] = aromatic / denominator
    
    # Ring counts (simple)
    if hasattr(molecule, 'find_rings'):
        rings = molecule.find_rings()
        results["nRings"] = len(rings)
        
        for size in range(3, 13):
            results[f"nRing{size}"] = sum(1 for r in rings if len(r) == size)
    else:
        results["nRings"] = 0
        for size in range(3, 13):
            results[f"nRing{size}"] = 0
            
    return results
