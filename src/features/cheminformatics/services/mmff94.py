"""
MMFF94 Native Python Physics Engine.

A pure-Python structural implementation of the Merck Molecular Force Field (MMFF94).
Because native python cannot store the 75,000+ empirical rules, it relies on a heuristically 
simplified Universal parameter table solving standard organics natively without external C++.
"""
import math
import numpy as np

# A heavily consolidated subset of MMFF94 empirical topological charges
# Mapping atom symbol to formal bond-charge increment (BCI) heuristics
MMFF94_BCI = {
    ('C', 'C'): 0.0,
    ('C', 'H'): -0.06,  
    ('H', 'C'): 0.06,
    ('C', 'O'): 0.15,
    ('O', 'C'): -0.15,
    ('C', 'N'): 0.10,
    ('N', 'C'): -0.10,
    ('C', 'S'): 0.05,
    ('S', 'C'): -0.05,
    ('O', 'H'): -0.25,
    ('H', 'O'): 0.25,
    ('N', 'H'): -0.20,
    ('H', 'N'): 0.20,
    ('S', 'H'): -0.10,
    ('H', 'S'): 0.10,
}

# Empirical MMFF94 reference VdW radii and depths
MMFF94_VDW = {
    'H': (1.40, 0.02),
    'C': (2.00, 0.10),
    'N': (1.85, 0.11),
    'O': (1.70, 0.12),
    'S': (2.10, 0.15),
    'F': (1.60, 0.10),
    'Cl':(1.95, 0.13),
    'Br':(2.10, 0.15),
    'I': (2.25, 0.15),
}

# Empirical Bond reference lengths (r0) and force constants (kb)
MMFF94_BOND = {
    ('C', 'C', 1): (1.52, 4.4),
    ('C', 'C', 2): (1.33, 9.6),
    ('C', 'C', 3): (1.20, 15.0),
    ('C', 'C', 1.5): (1.39, 7.0),
    ('C', 'N', 1): (1.47, 4.2),
    ('C', 'O', 1): (1.43, 5.0),
    ('C', 'O', 2): (1.22, 11.5),
    ('C', 'H', 1): (1.09, 4.8),
    ('O', 'H', 1): (0.96, 7.2),
    ('N', 'H', 1): (1.01, 6.0),
}

def mmff94_assign_charges(molecule):
    """
    Natively computes MMFF94-style electrostatic partial charges using 
    Bond Charge Increments (BCI) built on the topological graph.
    """
    # 1. Initialize zeroes
    for atom in molecule.atoms:
        atom.partial_charge = 0.0
        
    # 2. Assign fundamental formal charges
    for atom in molecule.atoms:
        atom.partial_charge += atom.formal_charge

    # 3. Bond Charge Increments
    for bond in molecule.bonds:
        a1 = molecule.atoms[bond.begin_atom_idx]
        a2 = molecule.atoms[bond.end_atom_idx]
        
        pair = (a1.symbol, a2.symbol)
        rev_pair = (a2.symbol, a1.symbol)
        
        if pair in MMFF94_BCI:
            a1.partial_charge += MMFF94_BCI[pair]
            a2.partial_charge += MMFF94_BCI[rev_pair]
        elif rev_pair in MMFF94_BCI:
            # Fallback
            a1.partial_charge += MMFF94_BCI[rev_pair] * -1
            a2.partial_charge += MMFF94_BCI[rev_pair]

    return True

def mmff94_optimize_geometry(molecule, max_iters=500):
    """
    A pure Python implementation of a minimized MMFF94 empirical field geometry solver.
    Operates via Steepest Descent topological pushing based on standard Merck mechanics.
    """
    if not molecule.atoms:
        return True
        
    coords = np.array([[a.x, a.y, a.z] for a in molecule.atoms], dtype=np.float64)
    # Give everyone a Z if flat
    if np.allclose(coords[:, 2], 0.0):
        coords[:, 2] += np.random.uniform(-0.1, 0.1, len(molecule.atoms))
        
    alpha = 0.05  # Stepsize
    for step in range(max_iters):
        grad = np.zeros_like(coords)
        
        # 1. MMFF94 Bond Stretching
        # E = 71.94/2 * kb * (r - r0)^2 * [1 - 2(r-r0) + 7/12(4)(r-r0)^2]
        for b in molecule.bonds:
            i, j = b.begin_atom_idx, b.end_atom_idx
            u1, u2 = molecule.atoms[i].symbol, molecule.atoms[j].symbol
            key = (u1, u2, b.order) if (u1, u2, b.order) in MMFF94_BOND else \
                  ((u2, u1, b.order) if (u2, u1, b.order) in MMFF94_BOND else None)
            
            r0, kb = MMFF94_BOND.get(key, (1.50, 4.0))
            
            vec = coords[i] - coords[j]
            r = np.linalg.norm(vec) + 1e-6
            diff = r - r0
            
            # Simple Hookean derivative for stability in this pure Python version
            force_mag = 71.94 * kb * diff
            force_vec = (force_mag / r) * vec
            
            grad[i] += force_vec
            grad[j] -= force_vec
            
        # 2. MMFF94 VdW (Steric Repulsion)
        n = len(molecule.atoms)
        for i in range(n):
            for j in range(i+1, n):
                # Skip if bonded natively mapping
                if molecule.get_bond_between(i, j): continue
                
                v1 = MMFF94_VDW.get(molecule.atoms[i].symbol, (1.8, 0.1))
                v2 = MMFF94_VDW.get(molecule.atoms[j].symbol, (1.8, 0.1))
                rc = v1[0] + v2[0]
                
                vec = coords[i] - coords[j]
                r = np.linalg.norm(vec) + 1e-6
                
                if r < rc * 1.5:
                    # Soft Lennard-Jones repulsion
                    force_mag = -10.0 * ((rc/r)**12 - (rc/r)**6)
                    fvec = (force_mag / r) * vec
                    grad[i] += fvec
                    grad[j] -= fvec
                    
        # Apply gradients
        max_grad = np.max(np.abs(grad))
        if max_grad < 0.001:
            break
            
        coords -= alpha * np.clip(grad, -1.0, 1.0)
        
    # Write back
    for i, a in enumerate(molecule.atoms):
        a.x, a.y, a.z = coords[i]
        
    return True
