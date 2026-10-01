"""
Lipophilicity Service — Calculate Molecular Lipophilicity Potential (MLP) based logP.

Provides a robust method for logP calculation based on fragmental values and 
surface-integrated potentials, ported from the MLP_Tools implementation.
"""

import math
import numpy as np
from typing import List, Tuple, Dict, Set, Optional
from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom
from src.core.domain.models.bond import BondType
from .lipophilicity_data import FRAGMENT_CORE, FRAGMENT_FB

# Constants from MLP_Tools (libmlpy2.py and mlp_main2.py)
LOGP_CONSTANTS = {
    'N_SPTS': 2000,
    'MINUS': 0.00072402,
    'PLUS': 0.00125992,
    'CONSTANT': -0.12439948,
    'F_CONST': 2.0,
    'F_CUTOFF': 3.25,
    'F_DELTA': 1.33,
    'D_CUTOFF': 20.0
}

# vdW Radii from MLP_Tools (asa_radii in libmlpy2.py)
ASA_RADII = {
    'AC': 2.0, 'AL': 2.0, 'AM': 2.0, 'SB': 2.0, 'AR': 1.88, 'AS': 1.85, 'AT': 2.0,
    'BA': 2.0, 'BK': 2.0, 'BE': 2.0, 'BI': 2.0, 'BH': 2.0, 'B': 2.0, 'BR': 1.85,
    'CD': 1.58, 'CS': 2.0, 'CA': 2.0, 'CF': 2.0, 'C': 1.7, 'CE': 2.0, 'CL': 1.75,
    'L': 1.75, 'CR': 2.0, 'CO': 2.0, 'CU': 1.4, 'CM': 2.0, 'DS': 2.0, 'DB': 2.0,
    'DY': 2.0, 'ES': 2.0, 'ER': 2.0, 'EU': 2.0, 'FM': 2.0, 'F': 1.47, 'FR': 2.0,
    'GD': 2.0, 'GA': 1.87, 'GE': 2.0, 'AU': 1.66, 'HF': 2.0, 'HS': 2.0, 'HE': 1.4,
    'HO': 2.0, 'H': 1.09, 'IN': 1.93, 'I': 1.98, 'IR': 2.0, 'FE': 2.0, 'KR': 2.02,
    'LA': 2.0, 'LR': 2.0, 'PB': 2.02, 'LI': 1.82, 'LU': 2.0, 'MG': 1.73, 'MN': 2.0,
    'MT': 2.0, 'MD': 2.0, 'HG': 1.55, 'MO': 2.0, 'ND': 2.0, 'NE': 1.54, 'NP': 2.0,
    'NI': 1.63, 'NB': 2.0, 'N': 1.55, 'NO': 2.0, 'OS': 2.0, 'O': 1.52, 'PD': 1.63,
    'P': 1.8, 'PT': 1.72, 'PU': 2.0, 'PO': 2.0, 'K': 2.75, 'PR': 2.0, 'PM': 2.0,
    'PA': 2.0, 'RA': 2.0, 'RN': 2.0, 'RE': 2.0, 'RH': 2.0, 'RB': 2.0, 'RU': 2.0,
    'RF': 2.0, 'SM': 2.0, 'SC': 2.0, 'SG': 2.0, 'SE': 1.9, 'SI': 2.1, 'AG': 1.72,
    'NA': 2.27, 'SR': 2.0, 'S': 1.8, 'TA': 2.0, 'TC': 2.0, 'TE': 2.06, 'TB': 2.0,
    'TL': 1.96, 'TH': 2.0, 'TM': 2.0, 'SN': 2.17, 'TI': 2.0, 'W': 2.0, 'U': 1.86,
    'V': 2.0, 'XE': 2.16, 'YB': 2.0, 'Y': 2.0, 'ZN': 1.39, 'ZR': 2.0, 'DU': 0.0,
    'LP': 0.85
}

class LipophilicityService:
    """
    Service for calculating Molecular Lipophilicity Potential (MLP) based logP.
    Re-implemented from MLP_Tools PyMOL plugin logic.
    """
    
    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        
        # Ensure aromaticity is perceived before typing
        from src.features.smiles_parser.rules.aromaticity import perceive_aromaticity
        
        # Always perceive aromaticity to ensure atoms and bonds are synced
        # and rings are correctly identified. This is fast and avoids bugs
        # where bonds are marked aromatic but atoms are not.
        perceive_aromaticity(molecule)
        
        # Propagate from bonds to atoms just in case (for MOL2/SDF type 4)
        molecule.propagate_aromaticity()
        
        self.molecule.assign_sybyl_types()
        self.atom_fragment_values = [0.0] * len(molecule.atoms)

    def calculate_logp(self) -> float:
        """
        Calculate logP using the fragment-based MLP method.
        
        Returns:
            float: Calculated logP value.
        """
        contributions = self.get_atomic_contributions()
        if not contributions:
            return 0.0
            
        return sum(contributions) + LOGP_CONSTANTS['CONSTANT']

    def get_atomic_contributions(self) -> List[float]:
        """
        Calculate atomic logP contributions.
        
        Returns:
            List[float]: LogP contribution for each atom.
        """
        if not self.molecule.atoms:
            return []
            
        # 1. Assign fragmental values to atoms
        self._assign_fragment_values()
        
        # 2. Generate Solvent Accessible Surface (ASA) points
        points = self._generate_asa_points()
        if len(points) == 0:
            return [0.0] * len(self.molecule.atoms)
            
        # 3. Calculate MLP at each surface point
        # We need the individual potentials for each atom at each point
        n_atoms = len(self.molecule.atoms)
        atom_coords = np.array([[a.x, a.y, a.z] for a in self.molecule.atoms])
        frag_vals = np.array(self.atom_fragment_values)
        
        f_const = LOGP_CONSTANTS['F_CONST']
        f_cutoff = LOGP_CONSTANTS['F_CUTOFF']
        f_delta = LOGP_CONSTANTS['F_DELTA']
        plus_factor = LOGP_CONSTANTS['PLUS']
        minus_factor = LOGP_CONSTANTS['MINUS']
        d_cutoff = LOGP_CONSTANTS['D_CUTOFF']
        
        # We'll accumulate contributions per atom
        atomic_contributions = np.zeros(n_atoms)
        
        chunk_size = 1000
        for i in range(0, len(points), chunk_size):
            end = min(i + chunk_size, len(points))
            chunk_points = points[i:end]
            
            # dists: (chunk_points, n_atoms)
            diffs = chunk_points[:, np.newaxis, :] - atom_coords[np.newaxis, :, :]
            dists = np.sqrt(np.sum(diffs**2, axis=2))
            
            # Fermi-type distance function
            f_num = 1 + np.exp((-f_const * f_cutoff) / f_delta)
            potentials = f_num / (1 + np.exp((f_const * (dists - f_cutoff)) / f_delta))
            potentials[dists > d_cutoff] = 0.0
            
            # mlp_p: (chunk_points,)
            mlp_p = np.sum(potentials * frag_vals, axis=1)
            
            # factor_p: (chunk_points,)
            factors = np.where(mlp_p > 0, plus_factor, minus_factor)
            
            # Contribution per point per atom: factor_p * potential_pi * frag_val_i
            # We want to sum over p: sum_p (factors * potentials) * frag_vals
            # (factors[:, np.newaxis] * potentials): (chunk_points, n_atoms)
            # sum over axis 0: (n_atoms,)
            chunk_contribs = np.sum(factors[:, np.newaxis] * potentials, axis=0) * frag_vals
            atomic_contributions += chunk_contribs
            
        return atomic_contributions.tolist()

    def _assign_fragment_values(self):
        """
        Assign fragmental logP values to each atom based on its environment.
        Uses the fingerprinting method from MLP_Tools.
        """
        atoms = self.molecule.atoms
        n_atoms = len(atoms)
        
        # Initialize fingerprints and neighbor data
        # fp_group1: 4 digits for atom's own bonds
        # fp_groups_2_5: 4 sets of 4 digits for neighbors
        fp_group1 = [[] for _ in range(n_atoms)]
        fp_groups_2_5 = [[] for _ in range(n_atoms)]
        
        ignore_atoms = {'H', 'DU', 'LP'}
        
        def is_carbon(sym):
            return sym.upper() == 'C'
            
        # First pass: Build Group 1 of fingerprint and handle heteroatom-heteroatom connections
        for bond in self.molecule.bonds:
            id1, id2 = bond.begin_atom_idx, bond.end_atom_idx
            a1, a2 = atoms[id1], atoms[id2]
            
            if a1.symbol.upper() in ignore_atoms or a2.symbol.upper() in ignore_atoms:
                continue
                
            # Bond type mapping for first group
            btype_val = 1 # Default single
            if bond.is_double: btype_val = 2
            elif bond.is_triple: btype_val = 3
            elif bond.is_aromatic: btype_val = 5
            elif bond.is_amide: btype_val = 1
            
            fp_group1[id1].append(btype_val)
            fp_group1[id2].append(btype_val)
            
            # Remaining groups for heteroatom-heteroatom connections
            if not is_carbon(a1.symbol) and not is_carbon(a2.symbol):
                # Bond type mapping for hetero-hetero
                hh_val = 6 # Default single
                if bond.is_double: hh_val = 7
                elif bond.is_triple: hh_val = 9
                elif bond.is_aromatic: hh_val = 8
                elif bond.is_amide: hh_val = 6
                
                fp_groups_2_5[id1].append([hh_val, 0, 0, 0])
                fp_groups_2_5[id2].append([hh_val, 0, 0, 0])

        # Normalize Group 1: pad with zeros, sort descending
        for i in range(n_atoms):
            while len(fp_group1[i]) < 4:
                fp_group1[i].append(0)
            fp_group1[i].sort(reverse=True)
            
        # Second pass: Build Group 2-5 for heteroatom-carbon connections
        for bond in self.molecule.bonds:
            id1, id2 = bond.begin_atom_idx, bond.end_atom_idx
            a1, a2 = atoms[id1], atoms[id2]
            
            if a1.symbol.upper() in ignore_atoms or a2.symbol.upper() in ignore_atoms:
                continue
                
            if is_carbon(a1.symbol) and not is_carbon(a2.symbol):
                fp_groups_2_5[id2].append(fp_group1[id1])
            if is_carbon(a2.symbol) and not is_carbon(a1.symbol):
                fp_groups_2_5[id1].append(fp_group1[id2])

        # Final fingerprint construction and lookup
        for i in range(n_atoms):
            atom = atoms[i]
            if atom.symbol.upper() in ignore_atoms:
                continue
                
            # Normalize Group 2-5: pad with [0,0,0,0], sort descending
            while len(fp_groups_2_5[i]) < 4:
                fp_groups_2_5[i].append([0, 0, 0, 0])
            fp_groups_2_5[i].sort(key=lambda x: "".join(map(str, x)), reverse=True)
            
            # Combine into 20-digit string
            fp_str = "".join(map(str, fp_group1[i]))
            for group in fp_groups_2_5[i]:
                fp_str += "".join(map(str, group))
                
            # Lookup in FRAGMENT_CORE
            val = FRAGMENT_CORE.get((atom.symbol, fp_str))
            
            # Fallback to FRAGMENT_FB (short fingerprint)
            if val is None:
                short_fp = "".join(map(str, fp_group1[i]))
                val = FRAGMENT_FB.get((atom.symbol, short_fp), 0.0)
                
            self.atom_fragment_values[i] = val

    def _generate_asa_points(self, probe: float = 1.4) -> np.ndarray:
        """
        Generate surface points for Solvent Accessible Surface Area.
        Uses Fibonacci sphere algorithm and steric occlusion checking.
        """
        n_points = LOGP_CONSTANTS['N_SPTS']
        all_points = []
        
        # Pre-generate Fibonacci sphere points (unit sphere)
        indices = np.arange(0, n_points, dtype=float) + 0.5
        phi = np.arccos(1 - 2*indices/n_points)
        theta = np.pi * (1 + 5**0.5) * indices
        
        unit_points = np.stack([
            np.cos(theta) * np.sin(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(phi)
        ], axis=1)
        
        atom_coords = np.array([[a.x, a.y, a.z] for a in self.molecule.atoms])
        atom_radii = np.array([ASA_RADII.get(a.symbol.upper(), 1.7) + probe 
                             for a in self.molecule.atoms])
        
        for i in range(len(self.molecule.atoms)):
            atom_pos = atom_coords[i]
            radius = atom_radii[i]
            
            # Scale unit sphere to atom's SAS radius
            sphere_points = unit_points * radius + atom_pos
            
            # Vectorized occlusion check
            # diffs: (n_points, n_atoms, 3)
            diffs = sphere_points[:, np.newaxis, :] - atom_coords[np.newaxis, :, :]
            dists_sq = np.sum(diffs**2, axis=2)
            
            # Mask out current atom from distance check
            dists_sq[:, i] = np.inf
            
            # Check if point is outside all other atoms' SAS radii
            is_exposed = np.all(dists_sq >= atom_radii**2 - 1e-6, axis=1)
            
            exposed_points = sphere_points[is_exposed]
            if len(exposed_points) > 0:
                all_points.append(exposed_points)
            
        if not all_points:
            return np.array([])
            
        return np.concatenate(all_points)

    def _calculate_mlp(self, points: np.ndarray) -> np.ndarray:
        """
        Calculate Molecular Lipophilicity Potential at each surface point.
        Uses Fermi-type distance function potential summation.
        """
        if len(points) == 0:
            return np.array([])
            
        n_atoms = len(self.molecule.atoms)
        atom_coords = np.array([[a.x, a.y, a.z] for a in self.molecule.atoms])
        frag_vals = np.array(self.atom_fragment_values)
        
        # Calculate point-atom distances
        # Optimization: use chunks if points are too many to avoid OOM
        chunk_size = 1000
        mlp_values = np.zeros(len(points))
        
        f_const = LOGP_CONSTANTS['F_CONST']
        f_cutoff = LOGP_CONSTANTS['F_CUTOFF']
        f_delta = LOGP_CONSTANTS['F_DELTA']
        
        for i in range(0, len(points), chunk_size):
            end = min(i + chunk_size, len(points))
            chunk_points = points[i:end]
            
            # diffs: (chunk_points, n_atoms, 3)
            diffs = chunk_points[:, np.newaxis, :] - atom_coords[np.newaxis, :, :]
            dists = np.sqrt(np.sum(diffs**2, axis=2))
            
            # Fermi-type distance function
            f_num = 1 + np.exp((-f_const * f_cutoff) / f_delta)
            potentials = f_num / (1 + np.exp((f_const * (dists - f_cutoff)) / f_delta))
            
            # Apply distance cutoff (optimization/consistency from MLP_Tools)
            d_cutoff = LOGP_CONSTANTS['D_CUTOFF']
            potentials[dists > d_cutoff] = 0.0
            
            # Weighted potentials sum
            mlp_values[i:end] = np.sum(potentials * frag_vals, axis=1)
            
        return mlp_values


def calculate_logp(molecule: Molecule) -> float:
    """
    Convenience function to calculate logP for a molecule.
    """
    service = LipophilicityService(molecule)
    return service.calculate_logp()
