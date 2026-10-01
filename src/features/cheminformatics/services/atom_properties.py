"""
Atom property definitions for molecular analysis.

Provides definitions for H-bond donors, acceptors, lipophilic atoms,
and selection/counting functionality.
"""

from typing import List, Set, Tuple
from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom

class AtomPropertyAnalyzer:
    """
    Analyzes atom properties for molecular selection and analysis.
    
    Provides methods to identify H-bond donors, acceptors, lipophilic atoms,
    and count atoms by various properties.
    """
    
    # H-bond donor definitions
    H_BOND_DONORS = {
        # H attached to N, O, S can donate H-bonds
        'N',  # Amines, amides, etc.
        'O',  # Alcohols, phenols, carboxylic acids
        'S',  # Thiols, sulfides
    }
    
    # H-bond acceptor definitions  
    H_BOND_ACCEPTORS = {
        'N',  # Amines, amides, nitriles, pyridines
        'O',  # Carbonyls, alcohols, ethers, carboxylates
        'S',  # Thiocarbonyls, sulfones
        'F',  # Fluorine (weak acceptor)
        'Cl', # Chlorine (weak acceptor)
        'Br', # Bromine (weak acceptor)
    }
    
    # Lipophilic atom definitions (hydrophobic)
    LIPOPHILIC_ATOMS = {
        'H',  # Hydrogen (when attached to non-polar atoms)
        'C',  # Carbon atoms (especially sp2/sp3)
        'S',  # Sulfur (except in polar groups)
        'P',  # Phosphorus
        'Si', # Silicon
        'Cl', # Halogens (moderately lipophilic)
        'Br', # Halogens (moderately lipophilic)
        'I',  # Iodine (lipophilic)
    }
    
    # Polar atom definitions (opposite of lipophilic)
    POLAR_ATOMS = {
        'N',  # Nitrogen
        'O',  # Oxygen
        'F',  # Fluorine
    }
    
    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        
    def is_hbond_donor(self, atom_index: int) -> bool:
        """
        Check if atom is a H-bond donor.
        
        H-bond donors are H atoms attached to electronegative atoms (N, O, S).
        """
        atom = self.molecule.atoms[atom_index]
        
        # Only hydrogen can be donors
        if atom.symbol != 'H':
            return False
            
        # Check if H is attached to donor atom
        for neighbor_idx, _ in self.molecule._adjacency.get(atom_index, []):
            neighbor = self.molecule.atoms[neighbor_idx]
            if neighbor.symbol in self.H_BOND_DONORS:
                return True
                
        return False
    
    def is_hbond_acceptor(self, atom_index: int) -> bool:
        """
        Check if atom is a H-bond acceptor.
        
        H-bond acceptors are electronegative atoms with lone pairs.
        """
        atom = self.molecule.atoms[atom_index]
        return atom.symbol in self.H_BOND_ACCEPTORS
    
    def is_lipophilic(self, atom_index: int) -> bool:
        """
        Check if atom is lipophilic (hydrophobic).
        
        Note: Carbon atoms attached directly to polar atoms are excluded from being lipophilic.
        Note: Hydrogen atoms attached directly to polar atoms are excluded from being lipophilic.
        """
        atom = self.molecule.atoms[atom_index]
        
        # Special cases for carbon - check polar neighbors first
        if atom.symbol == 'C':
            # Check if carbon is attached to any polar atoms
            # If attached to polar atoms, exclude from lipophilic
            for neighbor_idx in self.molecule.get_neighbors(atom_index):
                neighbor = self.molecule.atoms[neighbor_idx]
                if neighbor.symbol in self.POLAR_ATOMS:
                    return False  # Carbon attached to polar atom - not lipophilic
            
            # If not attached to polar atoms, carbon is lipophilic
            return True
            
        # Special cases for hydrogen - check only direct polar attachment
        if atom.symbol == 'H':
            # Check if hydrogen is attached to any polar atoms
            for neighbor_idx in self.molecule.get_neighbors(atom_index):
                neighbor = self.molecule.atoms[neighbor_idx]
                if neighbor.symbol in self.POLAR_ATOMS:
                    return False  # Hydrogen attached to polar atom - not lipophilic
            
            # If not attached to polar atoms, hydrogen is lipophilic
            return True
            
        # Basic check by element (for non-carbon, non-hydrogen atoms)
        if atom.symbol in self.LIPOPHILIC_ATOMS:
            return True
            
        return False
    
    def is_polar(self, atom_index: int) -> bool:
        """
        Check if atom is polar.
        """
        atom = self.molecule.atoms[atom_index]
        return atom.symbol in self.POLAR_ATOMS
    
    def get_atoms_by_property(self, property_name: str) -> List[int]:
        """
        Get list of atom indices that match a property.
        
        Args:
            property_name: 'donor', 'acceptor', 'lipophilic', 'polar'
            
        Returns:
            List of atom indices matching the property
        """
        atoms = []
        
        for i in range(len(self.molecule.atoms)):
            if property_name == 'donor' and self.is_hbond_donor(i):
                atoms.append(i)
            elif property_name == 'acceptor' and self.is_hbond_acceptor(i):
                atoms.append(i)
            elif property_name == 'lipophilic' and self.is_lipophilic(i):
                atoms.append(i)
            elif property_name == 'polar' and self.is_polar(i):
                atoms.append(i)
                
        return atoms
    
    def count_atoms_by_property(self, property_name: str) -> int:
        """
        Count atoms matching a property.
        
        Args:
            property_name: 'donor', 'acceptor', 'lipophilic', 'polar'
            
        Returns:
            Number of atoms matching the property
        """
        return len(self.get_atoms_by_property(property_name))
    
    def get_atom_properties_summary(self) -> dict:
        """
        Get summary of all atom properties in the molecule.
        
        Returns:
            Dictionary with counts for each property
        """
        return {
            'total_atoms': len(self.molecule.atoms),
            'hbond_donors': self.count_atoms_by_property('donor'),
            'hbond_acceptors': self.count_atoms_by_property('acceptor'),
            'lipophilic_atoms': self.count_atoms_by_property('lipophilic'),
            'polar_atoms': self.count_atoms_by_property('polar'),
            'donor_indices': self.get_atoms_by_property('donor'),
            'acceptor_indices': self.get_atoms_by_property('acceptor'),
            'lipophilic_indices': self.get_atoms_by_property('lipophilic'),
            'polar_indices': self.get_atoms_by_property('polar'),
        }
    
    def select_atoms_by_property(self, property_name: str) -> Set[int]:
        """
        Select atoms by property for GUI selection.
        
        Args:
            property_name: 'donor', 'acceptor', 'lipophilic', 'polar'
            
        Returns:
            Set of atom indices to be selected
        """
        return set(self.get_atoms_by_property(property_name))
    
    def get_atom_type_info(self, atom_index: int) -> dict:
        """
        Get detailed property information for a specific atom.
        
        Args:
            atom_index: Index of the atom
            
        Returns:
            Dictionary with all properties for the atom
        """
        atom = self.molecule.atoms[atom_index]
        
        return {
            'index': atom_index,
            'symbol': atom.symbol,
            'element': atom.element,
            'is_donor': self.is_hbond_donor(atom_index),
            'is_acceptor': self.is_hbond_acceptor(atom_index),
            'is_lipophilic': self.is_lipophilic(atom_index),
            'is_polar': self.is_polar(atom_index),
            'partial_charge': getattr(atom, 'partial_charge', None),
            'coordinates': (atom.x, atom.y, atom.z),
        }
    
    def find_hbond_pairs(self) -> List[Tuple[int, int]]:
        """
        Find potential H-bond donor-acceptor pairs.
        
        Returns:
            List of (donor_index, acceptor_index) tuples
        """
        pairs = []
        donors = self.get_atoms_by_property('donor')
        acceptors = self.get_atoms_by_property('acceptor')
        
        for donor_idx in donors:
            for acceptor_idx in acceptors:
                # Skip if same atom
                if donor_idx == acceptor_idx:
                    continue
                    
                # Check distance cutoff (simplified - should use 3D distance)
                donor_atom = self.molecule.atoms[donor_idx]
                acceptor_atom = self.molecule.atoms[acceptor_idx]
                
                distance = ((donor_atom.x - acceptor_atom.x)**2 + 
                           (donor_atom.y - acceptor_atom.y)**2 + 
                           (donor_atom.z - acceptor_atom.z)**2)**0.5
                
                # H-bond distance cutoff (3.5 Å is typical)
                if distance < 3.5:
                    pairs.append((donor_idx, acceptor_idx))
                    
        return pairs


def analyze_molecule_properties(molecule: Molecule) -> dict:
    """
    Convenience function to analyze all properties of a molecule.
    
    Args:
        molecule: Molecule object to analyze
        
    Returns:
        Dictionary with complete property analysis
    """
    analyzer = AtomPropertyAnalyzer(molecule)
    return analyzer.get_atom_properties_summary()


def select_atoms_by_property(molecule: Molecule, property_name: str) -> Set[int]:
    """
    Convenience function to select atoms by property.
    
    Args:
        molecule: Molecule object
        property_name: 'donor', 'acceptor', 'lipophilic', 'polar'
        
    Returns:
        Set of atom indices
    """
    analyzer = AtomPropertyAnalyzer(molecule)
    return analyzer.select_atoms_by_property(property_name)
