import math
import numpy as np
from typing import List, Set, Dict, Tuple, Optional
from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom

class InteractionType:
    H_BOND = "Hydrogen Bond"
    HYDROPHOBIC = "Hydrophobic"
    SALT_BRIDGE = "Salt Bridge"
    PI_STACKING = "Pi-Stacking"

class MolecularInteraction:
    def __init__(self, atom1_idx: int, atom2_idx: int, type: str, distance: float, 
                 energy: float = 0.0, details: str = ""):
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.type = type
        self.distance = distance
        self.energy = energy
        self.details = details

class DockingPoseService:
    def __init__(self, molecule: Molecule):
        self.molecule = molecule

    def find_ligands(self) -> List[Set[int]]:
        """Identify potential ligand fragments, prioritizing HETATMs."""
        fragments = self.molecule.get_fragments()
        
        def fragment_score(frag):
            # Higher score = more likely to be a ligand
            het_count = sum(1 for idx in frag if getattr(self.molecule.atoms[idx], 'is_hetatm', False))
            atom0 = self.molecule.atoms[next(iter(frag))]
            res_name = getattr(atom0, 'res_name', '').upper()
            is_water = res_name in ('HOH', 'WAT', 'SOL', 'DOD')
            
            if is_water: return -1 # Waters are last
            
            het_ratio = het_count / len(frag) if frag else 0
            # Prioritize small-to-medium fragments with high het ratio
            score = het_ratio * 100
            if 5 < len(frag) < 100: score += 50 
            return score

        # Return all fragments (except water) sorted by score
        scored_frags = []
        for frag in fragments:
            if fragment_score(frag) >= 0:
                scored_frags.append(set(frag))
                
        return sorted(scored_frags, key=fragment_score, reverse=True)

    def find_nearby_residues(self, ligand_indices: Set[int], distance: float = 5.0) -> Set[int]:
        """Find residue sequence numbers within distance of the ligand."""
        nearby_res = set()
        ligand_atoms = [self.molecule.atoms[i] for i in ligand_indices if self.molecule.atoms[i].has_coords]
        
        for atom in self.molecule.atoms:
            if atom.index in ligand_indices or not atom.has_coords:
                continue
            
            # Check distance to any ligand atom
            for l_atom in ligand_atoms:
                dx = atom.x - l_atom.x
                dy = atom.y - l_atom.y
                dz = (atom.z or 0) - (l_atom.z or 0)
                if dx*dx + dy*dy + dz*dz <= distance*distance:
                    if hasattr(atom, 'res_seq') and atom.res_seq is not None:
                        nearby_res.add(atom.res_seq)
                    break
        
        return nearby_res

    def detect_interactions(self, ligand_indices: Set[int], residue_indices: Set[int]) -> List[MolecularInteraction]:
        """Detect various molecular interactions between ligand and nearby residues."""
        interactions = []
        
        # 1. H-Bonds
        interactions.extend(self._detect_hbonds(ligand_indices, residue_indices))
        
        # 2. Salt Bridges
        interactions.extend(self._detect_salt_bridges(ligand_indices, residue_indices))
        
        # 3. Hydrophobic
        interactions.extend(self._detect_hydrophobic(ligand_indices, residue_indices))
        
        return interactions

    def _detect_hbonds(self, ligand_indices: Set[int], residue_indices: Set[int]) -> List[MolecularInteraction]:
        """Detect H-bonds between ligand and residues."""
        hbonds = []
        
        # Define donors and acceptors
        # Simplification: N, O, S are potential acceptors. 
        # N, O with attached H are donors.
        
        def is_acceptor(atom: Atom) -> bool:
            return atom.symbol in ('N', 'O', 'S', 'F')
        
        def is_donor(atom: Atom, mol: Molecule) -> bool:
            if atom.symbol not in ('N', 'O', 'S'):
                return False
            # Check for explicit H
            for nb_idx in mol.get_neighbors(atom.index):
                if mol.atoms[nb_idx].symbol == 'H':
                    return True
            # Check for implicit H
            if getattr(atom, 'num_implicit_h', 0) > 0:
                return True
            return False

        ligand_atoms = [self.molecule.atoms[i] for i in ligand_indices]
        protein_atoms = [a for a in self.molecule.atoms if a.res_seq in residue_indices]
        
        # Ligand as donor, protein as acceptor
        for l_idx in ligand_indices:
            l_atom = self.molecule.atoms[l_idx]
            if not is_donor(l_atom, self.molecule) or not l_atom.has_coords:
                continue
            
            for p_atom in protein_atoms:
                if not is_acceptor(p_atom) or not p_atom.has_coords:
                    continue
                
                dist = self._dist(l_atom, p_atom)
                if 2.5 <= dist <= 3.6:
                    hbonds.append(MolecularInteraction(l_idx, p_atom.index, InteractionType.H_BOND, dist))

        # Protein as donor, ligand as acceptor
        for p_atom in protein_atoms:
            if not is_donor(p_atom, self.molecule) or not p_atom.has_coords:
                continue
            
            for l_idx in ligand_indices:
                l_atom = self.molecule.atoms[l_idx]
                if not is_acceptor(l_atom) or not l_atom.has_coords:
                    continue
                
                dist = self._dist(p_atom, l_atom)
                if 2.5 <= dist <= 3.6:
                    hbonds.append(MolecularInteraction(p_atom.index, l_idx, InteractionType.H_BOND, dist))
                    
        return hbonds

    def _detect_salt_bridges(self, ligand_indices: Set[int], residue_indices: Set[int]) -> List[MolecularInteraction]:
        """Detect ionic interactions."""
        # Highly simplified: Oppositely charged heavy atoms within 4.0A
        bridges = []
        ligand_atoms = [self.molecule.atoms[i] for i in ligand_indices]
        protein_atoms = [a for a in self.molecule.atoms if a.res_seq in residue_indices]
        
        for l_atom in ligand_atoms:
            l_q = getattr(l_atom, 'formal_charge', 0)
            if l_q == 0: continue
            
            for p_atom in protein_atoms:
                p_q = getattr(p_atom, 'formal_charge', 0)
                if p_q == 0: continue
                
                if l_q * p_q < 0: # Opposite signs
                    dist = self._dist(l_atom, p_atom)
                    if dist <= 4.5:
                        bridges.append(MolecularInteraction(l_atom.index, p_atom.index, InteractionType.SALT_BRIDGE, dist))
        return bridges

    def _detect_hydrophobic(self, ligand_indices: Set[int], residue_indices: Set[int]) -> List[MolecularInteraction]:
        """Detect hydrophobic contacts."""
        # Non-polar carbons within 4.5A
        contacts = []
        
        def is_hydrophobic(atom: Atom) -> bool:
            if atom.symbol != 'C': return False
            # Carbon bonded only to C or H is hydrophobic
            for nb_idx in self.molecule.get_neighbors(atom.index):
                nb = self.molecule.atoms[nb_idx]
                if nb.symbol not in ('C', 'H'):
                    return False
            return True

        ligand_h = [i for i in ligand_indices if is_hydrophobic(self.molecule.atoms[i])]
        protein_h = [a.index for a in self.molecule.atoms if a.res_seq in residue_indices and is_hydrophobic(a)]
        
        for li in ligand_h:
            la = self.molecule.atoms[li]
            if not la.has_coords: continue
            for pi in protein_h:
                pa = self.molecule.atoms[pi]
                if not pa.has_coords: continue
                dist = self._dist(la, pa)
                if dist <= 4.5:
                    contacts.append(MolecularInteraction(li, pi, InteractionType.HYDROPHOBIC, dist))
        
        # Cull hydrophobic contacts to only keep the shortest one per residue-ligand atom pair if needed, 
        # but for now we just return all.
        return contacts

    def _dist(self, a1, a2):
        return math.sqrt((a1.x-a2.x)**2 + (a1.y-a2.y)**2 + ((a1.z or 0)-(a2.z or 0))**2)

    def generate_report_csv(self, interactions: List[MolecularInteraction]) -> str:
        """Generate a CSV string from interactions."""
        lines = ["Atom 1,Element 1,Residue 1,Atom 2,Element 2,Residue 2,Type,Distance (A)"]
        for inter in interactions:
            a1 = self.molecule.atoms[inter.atom1_idx]
            a2 = self.molecule.atoms[inter.atom2_idx]
            
            res1 = f"{getattr(a1, 'res_name', 'UNK')}{getattr(a1, 'res_seq', '')}"
            res2 = f"{getattr(a2, 'res_name', 'UNK')}{getattr(a2, 'res_seq', '')}"
            
            lines.append(f"{a1.pdb_name or a1.symbol}{a1.index},{a1.symbol},{res1},"
                         f"{a2.pdb_name or a2.symbol}{a2.index},{a2.symbol},{res2},"
                         f"{inter.type},{inter.distance:.3f}")
        return "\n".join(lines)
