"""
SMILES Generator — Converts molecules to SMILES strings using BKChem (OASA).

Delegates purely and fully to OASA's highly robust deterministic tree-growth
SMILES export algorithm, maintaining perfect canonization and
proper ring-closure tracking natively.
"""

from src.vendors.oasa_bridge import domain_to_oasa_mol
import src.vendors.oasa.smiles as oasa_smiles

class SMILESGenerator:
    """Generate SMILES strings from molecules using BKChem algorithms natively."""
    
    def __init__(self, molecule):
        self.molecule = molecule
        
    def generate(self):
        """Generate SMILES string completely seamlessly from a domain molecule."""
        if not self.molecule or not self.molecule.atoms:
            return ""
            
        try:
            # 1. Translate internal graph representation to exact OASA molecule
            o_mol, _ = domain_to_oasa_mol(self.molecule)
            
            # 2. Use OASA's smiles spanning conversion algorithm
            # (mol_to_text automatically calculates canonical tree-spanning closures natively)
            smiles = oasa_smiles.mol_to_text(o_mol)
            return smiles
        except Exception as e:
            print(f"[ERROR] Failed to generate SMILES natively via BKChem/OASA: {e}")
            return ""


def generate_smiles(molecule):
    """Generate SMILES string from a molecule."""
    generator = SMILESGenerator(molecule)
    return generator.generate()
