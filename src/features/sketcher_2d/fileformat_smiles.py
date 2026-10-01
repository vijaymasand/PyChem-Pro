# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3) and OASA
from src.vendors.oasa.smiles import smiles as OasaSmiles
from .molecule import Molecule
from src.vendors.oasa.coords_generator import coords_generator as OasaCoordsGenerator

class Smiles:
    def __init__(self):
        self._oasa = OasaSmiles()

    def get_molecule(self, smiles_text):
        self._oasa.read_smiles(smiles_text)
        oasa_mol = self._oasa.structure
        if not oasa_mol: return None
        return oasa_mol

    def generate(self, mol):
        # Ensure the molecule is properly initialized for OASA
        if hasattr(mol, '_flush_cache'):
            mol._flush_cache()
        
        # Try to mark aromatic bonds if available
        if hasattr(mol, 'mark_aromatic_bonds'):
            try:
                mol.mark_aromatic_bonds()
            except:
                pass
        
        return self._oasa.get_smiles(mol)

class CoordsGenerator:
    def __init__(self):
        self._gen = OasaCoordsGenerator()

    def calculate_coords(self, mol, bond_length=30):
        self._gen.calculate_coords(mol)
        # OASA coords might be small, scale them
        for atom in mol.atoms:
            atom.x *= bond_length
            atom.y *= bond_length
