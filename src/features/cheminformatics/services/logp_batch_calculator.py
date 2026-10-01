"""
LogP Batch Calculator — Standalone utility to calculate MLP LogP for multiple chemical files.
Supports SDF, MOL, and MOL2 formats. Automatically generates 3D coordinates if missing.
"""

import os
import sys
import argparse
import pandas as pd
from typing import List, Optional

# Ensure project root is in path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../")))

from src.features.io.loaders.file_reader import read_mol, read_sdf, read_mol2
from src.features.smiles_parser.services.parser import parse_smiles
from src.features.cheminformatics.services.lipophilicity_service import calculate_logp
from src.services.coordinates.coord_gen_service import CoordinateGeneratorService
from src.core.domain.models.molecule import Molecule

def get_molecules_from_sdf(filepath: str) -> List[Molecule]:
    """Helper to read all molecules from an SDF file."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    blocks = content.split('$$$$')
    molecules = []
    for block in blocks:
        if not block.strip():
            continue
        try:
            mol = read_mol(block.strip(), is_string=True)
            molecules.append(mol)
        except Exception as e:
            print(f"Warning: Could not parse a molecule block in {filepath}: {e}")
    return molecules

def get_molecules_from_smiles(filepath: str) -> List[Molecule]:
    """Helper to read molecules from a SMILES file (one per line)."""
    molecules = []
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            # Handle "SMILES Name" format
            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else f"{os.path.basename(filepath)}_{i}"
            try:
                mol = parse_smiles(smiles)
                mol.name = name
                molecules.append(mol)
            except Exception as e:
                print(f"Warning: Could not parse SMILES on line {i+1}: {e}")
    return molecules

def process_files(input_paths: List[str], output_csv: Optional[str] = None):
    """
    Process a list of chemical files and calculate LogP for each molecule.
    
    Args:
        input_paths: List of file paths (MOL, SDF, MOL2)
        output_csv: Path to save results as CSV
    """
    results = []
    coord_svc = CoordinateGeneratorService()
    
    for path in input_paths:
        if not os.path.exists(path):
            print(f"Error: File not found: {path}")
            continue
            
        ext = os.path.splitext(path)[1].lower()
        mols = []
        
        try:
            if ext == '.mol':
                mols.append(read_mol(path))
            elif ext == '.sdf':
                mols.extend(get_molecules_from_sdf(path))
            elif ext == '.mol2':
                mols.append(read_mol2(path))
            elif ext in ('.smi', '.smiles', '.txt'):
                mols.extend(get_molecules_from_smiles(path))
            else:
                print(f"Skipping unsupported file format: {ext} ({path})")
                continue
        except Exception as e:
            print(f"Error reading {path}: {e}")
            continue
            
        for i, mol in enumerate(mols):
            mol_id = mol.name or f"{os.path.basename(path)}_{i}"
            print(f"Processing {mol_id}...")
            
            # Ensure 3D coordinates for MLP surface calculation
            # A molecule is considered 2D if all Z-coordinates are zero
            has_3d = any(abs(a.z) > 1e-4 for a in mol.atoms)
            if not has_3d:
                print(f"  2D layout detected. Generating 3D coordinates for {mol_id}...")
                coord_svc.generate_3d(mol)
                
            try:
                # Use Service directly for verbose access
                from src.features.cheminformatics.services.lipophilicity_service import LipophilicityService
                service = LipophilicityService(mol)
                logp = service.calculate_logp()
                
                results.append({
                    'File': os.path.basename(path),
                    'Molecule': mol_id,
                    'Formula': mol.molecular_formula(),
                    'MLP_LogP': round(logp, 4)
                })
                print(f"  MLP LogP: {logp:.4f}")
                
                # Verbose output (idea addition)
                if getattr(process_files, 'verbose', False):
                    print("  Atom Breakdown:")
                    for idx, val in enumerate(service.atom_fragment_values):
                        atom = mol.atoms[idx]
                        if atom.symbol != 'H':
                            print(f"    Atom {idx+1:2d} ({atom.symbol:2s} - {atom.sybyl_type:5s}): {val:7.4f}")
            except Exception as e:
                print(f"  Error calculating LogP for {mol_id}: {e}")
                
    if results:
        df = pd.DataFrame(results)
        if output_csv:
            df.to_csv(output_csv, index=False)
            print(f"\nSuccess! Results saved to {output_csv}")
        else:
            print("\nCalculation Results:")
            print(df.to_string(index=False))
    else:
        print("No results generated.")

def main():
    parser = argparse.ArgumentParser(description="Calculate MLP LogP for SDF, MOL, and MOL2 files.")
    parser.add_argument("input", nargs="+", help="Input files or directories")
    parser.add_argument("--output", "-o", help="Output CSV file path")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print detailed atom breakdown")
    
    args = parser.parse_args()
    
    # Store verbose flag on the function object (simple way to pass it)
    process_files.verbose = args.verbose
    
    # Expand directories if any
    expanded_paths = []
    for p in args.input:
        if os.path.isdir(p):
            for f in os.listdir(p):
                if f.lower().endswith(('.sdf', '.mol', '.mol2')):
                    expanded_paths.append(os.path.join(p, f))
        else:
            expanded_paths.append(p)
            
    if not expanded_paths:
        print("No valid chemical files found.")
        sys.exit(1)
        
    process_files(expanded_paths, args.output)

if __name__ == "__main__":
    main()
