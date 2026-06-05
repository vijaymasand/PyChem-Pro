"""
PyDes Public API for Batch Processing

Allows easy integration with Jupyter Notebooks and external Python scripts.
"""
from typing import List, Dict, Any, Union
import os
import pandas as pd
from src.core.domain.models.molecule import Molecule
from .engine import PyDesEngine

def pydes(files_or_molecules: List[Union[str, Molecule]], n_jobs: int = -1, export_csv: str = None) -> pd.DataFrame:
    """
    Calculate PyDes descriptors for a list of file paths or Molecule objects.
    
    Args:
        files_or_molecules: List of file paths (MOL, SDF, MOL2) or PyChem Molecule objects.
        n_jobs: Number of CPU cores to use (-1 for all cores).
        export_csv: Optional file path to save the results directly to CSV.
        
    Returns:
        pandas.DataFrame containing all calculated descriptors with Molecule_Name as the first column.
    """
    molecules_to_process = []
    
    # Check if we need to load files
    if files_or_molecules and isinstance(files_or_molecules[0], str):
        from pychem.api import load
        for file_path in files_or_molecules:
            try:
                # If it's an SDF with multiple molecules, load might return a list or we might need a specific SDF loader
                # Using PyChem's standard loader
                mol = load(file_path, parallel=False)
                if mol:
                    # Set a name based on file if not present
                    if not getattr(mol, 'name', None):
                        mol.name = os.path.basename(file_path)
                    molecules_to_process.append(mol)
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
    else:
        molecules_to_process = files_or_molecules
        
    # Run batch calculation
    results = PyDesEngine.calculate_batch(molecules_to_process, n_jobs=n_jobs)
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Ensure Molecule_Name is the first column
    if 'Molecule_Name' in df.columns:
        cols = ['Molecule_Name'] + [col for col in df.columns if col != 'Molecule_Name']
        df = df[cols]
        
    # Export if requested
    if export_csv:
        df.to_csv(export_csv, index=False)
        print(f"Results exported to {export_csv}")
        
    return df
