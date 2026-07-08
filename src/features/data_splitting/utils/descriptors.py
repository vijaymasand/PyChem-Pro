import multiprocessing
from enum import Enum
import numpy as np
import pandas as pd

# Direct PyChem-Pro internal imports
from src.features.smiles_parser.services.parser import parse_smiles
from src.features.layout_3d import generate_3d_coordinates
from src.features.descriptor_calculator.pydes.api import pydes


class DescriptorMode(Enum):
    """How feature vectors are obtained before distance-based splitting."""
    USE_CSV        = "use_csv"        # Columns already present in the CSV
    FINGERPRINTS   = "fingerprints"   # Compute via PyChem-Pro fingerprint engine
    CALCULATE      = "calculate"      # Compute via PyDes (slow, full descriptor set)


def get_feature_matrix(df, smiles_list, names_list,
                       mode: DescriptorMode,
                       csv_descriptor_cols=None,
                       fp_type=None,
                       n_jobs=1):
    """
    Unified entry point that returns (desc_df, valid_indices) regardless of mode.

    Parameters
    ----------
    df                  : original full DataFrame
    smiles_list         : list of SMILES strings (already filtered for non-null)
    names_list          : list of molecule names / IDs
    mode                : DescriptorMode enum value
    csv_descriptor_cols : list[str] - column names to use when mode == USE_CSV
    fp_type             : str - fingerprint type constant from utils.fingerprints
    n_jobs              : int - parallel workers

    Returns
    -------
    desc_df      : pd.DataFrame  (one row per valid molecule, numeric columns)
    valid_indices: list[int]     (row indices into df that were successful)
    """
    if mode == DescriptorMode.USE_CSV:
        if not csv_descriptor_cols:
            raise ValueError("No descriptor columns selected from the CSV.")
        desc_df = df[csv_descriptor_cols].copy().reset_index(drop=True)
        valid_indices = list(range(len(df)))
        return desc_df, valid_indices

    elif mode == DescriptorMode.FINGERPRINTS:
        from .fingerprints import compute_fingerprints, FP_MORGAN
        fp = fp_type if fp_type else FP_MORGAN
        desc_df, valid_indices = compute_fingerprints(smiles_list, names_list,
                                                      fp_type=fp, n_jobs=n_jobs)
        return desc_df, valid_indices

    else:  # DescriptorMode.CALCULATE
        desc_df, valid_indices = compute_dataset_descriptors(smiles_list, names_list,
                                                             n_jobs=n_jobs)
        return desc_df, valid_indices


def _process_single_molecule(smiles_and_name):
    """Worker function to parse and generate 3D coordinates for a single SMILES."""
    smiles, name = smiles_and_name
    try:
        mol = parse_smiles(smiles)
        if mol:
            mol.name = str(name)
            # Generate 3D coordinates for 3D descriptors
            try:
                generate_3d_coordinates(mol, optimize=True)
            except Exception:
                pass
            return mol
    except Exception:
        pass
    return None


ACTIVE_POOLS = []

def compute_dataset_descriptors(smiles_list, names_list, n_jobs=None):
    """
    Computes molecular descriptors for a list of SMILES in parallel.
    Returns:
        desc_df: DataFrame of calculated descriptors.
        valid_indices: List of original indices that were successfully calculated.
    """
    if n_jobs is None or n_jobs < 1:
        n_jobs = max(1, multiprocessing.cpu_count() // 2)

    # 1. Parse and generate 3D coordinates in parallel
    pool_inputs = list(zip(smiles_list, names_list))

    # Use multiprocessing to parse and generate coordinates
    # On Windows, we should guard this when running in GUI, which we will handle.
    with multiprocessing.Pool(processes=n_jobs) as pool:
        ACTIVE_POOLS.append(pool)
        try:
            molecules = pool.map(_process_single_molecule, pool_inputs)
        finally:
            if pool in ACTIVE_POOLS:
                ACTIVE_POOLS.remove(pool)

    # Filter out failed ones
    valid_molecules = []
    valid_indices = []
    for idx, mol in enumerate(molecules):
        if mol is not None:
            valid_molecules.append(mol)
            valid_indices.append(idx)

    if not valid_molecules:
        return pd.DataFrame(), []

    # 2. Run PyDes engine in batch mode
    # pydes already handles multiprocessing internally using n_jobs
    desc_df = pydes(valid_molecules, n_jobs=n_jobs)

    return desc_df, valid_indices


def preprocess_descriptors_numpy(desc_df, standardize=True):
    """
    Cleans and standardizes the descriptor DataFrame using pure NumPy.
    - Selects only numeric columns.
    - Fills NaNs with column means.
    - Removes constant columns.
    - Standardizes to zero mean and unit variance (if standardize=True).
    """
    df = desc_df.copy()
    if 'Molecule_Name' in df.columns:
        df = df.drop(columns=['Molecule_Name'])

    # Filter for numeric columns only
    df = df.select_dtypes(include=[np.number])

    if df.empty:
        raise ValueError("No numeric descriptors found in the descriptor DataFrame.")

    X = df.values.astype(float)

    # Handle NaNs: replace with column means
    with np.errstate(all='ignore'):
        col_means = np.nanmean(X, axis=0)
        col_means = np.nan_to_num(col_means, nan=0.0)

    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_means, inds[1])

    # Remove constant columns (std < 1e-6)
    col_stds = np.std(X, axis=0)
    valid_cols = col_stds > 1e-6

    if not np.any(valid_cols):
        return np.zeros_like(X)

    X = X[:, valid_cols]
    col_stds = col_stds[valid_cols]
    col_means = col_means[valid_cols]

    # Standardize
    if standardize:
        X = (X - col_means) / col_stds
    return X
