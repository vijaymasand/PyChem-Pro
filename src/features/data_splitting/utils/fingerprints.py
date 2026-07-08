"""
Fingerprint computation using PyChem-Pro's native engine.

Supports:
  - Morgan (ECFP-like, radius=2, 2048 bits)
  - Topological (path-based Daylight-like, 2048 bits)
  - MACCS Keys (166 structural keys)
"""
import multiprocessing
import numpy as np
import pandas as pd

# Direct PyChem-Pro internal imports
from src.features.smiles_parser.services.parser import parse_smiles
from src.features.cheminformatics.topology.fingerprints import (
    get_morgan_fingerprint,
    get_topological_fingerprint,
    get_maccs_keys,
)

# -- Fingerprint type constants (used in GUI combo box) -----------------------
FP_MORGAN      = "Morgan (ECFP4, 2048 bits)"
FP_TOPOLOGICAL = "Topological (Path, 2048 bits)"
FP_MACCS       = "MACCS Keys (166 bits)"
FP_COMBINED    = "Combined (Morgan + Topological)"

ALL_FP_TYPES = [FP_MORGAN, FP_TOPOLOGICAL, FP_MACCS, FP_COMBINED]


def _compute_fp_single(args):
    """Worker: parse SMILES and return fingerprint vector for one molecule."""
    smiles, name, fp_type = args
    try:
        mol = parse_smiles(smiles)
        if mol is None:
            return None, name

        if fp_type == FP_MORGAN:
            vec = get_morgan_fingerprint(mol, radius=2, n_bits=2048)
        elif fp_type == FP_TOPOLOGICAL:
            vec = get_topological_fingerprint(mol, min_path=1, max_path=7, n_bits=2048)
        elif fp_type == FP_MACCS:
            vec = get_maccs_keys(mol)
        elif fp_type == FP_COMBINED:
            morgan = get_morgan_fingerprint(mol, radius=2, n_bits=2048)
            topo   = get_topological_fingerprint(mol, min_path=1, max_path=7, n_bits=2048)
            vec = morgan + topo          # concatenate -> 4096 bits
        else:
            vec = get_morgan_fingerprint(mol, radius=2, n_bits=2048)

        return vec, name
    except Exception:
        return None, name


ACTIVE_POOLS = []

def compute_fingerprints(smiles_list, names_list, fp_type=FP_MORGAN, n_jobs=None):
    """
    Compute fingerprints for a list of SMILES in parallel.

    Returns
    -------
    desc_df      : pd.DataFrame  - one row per valid molecule, columns FP_0 ... FP_N
    valid_indices: list[int]     - original row indices of successfully processed molecules
    """
    if n_jobs is None or n_jobs < 1:
        n_jobs = max(1, multiprocessing.cpu_count() // 2)

    args = [(smi, name, fp_type) for smi, name in zip(smiles_list, names_list)]

    with multiprocessing.Pool(processes=n_jobs) as pool:
        ACTIVE_POOLS.append(pool)
        try:
            results = pool.map(_compute_fp_single, args)
        finally:
            if pool in ACTIVE_POOLS:
                ACTIVE_POOLS.remove(pool)

    rows = []
    valid_indices = []
    for idx, (vec, name) in enumerate(results):
        if vec is not None:
            rows.append(vec)
            valid_indices.append(idx)

    if not rows:
        return pd.DataFrame(), []

    n_bits = len(rows[0])
    col_names = [f"FP_{i}" for i in range(n_bits)]
    desc_df = pd.DataFrame(rows, columns=col_names)
    return desc_df, valid_indices
