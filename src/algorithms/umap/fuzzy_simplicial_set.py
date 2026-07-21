"""
Fuzzy Simplicial Set Construction and Adaptive Local Scaling (sigma, rho) with Multiprocessing.
Follows UMAP mathematical formulations for Riemannian metric local adaptation.
"""
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from .metrics import resolve_n_jobs

def _smooth_knn_worker(args):
    knn_dists_chunk, chunk_idx, k, n_iter, target, bandwidth = args
    n_chunk = len(knn_dists_chunk)
    rho_chunk = np.zeros(n_chunk, dtype=np.float64)
    sigma_chunk = np.zeros(n_chunk, dtype=np.float64)

    for i in range(n_chunk):
        dists = knn_dists_chunk[i]
        non_zero_dists = dists[dists > 0]
        if len(non_zero_dists) > 0:
            rho_chunk[i] = non_zero_dists[0]
        else:
            rho_chunk[i] = 0.0

        lo = 1e-8
        hi = 1e8
        mid = 1.0

        for _ in range(n_iter):
            mid = 0.5 * (lo + hi)
            val = np.sum(np.exp(-np.maximum(0.0, dists - rho_chunk[i]) / mid))
            if val > target:
                hi = mid
            else:
                lo = mid

        sigma_chunk[i] = mid

    return chunk_idx, rho_chunk, sigma_chunk

def smooth_knn_dist(
    knn_dists: np.ndarray,
    k: int,
    n_iter: int = 64,
    bandwidth: float = 1.0,
    n_jobs: int = None
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes adaptive local scale parameters (rho, sigma) for each sample in parallel.
    
    Parameters
    ----------
    knn_dists : np.ndarray of shape (n_samples, n_neighbors)
    k : int
        Number of neighbors target.
    n_iter : int
        Number of binary search iterations.
    bandwidth : float
        Bandwidth factor.
    n_jobs : int, optional
        Number of parallel jobs (defaults to 50% CPU cores).
        
    Returns
    -------
    rho : np.ndarray of shape (n_samples,)
    sigma : np.ndarray of shape (n_samples,)
    """
    n_samples = knn_dists.shape[0]
    target = np.log2(k) * bandwidth
    n_jobs = resolve_n_jobs(n_jobs)

    if n_samples <= 200 or n_jobs == 1:
        rho = np.zeros(n_samples, dtype=np.float64)
        sigma = np.zeros(n_samples, dtype=np.float64)
        for i in range(n_samples):
            dists = knn_dists[i]
            non_zero_dists = dists[dists > 0]
            rho[i] = non_zero_dists[0] if len(non_zero_dists) > 0 else 0.0
            lo, hi = 1e-8, 1e8
            for _ in range(n_iter):
                mid = 0.5 * (lo + hi)
                val = np.sum(np.exp(-np.maximum(0.0, dists - rho[i]) / mid))
                if val > target:
                    hi = mid
                else:
                    lo = mid
            sigma[i] = mid
        return rho, sigma

    chunks = np.array_split(np.arange(n_samples), n_jobs)
    tasks = [(knn_dists[chunk], chunk, k, n_iter, target, bandwidth) for chunk in chunks]

    rho = np.zeros(n_samples, dtype=np.float64)
    sigma = np.zeros(n_samples, dtype=np.float64)

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(_smooth_knn_worker, tasks)
        for chunk_idx, sub_rho, sub_sigma in results:
            rho[chunk_idx] = sub_rho
            sigma[chunk_idx] = sub_sigma

    return rho, sigma

def fuzzy_simplicial_set(
    D: np.ndarray,
    knn_indices: np.ndarray,
    knn_dists: np.ndarray,
    n_neighbors: int,
    set_op_mix_ratio: float = 1.0,
    n_jobs: int = None
) -> np.ndarray:
    """
    Builds the high-dimensional fuzzy graph / fuzzy simplicial set matrix P in parallel.
    
    Parameters
    ----------
    D : np.ndarray of shape (n_samples, n_samples)
    knn_indices : np.ndarray of shape (n_samples, k)
    knn_dists : np.ndarray of shape (n_samples, k)
    n_neighbors : int
    set_op_mix_ratio : float
    n_jobs : int, optional
        
    Returns
    -------
    P : np.ndarray of shape (n_samples, n_samples)
        Symmetrized fuzzy membership graph.
    """
    n_samples = D.shape[0]
    n_jobs = resolve_n_jobs(n_jobs)
    rho, sigma = smooth_knn_dist(knn_dists, k=n_neighbors, n_jobs=n_jobs)

    A = np.zeros((n_samples, n_samples), dtype=np.float64)

    def _graph_worker(chunk_idx):
        for i in chunk_idx:
            indices = knn_indices[i]
            dists = D[i, indices]
            val = np.exp(-np.maximum(0.0, dists - rho[i]) / np.maximum(sigma[i], 1e-12))
            A[i, indices] = val

    if n_samples <= 200 or n_jobs == 1:
        _graph_worker(np.arange(n_samples))
    else:
        chunks = np.array_split(np.arange(n_samples), n_jobs)
        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            list(executor.map(_graph_worker, chunks))

    # Symmetrization via Fuzzy Set Union: A + A^T - A * A^T
    A_transpose = A.T
    fuzzy_union = A + A_transpose - A * A_transpose
    fuzzy_intersection = A * A_transpose

    P = set_op_mix_ratio * fuzzy_union + (1.0 - set_op_mix_ratio) * fuzzy_intersection
    np.fill_diagonal(P, 0.0)

    p_max = np.max(P)
    if p_max > 0:
        P /= p_max

    return P
