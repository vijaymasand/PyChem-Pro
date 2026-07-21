"""
Distance Metrics & Pairwise Distance Computer for UMAP with Multiprocessing.
Supports standard continuous metrics as well as cheminformatics bit-vector metrics.
"""
import os
import numpy as np
from scipy.spatial.distance import pdist, squareform, cdist
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

def resolve_n_jobs(n_jobs: int = None) -> int:
    """Resolves n_jobs. If None or <= 0, defaults to 50% of available CPU cores."""
    cpu_cnt = os.cpu_count() or 2
    if n_jobs is None or n_jobs <= 0:
        return max(1, cpu_cnt // 2)
    return min(n_jobs, cpu_cnt)

def pairwise_distances(X: np.ndarray, Y: np.ndarray = None, metric: str = 'euclidean', n_jobs: int = None) -> np.ndarray:
    """
    Computes pairwise distance matrix between X (and optionally Y) in parallel.
    
    Parameters
    ----------
    X : np.ndarray of shape (n_samples_X, n_features)
    Y : np.ndarray of shape (n_samples_Y, n_features), optional
    metric : str
        Distance metric name ('euclidean', 'manhattan', 'cosine', 'tanimoto', 'jaccard', 'hamming', 'dice')
    n_jobs : int, optional
        Number of parallel jobs. Defaults to 50% CPU cores.
        
    Returns
    -------
    dist_matrix : np.ndarray of shape (n_samples_X, n_samples_Y or n_samples_X)
    """
    X = np.asarray(X, dtype=np.float64)
    metric = metric.lower()
    n_jobs = resolve_n_jobs(n_jobs)

    if metric in ['jaccard', 'tanimoto']:
        return _tanimoto_distances(X, Y)
    elif metric == 'dice':
        return _dice_distances(X, Y)
    elif metric in ['euclidean', 'l2']:
        metric_name = 'euclidean'
    elif metric in ['manhattan', 'l1', 'cityblock']:
        metric_name = 'cityblock'
    elif metric == 'cosine':
        metric_name = 'cosine'
    elif metric == 'hamming':
        metric_name = 'hamming'
    else:
        metric_name = metric

    if Y is None:
        # If dataset is small, compute directly
        if len(X) <= 200 or n_jobs == 1:
            try:
                return squareform(pdist(X, metric=metric_name))
            except Exception:
                return cdist(X, X, metric=metric_name)
        
        # Parallel chunked cdist computation
        n_samples = len(X)
        dist_matrix = np.zeros((n_samples, n_samples), dtype=np.float64)
        chunks = np.array_split(np.arange(n_samples), n_jobs)

        def worker_dist(chunk_idx):
            return chunk_idx, cdist(X[chunk_idx], X, metric=metric_name)

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            results = executor.map(worker_dist, chunks)
            for chunk_idx, sub_dists in results:
                dist_matrix[chunk_idx] = sub_dists

        return dist_matrix
    else:
        Y = np.asarray(Y, dtype=np.float64)
        if len(X) <= 200 or n_jobs == 1:
            return cdist(X, Y, metric=metric_name)

        n_samples = len(X)
        dist_matrix = np.zeros((n_samples, len(Y)), dtype=np.float64)
        chunks = np.array_split(np.arange(n_samples), n_jobs)

        def worker_dist_y(chunk_idx):
            return chunk_idx, cdist(X[chunk_idx], Y, metric=metric_name)

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            results = executor.map(worker_dist_y, chunks)
            for chunk_idx, sub_dists in results:
                dist_matrix[chunk_idx] = sub_dists

        return dist_matrix

def _tanimoto_distances(X: np.ndarray, Y: np.ndarray = None) -> np.ndarray:
    """Computes continuous / binary Tanimoto (Jaccard) distance matrix using BLAS matrix operations."""
    if Y is None:
        Y = X

    dot_product = np.dot(X, Y.T)
    norm_x = np.sum(X ** 2, axis=1, keepdims=True)
    norm_y = np.sum(Y ** 2, axis=1, keepdims=True).T

    denominator = norm_x + norm_y - dot_product
    denominator = np.maximum(denominator, 1e-12)
    similarity = np.clip(dot_product / denominator, 0.0, 1.0)
    return 1.0 - similarity

def _dice_distances(X: np.ndarray, Y: np.ndarray = None) -> np.ndarray:
    """Computes Dice distance matrix."""
    if Y is None:
        Y = X

    dot_product = np.dot(X, Y.T)
    norm_x = np.sum(X ** 2, axis=1, keepdims=True)
    norm_y = np.sum(Y ** 2, axis=1, keepdims=True).T

    denominator = norm_x + norm_y
    denominator = np.maximum(denominator, 1e-12)
    similarity = np.clip(2.0 * dot_product / denominator, 0.0, 1.0)
    return 1.0 - similarity

def _knn_chunk_worker(args):
    D_chunk, chunk_idx, k = args
    n_chunk = len(D_chunk)
    indices = np.zeros((n_chunk, k), dtype=int)
    dists = np.zeros((n_chunk, k), dtype=np.float64)

    for i in range(n_chunk):
        global_idx = chunk_idx[i]
        row = D_chunk[i]
        sorted_idx = np.argsort(row)
        if sorted_idx[0] == global_idx:
            neighbors = sorted_idx[1:k + 1]
        else:
            neighbors = sorted_idx[sorted_idx != global_idx][:k]

        indices[i] = neighbors
        dists[i] = row[neighbors]

    return chunk_idx, indices, dists

def find_knn(D: np.ndarray, n_neighbors: int, n_jobs: int = None):
    """
    Finds k-nearest neighbors indices and distances given pairwise distance matrix D in parallel.
    
    Parameters
    ----------
    D : np.ndarray of shape (n_samples, n_samples)
    n_neighbors : int
    n_jobs : int, optional
    
    Returns
    -------
    knn_indices : np.ndarray of shape (n_samples, k)
    knn_dists : np.ndarray of shape (n_samples, k)
    """
    n_samples = D.shape[0]
    k = min(n_neighbors, n_samples - 1) if n_samples > 1 else 1
    n_jobs = resolve_n_jobs(n_jobs)

    if n_samples <= 200 or n_jobs == 1:
        knn_indices = np.zeros((n_samples, k), dtype=int)
        knn_dists = np.zeros((n_samples, k), dtype=np.float64)
        for i in range(n_samples):
            row = D[i]
            sorted_idx = np.argsort(row)
            if sorted_idx[0] == i:
                neighbors = sorted_idx[1:k + 1]
            else:
                neighbors = sorted_idx[sorted_idx != i][:k]
            knn_indices[i] = neighbors
            knn_dists[i] = row[neighbors]
        return knn_indices, knn_dists

    chunks = np.array_split(np.arange(n_samples), n_jobs)
    tasks = [(D[chunk], chunk, k) for chunk in chunks]

    knn_indices = np.zeros((n_samples, k), dtype=int)
    knn_dists = np.zeros((n_samples, k), dtype=np.float64)

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(_knn_chunk_worker, tasks)
        for chunk_idx, sub_indices, sub_dists in results:
            knn_indices[chunk_idx] = sub_indices
            knn_dists[chunk_idx] = sub_dists

    return knn_indices, knn_dists
