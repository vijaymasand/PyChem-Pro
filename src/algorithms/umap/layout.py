"""
UMAP Layout Optimization and Initialization (Spectral & SGD).
Provides graph spectral initialization and SGD manifold embedding optimizer.
"""
import numpy as np
from scipy.sparse.csgraph import connected_components

def spectral_layout(graph: np.ndarray, n_components: int = 2, random_state: int = 42) -> np.ndarray:
    """
    Computes spectral initialization from the graph Laplacian of fuzzy graph P.
    
    Parameters
    ----------
    graph : np.ndarray of shape (n_samples, n_samples)
    n_components : int
    random_state : int
    
    Returns
    -------
    embedding : np.ndarray of shape (n_samples, n_components)
    """
    n_samples = graph.shape[0]
    rng = np.random.default_rng(random_state)

    try:
        # Check connected components
        n_components_cc, _ = connected_components(graph, directed=False)
        if n_components_cc > 1:
            # Fallback to PCA / random if graph is disconnected
            return 0.0001 * rng.standard_normal((n_samples, n_components))

        # Normalized Laplacian: L = I - D^{-1/2} P D^{-1/2}
        degrees = np.sum(graph, axis=1)
        degrees_inv_sqrt = np.power(np.maximum(degrees, 1e-12), -0.5)
        D_inv_sqrt = np.diag(degrees_inv_sqrt)
        
        normalized_graph = D_inv_sqrt @ graph @ D_inv_sqrt
        
        # Symmetric eigendecomposition
        eigvals, eigvecs = np.linalg.eigh(normalized_graph)
        
        # Sort eigenvectors by descending eigenvalues (smallest eigenvalues of L)
        idx = np.argsort(eigvals)[::-1]
        eigvecs = eigvecs[:, idx]
        
        # Select components 1..n_components (skipping trivial first eigenvector)
        if eigvecs.shape[1] > n_components:
            embedding = eigvecs[:, 1:n_components + 1]
        else:
            embedding = eigvecs[:, :n_components]
            
        # Scale embedding
        std = np.std(embedding, axis=0, keepdims=True)
        std[std == 0] = 1.0
        embedding = (embedding / std) * 0.0001

        return embedding
    except Exception:
        return 0.0001 * rng.standard_normal((n_samples, n_components))

def optimize_layout(
    graph: np.ndarray,
    n_components: int = 2,
    n_epochs: int = 200,
    learning_rate: float = 1.0,
    a: float = 1.5769,
    b: float = 0.8951,
    negative_sample_rate: int = 5,
    init_embedding: np.ndarray = None,
    random_state: int = 42
) -> np.ndarray:
    """
    Optimizes low-dimensional coordinates Y using Stochastic Gradient Descent.
    
    Parameters
    ----------
    graph : np.ndarray of shape (n_samples, n_samples)
    n_components : int
    n_epochs : int
    learning_rate : float
    a : float
    b : float
    negative_sample_rate : int
    init_embedding : np.ndarray, optional
    random_state : int
    
    Returns
    -------
    Y : np.ndarray of shape (n_samples, n_components)
    """
    rng = np.random.default_rng(random_state)
    n_samples = graph.shape[0]

    if init_embedding is not None:
        Y = np.copy(init_embedding)
    else:
        Y = 0.0001 * rng.standard_normal((n_samples, n_components))

    # Extract non-zero edge indices and weights
    edge_mask = graph > 0
    edges = np.argwhere(edge_mask)
    weights = graph[edge_mask]

    n_edges = len(edges)
    if n_edges == 0:
        return Y

    # Pre-calculate edge sampling probabilities / repeat indices
    for epoch in range(n_epochs):
        alpha = learning_rate * (1.0 - (epoch / float(n_epochs)))

        # Shuffle edges per epoch
        order = rng.permutation(n_edges)
        
        for idx in order:
            i, j = edges[idx]
            w = weights[idx]

            diff = Y[i] - Y[j]
            dist2 = np.dot(diff, diff) + 1e-12

            # 1. Attractive force gradient
            attr_coeff = (-2.0 * a * b * (dist2 ** (b - 1.0))) / (1.0 + a * (dist2 ** b))
            attr_grad = np.clip(w * attr_coeff * diff, -4.0, 4.0)

            Y[i] += alpha * attr_grad
            Y[j] -= alpha * attr_grad

            # 2. Negative sampling (Repulsive force gradient)
            for _ in range(negative_sample_rate):
                k = rng.integers(0, n_samples)
                if k == i or k == j:
                    continue

                diff_n = Y[i] - Y[k]
                dist2_n = np.dot(diff_n, diff_n) + 1e-12

                rep_coeff = (2.0 * b) / ((1e-3 + dist2_n) * (1.0 + a * (dist2_n ** b)))
                rep_grad = np.clip((1.0 - w) * rep_coeff * diff_n, -4.0, 4.0)

                Y[i] += alpha * rep_grad

    return Y
