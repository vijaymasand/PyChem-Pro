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
    Optimized with NumPy vectorization for better performance.
    
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

    # Convert to sparse matrix for efficient edge operations
    from scipy.sparse import csr_matrix
    graph_sparse = csr_matrix(graph)
    
    # Get edges and weights in COO format for efficient iteration
    rows, cols = graph_sparse.nonzero()
    weights = graph_sparse.data
    
    # Stack edges
    edges = np.column_stack((rows, cols))
    n_edges = len(edges)
    
    if n_edges == 0:
        return Y

    # Pre-compute constants to avoid repeated calculations
    a_b = a * b
    two_b = 2.0 * b
    
    for epoch in range(n_epochs):
        alpha = learning_rate * (1.0 - (epoch / float(n_epochs)))

        # Shuffle edges per epoch
        order = rng.permutation(n_edges)
        
        # Process edges in batches for better cache locality
        batch_size = min(1000, n_edges)
        for batch_start in range(0, n_edges, batch_size):
            batch_end = min(batch_start + batch_size, n_edges)
            batch_indices = order[batch_start:batch_end]
            
            for idx in batch_indices:
                i, j = edges[idx]
                w = weights[idx]

                # Vectorized distance computation
                diff = Y[i] - Y[j]
                dist2 = np.dot(diff, diff) + 1e-12

                # Attractive force gradient (vectorized)
                dist2_pow_b = dist2 ** b
                attr_coeff = (-a_b * (dist2_pow_b / dist2)) / (1.0 + a * dist2_pow_b)
                attr_grad = np.clip(w * attr_coeff * diff, -4.0, 4.0)

                Y[i] += alpha * attr_grad
                Y[j] -= alpha * attr_grad

                # Negative sampling (Repulsive force gradient)
                # Sample multiple negatives at once for efficiency
                neg_samples = rng.integers(0, n_samples, size=negative_sample_rate)
                valid_neg = neg_samples[(neg_samples != i) & (neg_samples != j)]
                
                for k in valid_neg:
                    diff_n = Y[i] - Y[k]
                    dist2_n = np.dot(diff_n, diff_n) + 1e-12
                    
                    dist2_n_pow_b = dist2_n ** b
                    rep_coeff = two_b / ((1e-3 + dist2_n) * (1.0 + a * dist2_n_pow_b))
                    rep_grad = np.clip((1.0 - w) * rep_coeff * diff_n, -4.0, 4.0)
                    
                    Y[i] += alpha * rep_grad

    return Y
