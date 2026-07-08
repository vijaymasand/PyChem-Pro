"""
Data splitting algorithms for QSAR chemical space partitioning.

Available algorithms:
    - Random
    - Scaffold-based
    - CADEX (Kennard-Stone)
    - D-Optimal Design
    - Sphere Exclusion
    - Boruta-based
    - Duplex
"""
from .random import RandomSplitter
from .scaffold import ScaffoldSplitter
from .cadex import CADEXSplitter
from .d_optimal import DOptimalSplitter
from .sphere_exclusion import SphereExclusionSplitter
from .boruta import BorutaSplitter
from .duplex import DuplexSplitter

SPLITTERS = {
    "Random": RandomSplitter,
    "Scaffold-based": ScaffoldSplitter,
    "CADEX (Kennard-Stone)": CADEXSplitter,
    "D-Optimal Design": DOptimalSplitter,
    "Sphere Exclusion": SphereExclusionSplitter,
    "Boruta-based": BorutaSplitter,
    "Duplex": DuplexSplitter,
}

# Algorithms that require a feature matrix (descriptor/fingerprint based)
ALGOS_NEED_FEATURES = [
    "CADEX (Kennard-Stone)", "D-Optimal Design",
    "Sphere Exclusion", "Boruta-based", "Duplex",
]


def get_splitter(name: str, target_ratio: float, n_jobs: int = 1, **kwargs):
    """Instantiates a splitter by name.

    Parameters
    ----------
    name : str
        One of the keys in ``SPLITTERS``.
    target_ratio : float
        Fraction of data for the training set (e.g. 0.8).
    n_jobs : int
        Number of parallel workers for multiprocessing-enabled algorithms.
    **kwargs
        Extra keyword arguments forwarded to the splitter constructor
        (e.g. ``metric='euclidean'``).

    Returns
    -------
    BaseSplitter
        An initialized splitter instance.
    """
    if name not in SPLITTERS:
        raise ValueError(
            f"Unknown splitting algorithm: {name}. "
            f"Available: {list(SPLITTERS.keys())}"
        )
    return SPLITTERS[name](target_ratio=target_ratio, n_jobs=n_jobs, **kwargs)
