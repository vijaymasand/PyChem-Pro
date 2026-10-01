"""
Data Splitting Feature

Provides QSAR dataset partitioning algorithms with a PySide6 GUI
and a headless programmatic API.
"""

from .split_engine import DataSplitEngine, SplitResult, split_dataset
from .algorithms import SPLITTERS, ALGOS_NEED_FEATURES
from .utils.descriptors import DescriptorMode
from .utils.fingerprints import ALL_FP_TYPES

__all__ = [
    'DataSplitEngine',
    'SplitResult',
    'split_dataset',
    'SPLITTERS',
    'ALGOS_NEED_FEATURES',
    'DescriptorMode',
    'ALL_FP_TYPES',
]
