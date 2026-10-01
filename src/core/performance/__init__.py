"""
PyChem Performance Optimization Module

Provides multiprocessing-based optimizations for file loading and startup
to improve performance on systems with limited resources (8GB RAM, i7 3rd gen).

Usage:
    from src.core.performance import ParallelFileLoader, StartupOptimizer
    
    # Parallel file loading
    loader = ParallelFileLoader()
    molecule = loader.load_file("path/to/file.mol2")
    
    # Startup optimization
    optimizer = StartupOptimizer()
    optimizer.optimize_imports()
"""

from .parallel_loader import ParallelFileLoader
from .startup_optimizer import StartupOptimizer
from .profiling import (
    PerformanceProfiler,
    get_profiler,
    profile_operation,
    timed_function,
    print_performance_report
)

__all__ = [
    'ParallelFileLoader',
    'StartupOptimizer',
    'PerformanceProfiler',
    'get_profiler',
    'profile_operation',
    'timed_function',
    'print_performance_report'
]
