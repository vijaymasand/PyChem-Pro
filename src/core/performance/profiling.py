"""
Performance Profiling Module — Monitor and analyze PyChem performance.

Provides timing utilities, memory profiling, and bottleneck detection for
identifying slow operations in startup, file loading, and rendering.
"""

import time
import functools
import logging
from typing import Dict, List, Optional, Callable, Any, Tuple
from collections import defaultdict
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class PerformanceProfiler:
    """
    Central performance profiler for PyChem.
    
    Tracks timing for key operations:
    - Startup phases
    - File loading
    - Rendering
    - Cheminformatics calculations
    
    Usage:
        profiler = PerformanceProfiler()
        
        with profiler.time_operation("file_load"):
            molecule = load_mol2("large.mol2")
        
        print(profiler.get_report())
    """
    
    def __init__(self):
        self.timings: Dict[str, List[float]] = defaultdict(list)
        self.counts: Dict[str, int] = defaultdict(int)
        self.current_timers: Dict[str, float] = {}
        self.enabled = True
        
    def enable(self):
        """Enable profiling."""
        self.enabled = True
        
    def disable(self):
        """Disable profiling."""
        self.enabled = False
    
    @contextmanager
    def time_operation(self, name: str):
        """Context manager to time an operation."""
        if not self.enabled:
            yield
            return
            
        start = time.perf_counter()
        try:
            yield
        finally:
            elapsed = time.perf_counter() - start
            self.timings[name].append(elapsed)
            self.counts[name] += 1
            logger.debug(f"Operation '{name}' took {elapsed:.4f}s")
    
    def timed(self, name: Optional[str] = None):
        """Decorator to time a function."""
        def decorator(func: Callable) -> Callable:
            op_name = name or func.__name__
            
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                if not self.enabled:
                    return func(*args, **kwargs)
                    
                with self.time_operation(op_name):
                    return func(*args, **kwargs)
            
            return wrapper
        return decorator
    
    def record_timing(self, name: str, elapsed: float):
        """Manually record a timing."""
        if self.enabled:
            self.timings[name].append(elapsed)
            self.counts[name] += 1
    
    def get_stats(self, name: str) -> Dict[str, float]:
        """Get statistics for a timed operation."""
        times = self.timings.get(name, [])
        if not times:
            return {'count': 0, 'total': 0, 'mean': 0, 'min': 0, 'max': 0}
        
        return {
            'count': len(times),
            'total': sum(times),
            'mean': sum(times) / len(times),
            'min': min(times),
            'max': max(times)
        }
    
    def get_report(self) -> str:
        """Generate formatted performance report."""
        lines = [
            "PyChem Performance Report",
            "=" * 70,
            f"{'Operation':<40} {'Count':>8} {'Total':>10} {'Mean':>10}",
            "-" * 70
        ]
        
        # Sort by total time (most impactful first)
        sorted_ops = sorted(
            self.timings.keys(),
            key=lambda k: sum(self.timings[k]),
            reverse=True
        )
        
        for op in sorted_ops:
            stats = self.get_stats(op)
            lines.append(
                f"{op:<40} {stats['count']:>8} "
                f"{stats['total']:>9.3f}s {stats['mean']:>9.3f}s"
            )
        
        return "\n".join(lines)
    
    def get_bottlenecks(self, threshold: float = 0.5) -> List[Tuple[str, float]]:
        """Identify operations exceeding threshold time."""
        bottlenecks = []
        for name, times in self.timings.items():
            total = sum(times)
            if total >= threshold:
                bottlenecks.append((name, total))
        return sorted(bottlenecks, key=lambda x: x[1], reverse=True)
    
    def reset(self):
        """Clear all timing data."""
        self.timings.clear()
        self.counts.clear()
        self.current_timers.clear()


class MemoryProfiler:
    """
    Memory usage profiler for tracking memory consumption.
    
    Useful for identifying memory leaks and high-memory operations
    on systems with limited RAM (8GB).
    """
    
    def __init__(self):
        self.snapshots: List[Dict[str, Any]] = []
        self._has_psutil = False
        
        try:
            import psutil
            self._has_psutil = True
            self._process = psutil.Process()
        except ImportError:
            logger.warning("psutil not available, memory profiling limited")
    
    def snapshot(self, label: str = ""):
        """Take a memory snapshot."""
        import sys
        
        snapshot = {
            'label': label,
            'time': time.time(),
            'rss_mb': 0,
            'vms_mb': 0,
            'python_objects': len(gc.get_objects()) if 'gc' in dir() else 0
        }
        
        if self._has_psutil:
            mem_info = self._process.memory_info()
            snapshot['rss_mb'] = mem_info.rss / 1024 / 1024
            snapshot['vms_mb'] = mem_info.vms / 1024 / 1024
        
        self.snapshots.append(snapshot)
        logger.info(f"Memory snapshot '{label}': RSS={snapshot['rss_mb']:.1f}MB")
        
        return snapshot
    
    def get_usage_report(self) -> str:
        """Generate memory usage report."""
        if not self.snapshots:
            return "No memory snapshots available"
        
        lines = ["Memory Usage Report", "=" * 50]
        
        prev_rss = None
        for snap in self.snapshots:
            delta = ""
            if prev_rss is not None:
                change = snap['rss_mb'] - prev_rss
                delta = f" ({change:+.1f}MB)"
            
            lines.append(f"{snap['label']:<30} {snap['rss_mb']:>8.1f}MB{delta}")
            prev_rss = snap['rss_mb']
        
        return "\n".join(lines)
    
    def check_memory_pressure(self, threshold_mb: float = 6000) -> bool:
        """Check if memory usage is approaching threshold (for 8GB systems)."""
        if not self._has_psutil:
            return False
            
        mem_info = self._process.memory_info()
        rss_mb = mem_info.rss / 1024 / 1024
        
        if rss_mb > threshold_mb:
            logger.warning(f"High memory usage: {rss_mb:.1f}MB (threshold: {threshold_mb}MB)")
            return True
        return False


class FileLoadProfiler:
    """
    Specialized profiler for file loading operations.
    
    Tracks per-format statistics and identifies slow-loading files.
    """
    
    def __init__(self):
        self.file_timings: Dict[str, Dict[str, Any]] = {}
        self.format_stats: Dict[str, Dict[str, Any]] = defaultdict(lambda: {
            'count': 0,
            'total_time': 0,
            'total_size_kb': 0,
            'atoms': 0
        })
    
    def record_load(self, filepath: str, file_size_kb: float, 
                   num_atoms: int, elapsed: float, file_format: str):
        """Record a file loading operation."""
        self.file_timings[filepath] = {
            'size_kb': file_size_kb,
            'atoms': num_atoms,
            'time': elapsed,
            'format': file_format,
            'kb_per_second': file_size_kb / elapsed if elapsed > 0 else 0
        }
        
        stats = self.format_stats[file_format]
        stats['count'] += 1
        stats['total_time'] += elapsed
        stats['total_size_kb'] += file_size_kb
        stats['atoms'] += num_atoms
    
    def get_format_report(self) -> str:
        """Generate per-format loading statistics."""
        lines = ["File Loading Statistics by Format", "=" * 70]
        
        for fmt, stats in sorted(self.format_stats.items()):
            avg_time = stats['total_time'] / stats['count'] if stats['count'] > 0 else 0
            avg_size = stats['total_size_kb'] / stats['count'] if stats['count'] > 0 else 0
            avg_atoms = stats['atoms'] / stats['count'] if stats['count'] > 0 else 0
            
            lines.extend([
                f"\nFormat: {fmt.upper()}",
                f"  Files loaded: {stats['count']}",
                f"  Avg time: {avg_time:.3f}s",
                f"  Avg size: {avg_size:.1f}KB",
                f"  Avg atoms: {avg_atoms:.0f}",
                f"  Total time: {stats['total_time']:.3f}s"
            ])
        
        return "\n".join(lines)
    
    def get_slow_files(self, threshold: float = 1.0) -> List[Tuple[str, float]]:
        """Identify files taking longer than threshold to load."""
        slow = [
            (path, data['time'])
            for path, data in self.file_timings.items()
            if data['time'] > threshold
        ]
        return sorted(slow, key=lambda x: x[1], reverse=True)


# Global profiler instance for application-wide use
_global_profiler: Optional[PerformanceProfiler] = None
_global_memory_profiler: Optional[MemoryProfiler] = None


def get_profiler() -> PerformanceProfiler:
    """Get or create global performance profiler."""
    global _global_profiler
    if _global_profiler is None:
        _global_profiler = PerformanceProfiler()
    return _global_profiler


def get_memory_profiler() -> MemoryProfiler:
    """Get or create global memory profiler."""
    global _global_memory_profiler
    if _global_memory_profiler is None:
        _global_memory_profiler = MemoryProfiler()
    return _global_memory_profiler


@contextmanager
def profile_operation(name: str):
    """Convenience context manager using global profiler."""
    profiler = get_profiler()
    with profiler.time_operation(name):
        yield


def timed_function(name: Optional[str] = None):
    """Convenience decorator using global profiler."""
    profiler = get_profiler()
    return profiler.timed(name)


def print_performance_report():
    """Print the global performance report."""
    profiler = get_profiler()
    print(profiler.get_report())
