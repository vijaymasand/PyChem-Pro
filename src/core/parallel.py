"""
Centralized multiprocessing executor for PyChem.
Uses 50% of available CPU cores by default.

Cross-platform: forces 'spawn' start method on all platforms to avoid:
- macOS: fork + Qt = crash (CoreFoundation assertion)
- Windows: fork not available (spawn is the only option)
- Linux: fork works but spawn is safer with Qt
"""
import os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

# Force spawn on all platforms for Qt safety
_mp_context = mp.get_context('spawn')


class ParallelExecutor:
    def __init__(self, core_fraction: float = 0.5, max_workers: int = None):
        self._core_fraction = core_fraction
        self._max_workers = max_workers
        self._pool = None

    @property
    def num_workers(self) -> int:
        if self._max_workers:
            return self._max_workers
        cores = os.cpu_count() or 4
        return max(1, int(cores * self._core_fraction))

    @property
    def pool(self) -> ProcessPoolExecutor:
        if self._pool is None:
            cores = os.cpu_count() or 4
            workers = self.num_workers
            print(f"[ParallelExecutor] Initializing ProcessPoolExecutor: "
                  f"Total logical CPU cores detected = {cores}, "
                  f"Core allocation fraction = {self._core_fraction * 100}%, "
                  f"Workers spawned = {workers}")
            self._pool = ProcessPoolExecutor(
                max_workers=workers,
                mp_context=_mp_context
            )
        return self._pool

    def map(self, fn, chunks, timeout=60):
        if not chunks:
            return []
        if len(chunks) <= 1 or self.num_workers <= 1:
            return [fn(chunk) for chunk in chunks]
        try:
            futures = [self.pool.submit(fn, chunk) for chunk in chunks]
            return [f.result(timeout=timeout) for f in futures]
        except Exception:
            return [fn(chunk) for chunk in chunks]

    def map_unordered(self, fn, chunks, timeout=60):
        if not chunks:
            return
        if len(chunks) <= 1 or self.num_workers <= 1:
            for chunk in chunks:
                yield fn(chunk)
            return
        try:
            futures = {self.pool.submit(fn, c): i for i, c in enumerate(chunks)}
            for future in as_completed(futures, timeout=timeout):
                yield future.result()
        except Exception:
            for chunk in chunks:
                yield fn(chunk)

    def shutdown(self):
        if self._pool:
            self._pool.shutdown(wait=True)
            self._pool = None
