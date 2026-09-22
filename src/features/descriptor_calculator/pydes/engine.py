"""
PyDes Engine

Coordinates the calculation of all PyDes descriptors.
Provides methods for multiprocessing calculation of multiple molecules,
including a streaming mode that yields results one molecule at a time
so callers can write to disk incrementally without holding the entire
dataset in memory.
"""
from typing import Dict, Any, List, Callable, Optional
import multiprocessing
import threading
import traceback

from .constitutional import calculate_constitutional
from .topological import calculate_topological


# ---------------------------------------------------------------------------
# Module-level worker — must be at top level so multiprocessing can pickle it
# ---------------------------------------------------------------------------

def _worker_fn(molecule) -> Dict[str, Any]:
    """Top-level worker function for multiprocessing (must be picklable)."""
    return PyDesEngine.calculate_single(molecule)


class PyDesEngine:
    """Engine to calculate all PyDes custom descriptors."""

    @staticmethod
    def calculate_single(molecule) -> Dict[str, Any]:
        """Calculate descriptors for a single molecule."""
        results = {}
        try:
            # Name extraction based on PyChem-Pro Molecule object
            name = getattr(molecule, 'name', '') or getattr(molecule, 'title', 'Unknown')
            results['Molecule_Name'] = name

            # Combine all descriptor types
            const_res = calculate_constitutional(molecule)
            topo_res = calculate_topological(molecule)

            results.update(const_res)
            results.update(topo_res)

            # Derived ATS indices (scaled / squared variants)
            keys = list(topo_res.keys())
            for k in keys:
                if 'ATS' in k:
                    val = topo_res[k]
                    results[f"{k}_norm"] = val / (const_res.get('nHeavyAtoms', 1) or 1)
                    results[f"{k}_sq"] = val ** 2

            # Append the full PyChem-Pro descriptor set
            try:
                from ..descriptor_engine import DescriptorEngine
                engine = DescriptorEngine(enable_cache=False)
                original_results = engine.calculate_all(molecule)
                for k, v in original_results.items():
                    results[k] = v.value
            except Exception as inner_e:
                print(f"[PyDes] Warning: Could not calculate the descriptor engine set: {inner_e}")

        except Exception as e:
            print(f"[PyDes Error] Calculation failed for molecule: {e}")
            traceback.print_exc()

        return results

    # ------------------------------------------------------------------
    # Legacy batch API (kept for backwards-compat / scripting API)
    # ------------------------------------------------------------------

    @classmethod
    def _worker(cls, molecule):
        """Worker function for multiprocessing (legacy)."""
        return cls.calculate_single(molecule)

    @classmethod
    def calculate_batch(cls, molecules: List[Any], n_jobs: int = -1) -> List[Dict[str, Any]]:
        """Calculate descriptors for multiple molecules using multiprocessing.

        NOTE: This loads all results into memory at once. For large datasets,
        prefer ``calculate_batch_streaming`` which writes results incrementally.
        """
        if not molecules:
            return []

        if n_jobs < 1:
            n_jobs = multiprocessing.cpu_count()

        results = []
        if n_jobs == 1 or len(molecules) == 1:
            for m in molecules:
                results.append(cls.calculate_single(m))
            return results

        # Multiprocessing pool
        try:
            with multiprocessing.Pool(processes=n_jobs) as pool:
                results = pool.map(_worker_fn, molecules)
        except Exception as e:
            print(f"Multiprocessing error: {e}. Falling back to sequential.")
            for m in molecules:
                results.append(cls.calculate_single(m))

        return results

    # ------------------------------------------------------------------
    # Streaming batch API (new)
    # ------------------------------------------------------------------

    @classmethod
    def calculate_batch_streaming(
        cls,
        molecules: List[Any],
        result_callback: Callable[[Dict[str, Any], int, int], None],
        n_jobs: int = -1,
        stop_event: Optional[threading.Event] = None,
    ) -> int:
        """Calculate descriptors for multiple molecules, calling *result_callback*
        as soon as each molecule is finished.

        This is the memory-efficient replacement for ``calculate_batch``.
        Results are never accumulated — the caller owns each dict immediately
        after the callback returns.

        Args:
            molecules:        List of PyChem Molecule objects to process.
            result_callback:  Callable(result_dict, molecule_index, total_count).
                              Called from the calling thread (not a worker process).
                              The index is 0-based and counts completed molecules.
            n_jobs:           Number of worker processes (-1 = all CPUs).
            stop_event:       Optional threading.Event. When set, the loop exits
                              after the current in-flight results are consumed.

        Returns:
            Number of molecules actually processed (may be < len(molecules) if
            stop_event was set mid-run).
        """
        if not molecules:
            return 0

        total = len(molecules)

        if n_jobs < 1:
            n_jobs = multiprocessing.cpu_count()

        # Cap workers to the number of molecules (no point spawning more)
        n_jobs = min(n_jobs, total)

        completed = 0

        # Single-process path — avoids Pool overhead for small datasets or n_jobs==1
        if n_jobs == 1:
            for mol in molecules:
                if stop_event and stop_event.is_set():
                    break
                result = cls.calculate_single(mol)
                result_callback(result, completed, total)
                completed += 1
            return completed

        # Multiprocessing path
        try:
            with multiprocessing.Pool(processes=n_jobs) as pool:
                # imap_unordered yields results as soon as any worker finishes —
                # no waiting for the entire batch.
                for result in pool.imap_unordered(_worker_fn, molecules):
                    if stop_event and stop_event.is_set():
                        # Terminate remaining workers cleanly and exit
                        pool.terminate()
                        pool.join()
                        break
                    result_callback(result, completed, total)
                    completed += 1
        except Exception as e:
            print(f"[PyDes] Multiprocessing error: {e}. Falling back to sequential.")
            traceback.print_exc()
            # Finish remaining molecules sequentially
            for mol in molecules[completed:]:
                if stop_event and stop_event.is_set():
                    break
                result = cls.calculate_single(mol)
                result_callback(result, completed, total)
                completed += 1

        return completed
