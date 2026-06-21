"""
Loader service — wraps existing file readers with ParallelExecutor.

Implements the ILoader protocol: can_load / load / supported_extensions.
Large files are automatically routed through ParallelFileLoader for
chunked multiprocess parsing; small files use the standard readers.
"""
import os
from typing import Optional, List

from src.core.parallel import ParallelExecutor


class LoaderService:
    """File loading service with parallel support.

    Usage:
        svc = LoaderService()            # default executor
        mol = svc.load("protein.pdb")    # auto-selects reader
    """

    SUPPORTED: List[str] = [".pdb", ".ent", ".mol", ".mol2", ".sdf", ".sd",
                             ".cif", ".mmcif"]

    def __init__(self, executor: Optional[ParallelExecutor] = None):
        self.executor = executor or ParallelExecutor()

    # ------------------------------------------------------------------
    # ILoader protocol
    # ------------------------------------------------------------------

    def supported_extensions(self) -> List[str]:
        """Return list of file extensions this service can handle."""
        return list(self.SUPPORTED)

    def can_load(self, path: str) -> bool:
        """Return True if *path* has a supported extension."""
        ext = os.path.splitext(path)[1].lower()
        return ext in self.SUPPORTED

    def load(self, path: str, parallel: bool = True):
        """Load a molecular file, returning a Molecule.

        For large files (>100 KB) parallel chunk parsing is used
        automatically when *parallel* is True.

        Args:
            path: Path to a molecular file.
            parallel: Allow parallel processing for large files.

        Returns:
            Molecule instance.

        Raises:
            FileNotFoundError: If *path* does not exist.
            ValueError: If the file extension is not supported.
        """
        if not os.path.isfile(path):
            raise FileNotFoundError(f"File not found: {path}")

        ext = os.path.splitext(path)[1].lower()

        if ext not in self.SUPPORTED:
            raise ValueError(f"Unsupported file format: {ext}")

        if parallel:
            # Delegate to ParallelFileLoader which handles size-based
            # decisions and falls back to standard readers internally.
            from src.core.performance.parallel_loader import ParallelFileLoader

            loader = ParallelFileLoader(executor=self.executor)
            return loader.load_file(path)

        # Non-parallel path — call readers directly.
        return self._load_standard(path, ext)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _load_standard(path: str, ext: str):
        """Single-threaded load via existing file_reader functions."""
        from src.features.io.loaders.file_reader import (
            read_mol, read_sdf, read_mol2, read_pdb,
        )

        if ext == ".pdb" or ext == ".ent":
            return read_pdb(path)
        elif ext == ".mol":
            return read_mol(path)
        elif ext == ".mol2":
            return read_mol2(path)
        elif ext in (".sdf", ".sd"):
            return read_sdf(path)
        elif ext in (".cif", ".mmcif"):
            from pdbx.mmcif_molecule import read_mmcif
            return read_mmcif(path)
        else:
            raise ValueError(f"Unsupported file format: {ext}")
