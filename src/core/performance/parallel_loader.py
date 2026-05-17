"""
Parallel File Loader — Multiprocessing-optimized molecular file loading.

Uses 50% of available CPU cores to speed up file parsing for large molecules.
Implements chunked reading and parallel processing for MOL, MOL2, SDF, and PDB files.

Delegates all multiprocessing to the centralized ParallelExecutor
(src/core/parallel.py) so every subsystem shares one process pool
capped at 50 % of CPU cores.
"""

import os
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any
import logging

from src.core.parallel import ParallelExecutor

# Set up logging
logger = logging.getLogger(__name__)


class ParallelFileLoader:
    """
    High-performance file loader using multiprocessing for large files.

    Automatically detects file size and uses parallel processing for files
    above a threshold (100KB default). Smaller files use standard loading.

    Attributes:
        executor: ParallelExecutor instance used for multiprocessing
        parallel_threshold_kb: File size threshold for parallel processing
    """

    def __init__(self, executor: Optional[ParallelExecutor] = None,
                 parallel_threshold_kb: float = 100.0,
                 timeout_seconds: float = 30.0,
                 # Legacy kwargs kept for backward compatibility
                 workers: Optional[int] = None,
                 use_parallel: Optional[bool] = None):
        """
        Initialize the parallel file loader.

        Args:
            executor: ParallelExecutor instance (creates a default one if None)
            parallel_threshold_kb: Minimum file size (KB) to use parallel processing
            timeout_seconds: Maximum time to wait for parallel processing
            workers: *Deprecated* — ignored when *executor* is supplied.
                     Kept for backward compatibility; if no executor is given
                     and workers is set, a ParallelExecutor with that many
                     workers is created.
            use_parallel: *Deprecated* — set to False to disable parallelism.
        """
        self.parallel_threshold_kb = parallel_threshold_kb
        self.timeout_seconds = timeout_seconds
        self._import_cache = {}

        # Determine if parallel processing is explicitly disabled
        if use_parallel is False:
            self._use_parallel = False
            self.executor = executor or ParallelExecutor(max_workers=1)
            logger.info("ParallelFileLoader: parallel processing explicitly disabled")
            return

        self._use_parallel = True

        # Build or accept the executor
        if executor is not None:
            self.executor = executor
        elif workers is not None:
            self.executor = ParallelExecutor(max_workers=workers)
        else:
            self.executor = ParallelExecutor()

        logger.info(
            f"ParallelFileLoader initialized: workers={self.executor.num_workers}, "
            f"threshold={parallel_threshold_kb}KB"
        )
    
    def load_file(self, filepath: str) -> Any:
        """
        Load a molecular file with automatic parallel processing for large files.
        
        Args:
            filepath: Path to the molecular file (MOL, MOL2, SDF, PDB)
            
        Returns:
            Molecule object
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is unsupported
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        file_size_kb = filepath.stat().st_size / 1024
        ext = filepath.suffix.lower()
        
        # Use parallel loading for large files (if available)
        if (file_size_kb > self.parallel_threshold_kb and
            self.executor.num_workers > 1 and
            self._use_parallel):
            logger.info(f"Using parallel loading for {filepath.name} ({file_size_kb:.1f}KB)")
            try:
                return self._load_parallel(filepath, ext)
            except Exception as e:
                logger.warning(f"Parallel loading failed: {e}. Falling back to standard loading.")
                return self._load_standard(filepath, ext)
        else:
            if self.executor.num_workers == 1:
                logger.info(f"Using standard loading for {filepath.name} ({file_size_kb:.1f}KB)")
            return self._load_standard(filepath, ext)
    
    def _load_standard(self, filepath: Path, ext: str) -> Any:
        """Standard single-threaded file loading."""
        from src.features.io.loaders.file_reader import read_mol, read_sdf, read_mol2, read_pdb
        
        if ext in ('.mol',):
            return read_mol(str(filepath))
        elif ext in ('.sdf', '.sd'):
            return read_sdf(str(filepath))
        elif ext in ('.mol2',):
            return read_mol2(str(filepath))
        elif ext in ('.pdb', '.ent'):
            return read_pdb(str(filepath))
        else:
            raise ValueError(f"Unsupported file format: {ext}")
    
    def _load_parallel(self, filepath: Path, ext: str) -> Any:
        """
        Parallel file loading using multiprocessing.
        
        For PDB files: Parallel atom record processing
        For MOL2 files: Parallel section parsing
        For SDF files: Parallel molecule block processing
        """
        if ext in ('.pdb', '.ent'):
            return self._load_pdb_parallel(filepath)
        elif ext in ('.mol2',):
            return self._load_mol2_parallel(filepath)
        elif ext in ('.sdf', '.sd'):
            return self._load_sdf_parallel(filepath)
        else:
            # Fall back to standard loading for unsupported formats
            return self._load_standard(filepath, ext)
    
    def _load_pdb_parallel(self, filepath: Path) -> Any:
        """
        Load PDB file using parallel processing.
        
        Strategy:
        1. First pass: Count and partition atoms into chunks
        2. Process chunks in parallel to extract atom data
        3. Merge results and build bonds
        """
        from src.core.domain.models.molecule import Molecule
        from src.core.domain.models.atom import Atom
        from src.core.domain.models.bond import Bond, BondType
        
        # Read file content
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Find ATOM/HETATM records and partition into chunks
        atom_records = []
        for i, line in enumerate(lines):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_records.append((i, line))
        
        if not atom_records:
            raise ValueError("No ATOM/HETATM records found in PDB file")
        
        # Partition atoms into chunks for parallel processing
        num_workers = self.executor.num_workers
        chunk_size = max(1, len(atom_records) // num_workers)
        chunks = [
            atom_records[i:i + chunk_size]
            for i in range(0, len(atom_records), chunk_size)
        ]

        logger.info(f"Processing {len(atom_records)} atoms in {len(chunks)} chunks")

        # Process chunks via centralized ParallelExecutor
        results = self.executor.map(
            _parse_pdb_chunk, chunks, timeout=self.timeout_seconds
        )

        all_atoms = []
        for chunk_atoms in results:
            all_atoms.extend(chunk_atoms)
        
        # Sort by original PDB serial number to maintain order
        all_atoms.sort(key=lambda x: x[0])  # x[0] is serial number
        
        # Build molecule
        mol = Molecule(name=filepath.stem)
        serial_to_idx = {}
        mol.begin_bulk_load()
        for serial, d in all_atoms:
            a = Atom(d['symbol'])
            a.x, a.y, a.z = d['x'], d['y'], d['z']
            a.pdb_name, a.res_name, a.chain_id, a.res_seq = d['name'], d['res_name'], d['chain_id'], d['res_seq']
            a.b_factor, a.is_hetatm = d['b_factor'], d['is_hetatm']
            
            idx = mol.add_atom(a)
            serial_to_idx[serial] = idx
        
        # Parse CONECT records
        self._parse_conect_records(lines, mol, serial_to_idx)
        
        # Auto-bond if needed
        from src.features.io.loaders.file_reader import _auto_bond_pdb
        if len(mol.atoms) > 0 and len(mol.bonds) < len(mol.atoms) * 0.8:
            _auto_bond_pdb(mol)
            
        mol.end_bulk_load()
        
        mol.properties['source'] = 'file'
        mol.properties['source_format'] = 'pdb'
        
        logger.info(f"PDB loaded: {len(mol.atoms)} atoms, {len(mol.bonds)} bonds")
        return mol
    
    def _parse_conect_records(self, lines: List[str], mol: Any, serial_to_idx: Dict[int, int]):
        """Parse CONECT records to build bonds."""
        from src.core.domain.models.bond import BondType
        
        added_bonds = set()
        for line in lines:
            if line.startswith('CONECT'):
                try:
                    parts = line[6:].split()
                    if len(parts) >= 2:
                        origin = int(parts[0])
                        for target_str in parts[1:]:
                            target = int(target_str)
                            if origin in serial_to_idx and target in serial_to_idx:
                                idx1 = serial_to_idx[origin]
                                idx2 = serial_to_idx[target]
                                bond_key = tuple(sorted((idx1, idx2)))
                                if bond_key not in added_bonds:
                                    mol.add_bond(idx1, idx2, BondType.SINGLE)
                                    added_bonds.add(bond_key)
                except (ValueError, IndexError):
                    continue
    
    def _load_mol2_parallel(self, filepath: Path) -> Any:
        """
        Load MOL2 file using parallel processing.
        
        Strategy:
        1. Parse header to find section boundaries
        2. Process ATOM and BOND sections in parallel
        3. Merge results
        """
        from src.core.domain.models.molecule import Molecule
        from src.core.domain.models.atom import Atom
        from src.core.domain.models.bond import BondType
        
        with open(filepath, 'r') as f:
            content = f.read()
        
        lines = content.split('\n')
        
        # Find section boundaries
        sections = {}
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith('@<TRIPOS>'):
                section_name = stripped.replace('@<TRIPOS>', '')
                sections[section_name] = i
        
        # Get molecule name
        mol_name = ""
        if 'MOLECULE' in sections:
            mol_start = sections['MOLECULE']
            if mol_start + 1 < len(lines):
                mol_name = lines[mol_start + 1].strip()
        
        mol = Molecule(name=mol_name)
        mol.properties['source'] = 'file'
        
        # Parse ATOM section in parallel if large
        if 'ATOM' in sections:
            atom_start = sections['ATOM'] + 1
            atom_end = len(lines)
            
            # Find end of atom section
            for i in range(atom_start, len(lines)):
                if lines[i].strip().startswith('@<TRIPOS>'):
                    atom_end = i
                    break
            
            atom_lines = lines[atom_start:atom_end]
            
            num_workers = self.executor.num_workers
            if len(atom_lines) > 500 and num_workers > 1:
                # Parallel atom parsing via centralized ParallelExecutor
                chunk_size = max(1, len(atom_lines) // num_workers)
                chunks = [
                    atom_lines[i:i + chunk_size]
                    for i in range(0, len(atom_lines), chunk_size)
                ]

                results = self.executor.map(
                    _parse_mol2_atom_chunk_simple, chunks,
                    timeout=self.timeout_seconds
                )

                all_atoms = []
                for chunk_atoms in results:
                    all_atoms.extend(chunk_atoms)

                # Sort by atom_id and add to molecule
                all_atoms.sort(key=lambda x: x[0])
                for _, atom in all_atoms:
                    mol.add_atom(atom)
            else:
                # Standard parsing for small files
                for i, line in enumerate(atom_lines):
                    atom = _parse_mol2_atom_line(line, i)
                    if atom:
                        mol.add_atom(atom)
        
        # Parse BOND section
        if 'BOND' in sections:
            bond_start = sections['BOND'] + 1
            bond_end = len(lines)
            
            for i in range(bond_start, len(lines)):
                if lines[i].strip().startswith('@<TRIPOS>'):
                    bond_end = i
                    break
            
            for i in range(bond_start, bond_end):
                line = lines[i].strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                
                a1 = int(parts[1]) - 1
                a2 = int(parts[2]) - 1
                bond_type_str = parts[3].lower()
                
                if bond_type_str == '2':
                    bt = BondType.DOUBLE
                elif bond_type_str == '3':
                    bt = BondType.TRIPLE
                elif bond_type_str in ('ar', 'am'):
                    bt = BondType.AROMATIC
                else:
                    bt = BondType.SINGLE
                
                if 0 <= a1 < len(mol.atoms) and 0 <= a2 < len(mol.atoms):
                    mol.add_bond(a1, a2, bt)
        
        logger.info(f"MOL2 loaded: {len(mol.atoms)} atoms, {len(mol.bonds)} bonds")
        return mol
    
    def _load_sdf_parallel(self, filepath: Path) -> Any:
        """
        Load SDF file with parallel molecule processing for multi-molecule files.
        """
        from src.features.io.loaders.file_reader import read_sdf
        
        # For now, use standard loading (SDF files are typically small)
        # Future: implement parallel molecule block processing
        return read_sdf(str(filepath))


# Worker functions for multiprocessing (must be module-level for pickle)
def _parse_pdb_chunk(chunk: List[Tuple[int, str]]) -> List[Tuple[int, Dict]]:
    """Parse a chunk of PDB ATOM/HETATM records."""
    atoms = []
    
    for _, line in chunk:
        try:
            serial = int(line[6:11])
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_seq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            # Extract element
            element = ''
            if len(line) >= 78:
                element = line[76:78].strip()
            if not element:
                name_clean = atom_name.lstrip('0123456789')
                if len(name_clean) >= 2:
                    first_two = name_clean[:2].upper()
                    if first_two in ('CL', 'BR', 'FE', 'ZN', 'MG', 'CA', 'NA', 'MN', 'CU', 'CO', 'NI'):
                        element = first_two
                    else:
                        element = name_clean[0].upper()
                elif name_clean:
                    element = name_clean[0].upper()
            
            b_factor = 0.0
            if len(line) >= 66:
                try:
                    b_factor = float(line[60:66])
                except ValueError:
                    pass
            
            is_hetatm = line.startswith('HETATM')
            
            atom_data = {
                'name': atom_name,
                'res_name': res_name,
                'chain_id': chain_id,
                'res_seq': res_seq,
                'x': x,
                'y': y,
                'z': z,
                'symbol': element or 'C',
                'b_factor': b_factor,
                'is_hetatm': is_hetatm
            }
            
            atoms.append((serial, atom_data))
        except (ValueError, IndexError):
            continue
    
    return atoms


def _parse_mol2_atom_line(line: str, index: int) -> Optional[Tuple[int, Any]]:
    """Parse a single MOL2 ATOM record."""
    from src.core.domain.models.atom import Atom
    
    line = line.strip()
    if not line:
        return None
    
    parts = line.split()
    if len(parts) < 6:
        return None
    
    try:
        atom_id = int(parts[0])
        atom_name = parts[1]
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        atom_type = parts[5]
        
        # Extract element symbol
        symbol = atom_type.split('.')[0]
        if symbol.upper() in ('CL', 'BR'):
            symbol = symbol[0].upper() + symbol[1].lower()
        elif len(symbol) == 1:
            symbol = symbol.upper()
        else:
            symbol = symbol[0].upper() + symbol[1:].lower()
        
        atom = Atom(symbol)
        atom.x = x
        atom.y = y
        atom.z = z
        atom.sybyl_type = atom_type
        
        # Parse partial charge
        if len(parts) >= 9:
            try:
                atom.partial_charge = float(parts[8])
            except ValueError:
                pass
        
        return (atom_id, atom)
    except (ValueError, IndexError):
        return None


def _parse_mol2_atom_chunk(lines: List[str], base_index: int) -> List[Tuple[int, Any]]:
    """Parse a chunk of MOL2 ATOM records (legacy two-arg form)."""
    atoms = []
    for i, line in enumerate(lines):
        result = _parse_mol2_atom_line(line, base_index + i)
        if result:
            atoms.append(result)
    return atoms


def _parse_mol2_atom_chunk_simple(lines: List[str]) -> List[Tuple[int, Any]]:
    """Parse a chunk of MOL2 ATOM records (single-arg for ParallelExecutor.map)."""
    atoms = []
    for i, line in enumerate(lines):
        result = _parse_mol2_atom_line(line, i)
        if result:
            atoms.append(result)
    return atoms
