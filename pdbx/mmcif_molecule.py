# PyChem-Pro — pdbx/mmcif_molecule.py
# Written for Python 3. No Python 2 compatibility required.
#

"""
pdbx/mmcif_molecule.py
======================
PyChem-Pro bridge between PDBx/mmCIF CIF containers and the
``Molecule / Atom / Bond`` domain model.

Public API
----------
::

    from pdbx import read_mmcif, write_mmcif

    mol = read_mmcif("1abc.cif")          # -> Molecule
    write_mmcif(mol, "out.cif")           # -> None  (also returns str)

    # multi-model files
    mol = read_mmcif("1abc.cif", model=2)

    # string input
    mol = read_mmcif(cif_string, from_string=True)

Design
------
Reading
~~~~~~~
1. Parse the CIF file with ``PdbxReader`` into a list of ``DataContainer``
   objects (the raw container model from the WWPDB pdbx library).
2. ``MmcifExtractor`` picks the first ``DataContainer`` and extracts:
   - ``_atom_site``         -> Atom coordinates (auth_* preferred over label_*)
   - ``_struct_conf``       -> Secondary structure (HELX, STRN -> H / E / C)
   - ``_struct_conn``       -> SSBOND / LINK bonds
   - ``_entity``            -> Molecule name / type
   - ``_entry``             -> Entry ID used as molecule name fallback
3. A ``Molecule`` is assembled and returned.

Writing
~~~~~~~
1. A ``DataContainer`` is built from a ``Molecule``.
2. Categories written:
   - ``_entry``             (id)
   - ``_atom_site``         (all atoms)
   - ``_struct_conn``       (covalent bonds for HETATM, SSBOND for CYS)
   - ``_pdbx_software``     (credit: PyChem-Pro)
3. ``PdbxWriter`` serialises to file.

No new dependencies — uses only stdlib ``io``, ``os``, ``re``, and the
bundled ``pdbx`` package.
"""

from __future__ import annotations

import io
import os
import re
from typing import Dict, List, Optional, Tuple

from pdbx.reader.PdbxReader import PdbxReader
from pdbx.writer.PdbxWriter import PdbxWriter
from pdbx.reader.PdbxContainers import DataContainer, DataCategory


# ──────────────────────────────────────────────────────────────────────────────
# Sentinel for "column not present in this CIF"
# ──────────────────────────────────────────────────────────────────────────────
_MISSING = object()

# Null-value tokens used by the CIF standard
_CIF_NULL = {".", "?", None, ""}

# Standard 20 amino acids + common non-standard residues
_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEC', 'PYL', 'HIP', 'HIE', 'HID', 'CYX', 'CYM',
}

# Nucleic acid residues
_NUCLEOTIDES = {
    'A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT', 'DU',
    'ADE', 'CYT', 'GUA', 'THY', 'URA',
}


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _null(val) -> bool:
    """Return True if *val* is a CIF null/unknown token."""
    return val in _CIF_NULL


def _safe_float(val, default=0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _safe_int(val, default=0) -> int:
    try:
        return int(val)
    except (TypeError, ValueError):
        return default


def _safe_str(val, default="") -> str:
    return default if _null(val) else str(val).strip()


# ──────────────────────────────────────────────────────────────────────────────
# MmcifExtractor — fast-path category accessor
# ──────────────────────────────────────────────────────────────────────────────

class MmcifExtractor:
    """Lightweight accessor for mmCIF categories within a DataContainer.

    Provides index-based column access so that row iteration is O(1) per
    cell rather than O(n_attributes) per call.
    """

    def __init__(self, container: DataContainer):
        self.container = container

    # ── Generic helpers ────────────────────────────────────────────────

    def has_category(self, cat_name: str) -> bool:
        return self.container.exists(cat_name)

    def get_category(self, cat_name: str) -> Optional[DataCategory]:
        return self.container.getObj(cat_name)

    def col_index(self, category: DataCategory, col: str) -> int:
        """Return column index or -1 if the column is not present."""
        return category.getAttributeIndex(col)

    def cell(self, row: list, col_idx: int, default=None):
        """Safe row[col_idx] with null-checking."""
        if col_idx < 0 or col_idx >= len(row):
            return default
        val = row[col_idx]
        return default if _null(val) else val

    # ── _atom_site ─────────────────────────────────────────────────────

    def extract_atom_site(self, model: int = 1) -> List[dict]:
        """Return a list of atom dicts from _atom_site, filtered by *model*.

        Auth columns are preferred over label columns (PDB convention).
        """
        cat = self.get_category("atom_site")
        if cat is None:
            return []

        # Build column-index map once
        cols = {
            name: self.col_index(cat, name)
            for name in [
                "id", "group_PDB", "type_symbol",
                "label_atom_id", "auth_atom_id",
                "label_comp_id", "auth_comp_id",
                "label_asym_id", "auth_asym_id",
                "label_seq_id", "auth_seq_id",
                "pdbx_PDB_ins_code",
                "Cartn_x", "Cartn_y", "Cartn_z",
                "occupancy", "B_iso_or_equiv",
                "pdbx_PDB_model_num",
                "pdbx_formal_charge",
                "label_alt_id",
            ]
        }

        c = cols  # shorthand
        atoms = []

        for row in cat.getRowList():
            # Model filter
            model_val = self.cell(row, c["pdbx_PDB_model_num"])
            if model_val is not None and _safe_int(model_val, 1) != model:
                continue

            # Skip alternate conformers (keep first / blank)
            alt = self.cell(row, c["label_alt_id"])
            if alt is not None and alt not in (".", "A", "1"):
                continue

            # Auth values take priority over label values
            atom_name = (self.cell(row, c["auth_atom_id"])
                         or self.cell(row, c["label_atom_id"])
                         or "")
            res_name = (self.cell(row, c["auth_comp_id"])
                        or self.cell(row, c["label_comp_id"])
                        or "")
            chain_id = (self.cell(row, c["auth_asym_id"])
                        or self.cell(row, c["label_asym_id"])
                        or "A")
            res_seq_raw = (self.cell(row, c["auth_seq_id"])
                           or self.cell(row, c["label_seq_id"])
                           or "0")

            atoms.append({
                "serial":       _safe_int(self.cell(row, c["id"]), 0),
                "group":        _safe_str(self.cell(row, c["group_PDB"]), "ATOM"),
                "element":      _safe_str(self.cell(row, c["type_symbol"]), "X"),
                "pdb_name":     atom_name.strip(),
                "res_name":     res_name.strip(),
                "chain_id":     chain_id.strip(),
                "res_seq":      _safe_int(res_seq_raw, 0),
                "ins_code":     _safe_str(self.cell(row, c["pdbx_PDB_ins_code"]), ""),
                "x":            _safe_float(self.cell(row, c["Cartn_x"]), 0.0),
                "y":            _safe_float(self.cell(row, c["Cartn_y"]), 0.0),
                "z":            _safe_float(self.cell(row, c["Cartn_z"]), 0.0),
                "occupancy":    _safe_float(self.cell(row, c["occupancy"]), 1.0),
                "b_factor":     _safe_float(self.cell(row, c["B_iso_or_equiv"]), 0.0),
                "formal_charge": _safe_int(self.cell(row, c["pdbx_formal_charge"]), 0),
                "is_hetatm":    _safe_str(self.cell(row, c["group_PDB"]), "ATOM") == "HETATM",
            })

        return atoms

    # ── _struct_conf ───────────────────────────────────────────────────

    def extract_struct_conf(self) -> List[dict]:
        """Return secondary structure range records from _struct_conf."""
        cat = self.get_category("struct_conf")
        if cat is None:
            return []

        cols = {n: self.col_index(cat, n) for n in [
            "conf_type_id",
            "beg_auth_asym_id", "beg_auth_seq_id",
            "end_auth_asym_id", "end_auth_seq_id",
        ]}
        c = cols

        records = []
        for row in cat.getRowList():
            conf_type = _safe_str(self.cell(row, c["conf_type_id"]), "").upper()
            if conf_type.startswith("HELX"):
                ss = "H"
            elif conf_type.startswith("STRN"):
                ss = "E"
            else:
                continue

            beg_chain = _safe_str(self.cell(row, c["beg_auth_asym_id"]), "A")
            end_chain = _safe_str(self.cell(row, c["end_auth_asym_id"]), "A")
            beg_seq = _safe_int(self.cell(row, c["beg_auth_seq_id"]), 0)
            end_seq = _safe_int(self.cell(row, c["end_auth_seq_id"]), 0)
            records.append({
                "ss": ss,
                "beg_chain": beg_chain, "beg_seq": beg_seq,
                "end_chain": end_chain, "end_seq": end_seq,
            })

        return records

    # ── _struct_sheet_range ────────────────────────────────────────────

    def extract_struct_sheet_range(self) -> List[dict]:
        """Return beta-sheet ranges from _struct_sheet_range."""
        cat = self.get_category("struct_sheet_range")
        if cat is None:
            return []

        cols = {n: self.col_index(cat, n) for n in [
            "beg_auth_asym_id", "beg_auth_seq_id",
            "end_auth_asym_id", "end_auth_seq_id",
        ]}
        c = cols

        records = []
        for row in cat.getRowList():
            records.append({
                "ss": "E",
                "beg_chain": _safe_str(self.cell(row, c["beg_auth_asym_id"]), "A"),
                "beg_seq":   _safe_int(self.cell(row, c["beg_auth_seq_id"]), 0),
                "end_chain": _safe_str(self.cell(row, c["end_auth_asym_id"]), "A"),
                "end_seq":   _safe_int(self.cell(row, c["end_auth_seq_id"]), 0),
            })
        return records

    # ── _struct_conn ───────────────────────────────────────────────────

    def extract_struct_conn(self) -> List[dict]:
        """Return inter-residue bond records from _struct_conn.

        Captures SSBOND (disulfide), LINK (metal coordination, covalent), etc.
        """
        cat = self.get_category("struct_conn")
        if cat is None:
            return []

        cols = {n: self.col_index(cat, n) for n in [
            "conn_type_id",
            "ptnr1_auth_asym_id", "ptnr1_auth_comp_id",
            "ptnr1_auth_seq_id",  "ptnr1_label_atom_id",
            "ptnr2_auth_asym_id", "ptnr2_auth_comp_id",
            "ptnr2_auth_seq_id",  "ptnr2_label_atom_id",
        ]}
        c = cols

        records = []
        for row in cat.getRowList():
            conn_type = _safe_str(self.cell(row, c["conn_type_id"]), "").lower()
            records.append({
                "conn_type": conn_type,
                "chain1":    _safe_str(self.cell(row, c["ptnr1_auth_asym_id"]), ""),
                "res_name1": _safe_str(self.cell(row, c["ptnr1_auth_comp_id"]), ""),
                "seq1":      _safe_int(self.cell(row, c["ptnr1_auth_seq_id"]), 0),
                "atom1":     _safe_str(self.cell(row, c["ptnr1_label_atom_id"]), "SG"),
                "chain2":    _safe_str(self.cell(row, c["ptnr2_auth_asym_id"]), ""),
                "res_name2": _safe_str(self.cell(row, c["ptnr2_auth_comp_id"]), ""),
                "seq2":      _safe_int(self.cell(row, c["ptnr2_auth_seq_id"]), 0),
                "atom2":     _safe_str(self.cell(row, c["ptnr2_label_atom_id"]), "SG"),
            })
        return records

    # ── _entity ────────────────────────────────────────────────────────

    def extract_entity_info(self) -> Dict[str, str]:
        """Return a dict mapping entity_id -> description from _entity."""
        cat = self.get_category("entity")
        if cat is None:
            return {}

        id_idx = self.col_index(cat, "id")
        desc_idx = self.col_index(cat, "pdbx_description")
        type_idx = self.col_index(cat, "type")

        info = {}
        for row in cat.getRowList():
            eid = self.cell(row, id_idx, "")
            desc = self.cell(row, desc_idx, "")
            etype = self.cell(row, type_idx, "")
            if eid:
                info[str(eid)] = {"description": str(desc or ""), "type": str(etype or "")}
        return info

    # ── _entry ─────────────────────────────────────────────────────────

    def get_entry_id(self) -> str:
        """Return the PDB entry ID from _entry.id, or empty string."""
        cat = self.get_category("entry")
        if cat is None:
            return ""
        id_idx = self.col_index(cat, "id")
        if cat.getRowCount() > 0:
            return self.cell(cat.getRow(0), id_idx, "") or ""
        return ""


# ──────────────────────────────────────────────────────────────────────────────
# read_mmcif
# ──────────────────────────────────────────────────────────────────────────────

def read_mmcif(source, from_string: bool = False, model: int = 1):
    """Read a PDBx/mmCIF file or string and return a ``Molecule``.

    Args:
        source:      File path (str) **or** CIF text (str) when
                     *from_string* is True.
        from_string: If True, treat *source* as raw CIF text rather
                     than a file path.
        model:       Model number to load (default=1). Relevant for
                     NMR ensembles that contain multiple models.

    Returns:
        A ``Molecule`` instance with atoms, bonds, secondary structure,
        and metadata stored in ``molecule.properties``.

    Raises:
        FileNotFoundError: If *source* is a path and does not exist.
        ValueError:        If no ``data_`` block is found in the CIF.
    """
    # ── Lazy domain imports (avoids circular imports at module level) ──
    from src.core.domain.models.molecule import Molecule
    from src.core.domain.models.atom import Atom
    from src.core.domain.models.bond import Bond, BondType

    # ── Parse CIF ─────────────────────────────────────────────────────
    container_list: List[DataContainer] = []

    if from_string:
        fh = io.StringIO(source)
        reader = PdbxReader(fh)
        reader.read(container_list)
    else:
        if not os.path.isfile(source):
            raise FileNotFoundError(f"mmCIF file not found: {source}")
        with open(source, "r", encoding="utf-8", errors="replace") as fh:
            reader = PdbxReader(fh)
            reader.read(container_list)

    # Pick the first DataContainer
    data_block = None
    for c in container_list:
        if hasattr(c, "getGlobal") and not c.getGlobal():
            data_block = c
            break

    if data_block is None:
        raise ValueError("No data_ block found in CIF input")

    extractor = MmcifExtractor(data_block)

    # ── Extract entry / entity info ────────────────────────────────────
    entry_id = extractor.get_entry_id()
    entity_info = extractor.extract_entity_info()

    # Determine molecule name
    mol_name = entry_id or (os.path.splitext(os.path.basename(source))[0]
                            if not from_string else "mmcif_molecule")

    mol = Molecule(name=mol_name)
    mol.properties["entry_id"] = entry_id
    mol.properties["entity_info"] = entity_info
    mol.properties["source_format"] = "mmcif"
    mol.properties["model"] = model

    # ── Extract atoms ──────────────────────────────────────────────────
    atom_dicts = extractor.extract_atom_site(model=model)

    if not atom_dicts:
        # Return empty molecule rather than crashing
        return mol

    mol.begin_bulk_load()

    # key: (chain_id, res_seq, ins_code, pdb_name) -> atom index
    atom_key_to_idx: Dict[Tuple, int] = {}

    for ad in atom_dicts:
        atom = Atom(ad["element"] or "X")
        atom.x = ad["x"]
        atom.y = ad["y"]
        atom.z = ad["z"]
        atom.pdb_name = ad["pdb_name"]
        atom.res_name = ad["res_name"]
        atom.chain_id = ad["chain_id"]
        atom.res_seq = ad["res_seq"]
        atom.b_factor = ad["b_factor"]
        atom.is_hetatm = ad["is_hetatm"]
        atom.formal_charge = ad["formal_charge"]

        # Store occupancy in a dict we'll attach later
        idx = mol.add_atom(atom)
        key = (ad["chain_id"], ad["res_seq"], ad["ins_code"], ad["pdb_name"])
        atom_key_to_idx[key] = idx

    mol.properties["occupancies"] = {
        i: ad["occupancy"]
        for i, ad in enumerate(atom_dicts)
    }

    # ── Secondary structure ────────────────────────────────────────────
    ss_records = extractor.extract_struct_conf()
    ss_records += extractor.extract_struct_sheet_range()

    # Build lookup: chain -> list of (start, end, ss_code)
    ss_lookup: Dict[str, List[Tuple[int, int, str]]] = {}
    for rec in ss_records:
        chain = rec["beg_chain"]
        ss_lookup.setdefault(chain, []).append(
            (rec["beg_seq"], rec["end_seq"], rec["ss"])
        )

    for atom in mol.atoms:
        ss = "C"  # default: coil
        chain = getattr(atom, "chain_id", "")
        rseq = getattr(atom, "res_seq", 0)
        if chain in ss_lookup:
            for start, end, ss_code in ss_lookup[chain]:
                if start <= rseq <= end:
                    ss = ss_code
                    break
        atom.ss_type = ss

    # ── Struct_conn bonds (SSBOND / LINK) ─────────────────────────────
    conn_records = extractor.extract_struct_conn()
    added_bonds: set = set()

    for conn in conn_records:
        key1 = (conn["chain1"], conn["seq1"], "", conn["atom1"])
        key2 = (conn["chain2"], conn["seq2"], "", conn["atom2"])
        idx1 = atom_key_to_idx.get(key1)
        idx2 = atom_key_to_idx.get(key2)

        if idx1 is None or idx2 is None:
            continue

        bond_key = tuple(sorted((idx1, idx2)))
        if bond_key in added_bonds:
            continue
        added_bonds.add(bond_key)

        conn_type = conn["conn_type"]
        if "disulf" in conn_type or "ssbond" in conn_type:
            bt = BondType.SINGLE
        elif "metal" in conn_type:
            bt = BondType.SINGLE
        else:
            bt = BondType.SINGLE

        mol.add_bond(idx1, idx2, bt)

    mol.end_bulk_load()

    # ── Auto-bond ──────────────────────────────────────────────────────
    from src.features.io.loaders.file_reader import _auto_bond_pdb
    if len(mol.atoms) > 0 and len(mol.bonds) < len(mol.atoms) * 0.8:
        _auto_bond_pdb(mol)

    # ── Derive is_protein / is_nucleic flags ───────────────────────────
    res_names = {a.res_name for a in mol.atoms}
    mol.properties["is_protein"] = bool(res_names & _AMINO_ACIDS)
    mol.properties["is_nucleic"] = bool(res_names & _NUCLEOTIDES)
    mol.properties["helix_ranges"] = [
        (r["beg_chain"], r["beg_seq"], r["end_seq"])
        for r in ss_records if r["ss"] == "H"
    ]
    mol.properties["sheet_ranges"] = [
        (r["beg_chain"], r["beg_seq"], r["end_seq"])
        for r in ss_records if r["ss"] == "E"
    ]

    return mol


# ──────────────────────────────────────────────────────────────────────────────
# write_mmcif
# ──────────────────────────────────────────────────────────────────────────────

def write_mmcif(molecule, filepath: Optional[str] = None,
                entry_id: Optional[str] = None) -> str:
    """Write a ``Molecule`` to PDBx/mmCIF format.

    All atoms are written to the ``_atom_site`` category.  Covalent bonds
    between HETATM atoms are recorded in ``_struct_conn``.

    Args:
        molecule:  The Molecule to serialise.
        filepath:  If given, write to this path (UTF-8).  If None, return
                   the CIF text as a string.
        entry_id:  PDB-style entry identifier (e.g. ``"1ABC"``).  Defaults to
                   ``molecule.name`` (truncated to 4 chars if longer).

    Returns:
        The CIF text as a string.
    """
    from src.core.domain.models.bond import BondType

    eid = entry_id or (molecule.name[:4].upper() if molecule.name else "XXXX")

    container = DataContainer(eid)

    # ── _entry ─────────────────────────────────────────────────────────
    entry_cat = DataCategory("entry")
    entry_cat.appendAttribute("id")
    entry_cat.append([eid])
    container.append(entry_cat)

    # ── _atom_site ──────────────────────────────────────────────────────
    atom_cols = [
        "group_PDB", "id", "type_symbol",
        "label_atom_id",  "label_alt_id",
        "label_comp_id",  "label_asym_id", "label_entity_id", "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x", "Cartn_y", "Cartn_z",
        "occupancy", "B_iso_or_equiv",
        "pdbx_formal_charge",
        "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id",
        "pdbx_PDB_model_num",
    ]

    atom_cat = DataCategory("atom_site")
    for col in atom_cols:
        atom_cat.appendAttribute(col)

    for i, atom in enumerate(molecule.atoms):
        serial = i + 1
        group = "HETATM" if getattr(atom, "is_hetatm", False) else "ATOM"
        element = atom.symbol.upper()

        # PDB atom name
        pdb_name = getattr(atom, "pdb_name", None) or element
        res_name = getattr(atom, "res_name", None) or "LIG"
        chain_id = getattr(atom, "chain_id", None) or "A"
        res_seq = getattr(atom, "res_seq", None)
        if res_seq is None:
            res_seq = 1
        b_factor = getattr(atom, "b_factor", None)
        b_factor = 0.0 if b_factor is None else b_factor

        x = atom.x if atom.x is not None else 0.0
        y = atom.y if atom.y is not None else 0.0
        z = atom.z if atom.z is not None else 0.0

        formal_charge = getattr(atom, "formal_charge", 0) or 0

        row = [
            group,                     # group_PDB
            str(serial),               # id
            element,                   # type_symbol
            pdb_name,                  # label_atom_id
            ".",                       # label_alt_id
            res_name,                  # label_comp_id
            chain_id,                  # label_asym_id
            "1",                       # label_entity_id
            str(res_seq),              # label_seq_id
            "?",                       # pdbx_PDB_ins_code
            f"{x:.3f}",               # Cartn_x
            f"{y:.3f}",               # Cartn_y
            f"{z:.3f}",               # Cartn_z
            "1.00",                    # occupancy
            f"{b_factor:.2f}",        # B_iso_or_equiv
            str(formal_charge) if formal_charge != 0 else "?",   # pdbx_formal_charge
            str(res_seq),              # auth_seq_id
            res_name,                  # auth_comp_id
            chain_id,                  # auth_asym_id
            pdb_name,                  # auth_atom_id
            "1",                       # pdbx_PDB_model_num
        ]
        atom_cat.append(row)

    container.append(atom_cat)

    # ── _struct_conn (for HETATM covalent bonds + disulfides) ──────────
    conn_bonds = []
    for bond in molecule.bonds:
        a1 = molecule.atoms[bond.begin_atom_idx]
        a2 = molecule.atoms[bond.end_atom_idx]

        is_ssbond = (
            getattr(a1, "res_name", "") == "CYS"
            and getattr(a2, "res_name", "") == "CYS"
            and getattr(a1, "pdb_name", "") == "SG"
            and getattr(a2, "pdb_name", "") == "SG"
        )

        # Only emit struct_conn for HETATM or disulfide bonds
        if getattr(a1, "is_hetatm", False) or getattr(a2, "is_hetatm", False) or is_ssbond:
            conn_bonds.append((bond, a1, a2, is_ssbond))

    if conn_bonds:
        conn_cols = [
            "id", "conn_type_id",
            "ptnr1_label_comp_id", "ptnr1_label_asym_id",
            "ptnr1_label_seq_id",  "ptnr1_label_atom_id",
            "ptnr1_auth_comp_id",  "ptnr1_auth_asym_id",
            "ptnr1_auth_seq_id",   "ptnr1_label_atom_id",
            "ptnr2_label_comp_id", "ptnr2_label_asym_id",
            "ptnr2_label_seq_id",  "ptnr2_label_atom_id",
            "ptnr2_auth_comp_id",  "ptnr2_auth_asym_id",
            "ptnr2_auth_seq_id",   "ptnr2_label_atom_id",
        ]
        conn_cat = DataCategory("struct_conn")
        for col in conn_cols:
            conn_cat.appendAttribute(col)

        for idx, (bond, a1, a2, is_ssbond) in enumerate(conn_bonds, start=1):
            conn_type = "disulf" if is_ssbond else "covale"
            pdb1 = getattr(a1, "pdb_name", a1.symbol)
            pdb2 = getattr(a2, "pdb_name", a2.symbol)
            res1 = getattr(a1, "res_name", "LIG")
            res2 = getattr(a2, "res_name", "LIG")
            chain1 = getattr(a1, "chain_id", "A")
            chain2 = getattr(a2, "chain_id", "A")
            seq1 = str(getattr(a1, "res_seq", 1) or 1)
            seq2 = str(getattr(a2, "res_seq", 1) or 1)

            row = [
                f"conn{idx}", conn_type,
                res1, chain1, seq1, pdb1, res1, chain1, seq1, pdb1,
                res2, chain2, seq2, pdb2, res2, chain2, seq2, pdb2,
            ]
            conn_cat.append(row)

        container.append(conn_cat)

    # ── _software ───────────────────────────────────────────────────────
    sw_cat = DataCategory("software")
    sw_cat.appendAttribute("pdbx_ordinal")
    sw_cat.appendAttribute("name")
    sw_cat.appendAttribute("version")
    sw_cat.appendAttribute("classification")
    sw_cat.appendAttribute("date")
    sw_cat.append(["1", "PyChem-Pro", "1.0.0", "data collection", "?"])
    container.append(sw_cat)

    # ── Serialise ───────────────────────────────────────────────────────
    buf = io.StringIO()
    writer = PdbxWriter(buf)
    writer.write([container])
    cif_text = buf.getvalue()

    if filepath:
        with open(filepath, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(cif_text)

    return cif_text


# ──────────────────────────────────────────────────────────────────────────────
# Convenience: read multiple models from a multi-model CIF
# ──────────────────────────────────────────────────────────────────────────────

def read_mmcif_models(source, from_string: bool = False) -> List:
    """Read all models from a multi-model mmCIF file (e.g. NMR ensembles).

    Args:
        source:      File path or CIF text string.
        from_string: If True, treat *source* as raw CIF text.

    Returns:
        A list of Molecule objects, one per model number found.
    """
    from src.core.domain.models.molecule import Molecule

    # First pass: discover which model numbers are present
    container_list: List[DataContainer] = []
    if from_string:
        fh = io.StringIO(source)
        PdbxReader(fh).read(container_list)
    else:
        if not os.path.isfile(source):
            raise FileNotFoundError(f"mmCIF file not found: {source}")
        with open(source, "r", encoding="utf-8", errors="replace") as fh:
            PdbxReader(fh).read(container_list)

    data_block = next(
        (c for c in container_list
         if hasattr(c, "getGlobal") and not c.getGlobal()),
        None
    )
    if data_block is None:
        return []

    extractor = MmcifExtractor(data_block)
    cat = extractor.get_category("atom_site")
    if cat is None:
        return []

    model_col = extractor.col_index(cat, "pdbx_PDB_model_num")
    model_nums_seen = set()
    for row in cat.getRowList():
        m = extractor.cell(row, model_col)
        if m is not None:
            model_nums_seen.add(_safe_int(m, 1))

    if not model_nums_seen:
        model_nums_seen = {1}

    molecules = []
    for model_num in sorted(model_nums_seen):
        mol = read_mmcif(source, from_string=from_string, model=model_num)
        mol.properties["model"] = model_num
        molecules.append(mol)

    return molecules
