"""
MMFF94 parameter loader -- type-keyed, sourced from Jmol's mmff94.par.txt.

Public API: `get_bond_params`, `get_angle_params`, `get_sb_params`,
`get_torsion_params`, `get_oop_params`, `get_vdw_params`, `get_bci`,
`get_pbci`, `get_atom_props`. All keyed by MMFF94 atom type numbers (int).

Loading is lazy: the first call triggers a parse + pickle-cache write;
subsequent calls reuse the in-memory dicts.
"""
from __future__ import annotations

import os
import pickle
import re
import threading
from pathlib import Path
from typing import Dict, Tuple, Optional

_DATA_DIR = Path(__file__).parent / "data"
_PARAM_FILE = _DATA_DIR / "mmff94.par.txt"
_ATOM_TYPES_FILE = _DATA_DIR / "mmff94_atom_types.txt"
_CACHE_DIR = _DATA_DIR / "_cache"
_CACHE_FILE = _CACHE_DIR / "mmff94_params.pkl"

_load_lock = threading.Lock()
_PARAMS: Optional[Dict[str, dict]] = None

# Sentinel for the lazy-loaded atom-type properties dict.
_ATOM_PROPS: Optional[Dict[int, dict]] = None
_atom_props_lock = threading.Lock()


# Section header pattern. Matches lines like:
#   "8.  MMFFANG.PAR: This file supplies parameters for angle-bending..."
#   "6.  MMFFBOND.PAR: ..."
_SECTION_HEADER = re.compile(
    r"^\s*\d+\.\s+(MMFF[A-Z]+)\.PAR\b", re.IGNORECASE
)

# Whitelist of sections we care about. Jmol's mmff94.par.txt also contains
# MMFFBNDK (empirical-rule bond-length defaults) which we ignore in Phase 1.
_EXPECTED_SECTIONS = frozenset({
    "MMFFPBCI", "MMFFCHG", "MMFFANG", "MMFFBOND",
    "MMFFOOP", "MMFFSTBN", "MMFFDFSB", "MMFFTOR", "MMFFVDW",
})


def _read_sections() -> Dict[str, list]:
    """Read the Jmol parameter file and split into named sections.

    Returns a dict {section_name: [data_lines]} where section_name is the
    bare section identifier (e.g., 'MMFFBOND') and data_lines is the list
    of non-comment, non-blank lines within that section.

    Only sections listed in `_EXPECTED_SECTIONS` are retained; auxiliary
    sections such as MMFFBNDK (empirical-rule defaults) are skipped.
    """
    sections: Dict[str, list] = {name: [] for name in _EXPECTED_SECTIONS}
    current: Optional[str] = None
    with _PARAM_FILE.open("r", encoding="ascii", errors="replace") as f:
        for raw in f:
            m = _SECTION_HEADER.match(raw)
            if m:
                name = m.group(1).upper()
                current = name if name in _EXPECTED_SECTIONS else None
                continue
            if current is None:
                continue
            stripped = raw.strip()
            if not stripped or stripped.startswith(("*", "#")):
                continue
            sections[current].append(stripped)
    return sections


# ---------------------------------------------------------------------------
# Bond stretching: MMFFBOND
# ---------------------------------------------------------------------------

def _parse_bond_table(lines: list) -> Dict[Tuple[int, int, int], Tuple[float, float]]:
    """Parse the MMFFBOND section.

    Each data line: bond_type i j kb r0
    bond_type maps to bond_order roughly as:
      0 -> 1 (normal single), 1 -> 1 (aromatic-ish single), etc.
    We key the dict by (min(i,j), max(i,j), bond_type) so callers
    look up by the same bond_type the table uses (i.e., we keep
    Jmol's bond_type semantics).
    """
    table: Dict[Tuple[int, int, int], Tuple[float, float]] = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 5:
            continue
        try:
            bt = int(parts[0])
            i = int(parts[1])
            j = int(parts[2])
            kb = float(parts[3])
            r0 = float(parts[4])
        except ValueError:
            continue
        key = (min(i, j), max(i, j), bt)
        table[key] = (kb, r0)
    return table


def get_bond_params(type_i: int, type_j: int, bond_type: int = 0) -> Optional[Tuple[float, float]]:
    """Look up MMFF94 bond stretching parameters.

    Args:
        type_i, type_j: MMFF94 atom type numbers (1-99).
        bond_type: MMFF94 bond-type code (0=normal single, see Halgren 1996).

    Returns:
        (kb, r0) tuple or None if no entry.
    """
    table = _load()["bond"]
    key = (min(type_i, type_j), max(type_i, type_j), bond_type)
    return table.get(key)


# ---------------------------------------------------------------------------
# Angle bending: MMFFANG
# ---------------------------------------------------------------------------

def _parse_angle_table(lines: list) -> Dict[Tuple[int, int, int, int], Tuple[float, float]]:
    """Parse MMFFANG section. Format: angle_type i j k ka theta0.

    Keyed by (i_outer_min, j_center, k_outer_max, angle_type).
    """
    table: Dict[Tuple[int, int, int, int], Tuple[float, float]] = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 6:
            continue
        try:
            at = int(parts[0])
            i = int(parts[1])
            j = int(parts[2])
            k = int(parts[3])
            ka = float(parts[4])
            t0 = float(parts[5])
        except ValueError:
            continue
        i_min, k_max = (i, k) if i <= k else (k, i)
        table[(i_min, j, k_max, at)] = (ka, t0)
    return table


def get_angle_params(type_i: int, type_j: int, type_k: int,
                     angle_type: int = 0) -> Optional[Tuple[float, float]]:
    """Look up MMFF94 angle bending parameters.

    Args:
        type_i, type_j (center), type_k: MMFF94 atom type numbers.
        angle_type: MMFF94 angle-type code (0=normal).

    Returns: (ka, theta0_deg) or None.
    """
    table = _load()["angle"]
    i_min, k_max = (type_i, type_k) if type_i <= type_k else (type_k, type_i)
    return table.get((i_min, type_j, k_max, angle_type))


# ---------------------------------------------------------------------------
# Stretch-bend: MMFFSTBN
# ---------------------------------------------------------------------------

def _parse_sb_table(lines: list) -> Dict[Tuple[int, int, int, int], Tuple[float, float]]:
    """Parse MMFFSTBN. Format: sb_type i j k kbai kbak.

    NOTE: SB is NOT symmetric in i<->k because kbai applies to the i-j bond
    and kbak applies to the j-k bond. We key by (i, j, k, sb_type) literal.
    """
    table = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 6:
            continue
        try:
            sb = int(parts[0])
            i = int(parts[1]); j = int(parts[2]); k = int(parts[3])
            kbai = float(parts[4]); kbak = float(parts[5])
        except ValueError:
            continue
        table[(i, j, k, sb)] = (kbai, kbak)
    return table


def get_sb_params(type_i: int, type_j: int, type_k: int,
                  sb_type: int = 0) -> Optional[Tuple[float, float]]:
    """Look up MMFF94 stretch-bend parameters.

    Returns: (kbai, kbak) where kbai is for bond i-j and kbak is for j-k.
    Not symmetric in i<->k.
    """
    table = _load()["sb"]
    direct = table.get((type_i, type_j, type_k, sb_type))
    if direct is not None:
        return direct
    # Try swapped order with kbai<->kbak
    swapped = table.get((type_k, type_j, type_i, sb_type))
    if swapped is not None:
        return (swapped[1], swapped[0])
    return None


# ---------------------------------------------------------------------------
# Torsion: MMFFTOR
# ---------------------------------------------------------------------------

def _parse_torsion_table(lines: list) -> Dict[Tuple[int, int, int, int, int], Tuple[float, float, float]]:
    """Parse MMFFTOR. Format: tor_type i j k l V1 V2 V3."""
    table = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 8:
            continue
        try:
            tt = int(parts[0])
            i = int(parts[1]); j = int(parts[2]); k = int(parts[3]); l = int(parts[4])
            v1 = float(parts[5]); v2 = float(parts[6]); v3 = float(parts[7])
        except ValueError:
            continue
        # Canonical order: keep central bond j-k; sort by (i,l) tuple-wise
        if (j, k) > (k, j):
            j, k = k, j
            i, l = l, i
        table[(i, j, k, l, tt)] = (v1, v2, v3)
    return table


def get_torsion_params(type_i: int, type_j: int, type_k: int, type_l: int,
                       tor_type: int = 0) -> Optional[Tuple[float, float, float]]:
    """Look up MMFF94 torsion parameters. Returns (V1, V2, V3) or None."""
    table = _load()["torsion"]
    # Canonicalize the same way as the parser
    if (type_j, type_k) > (type_k, type_j):
        type_j, type_k = type_k, type_j
        type_i, type_l = type_l, type_i
    return table.get((type_i, type_j, type_k, type_l, tor_type))


# ---------------------------------------------------------------------------
# Out-of-plane: MMFFOOP
# ---------------------------------------------------------------------------

def _parse_oop_table(lines: list) -> Dict[Tuple[int, int, int, int], float]:
    """Parse MMFFOOP. Format: i j(center) k l koop.

    Keyed by (i_min, j_center, k_mid, l_max) -- canonicalize the outers.
    """
    table = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 5:
            continue
        try:
            i = int(parts[0]); j = int(parts[1]); k = int(parts[2]); l = int(parts[3])
            koop = float(parts[4])
        except ValueError:
            continue
        outers = tuple(sorted([i, k, l]))
        table[(outers[0], j, outers[1], outers[2])] = koop
    return table


def get_oop_params(type_i: int, type_j: int, type_k: int, type_l: int) -> Optional[float]:
    """Look up MMFF94 out-of-plane force constant koop. j is the center."""
    table = _load()["oop"]
    outers = tuple(sorted([type_i, type_k, type_l]))
    return table.get((outers[0], type_j, outers[1], outers[2]))


# ---------------------------------------------------------------------------
# Van der Waals: MMFFVDW
# ---------------------------------------------------------------------------

def _parse_vdw_table(lines: list) -> Dict[int, Tuple[float, float, float, float, int]]:
    """Parse MMFFVDW. Format: type alpha N A G DA_symbol.

    DA_symbol is encoded as int(ord('-')|'D'|'A'). The downstream
    VdwCalc uses int comparison (matches Jmol's `CalculationsMMFF.DA_D`
    and `DA_DA` integer constants).
    """
    table = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 6:
            continue
        try:
            t = int(parts[0])
            alpha = float(parts[1]); N = float(parts[2])
            A = float(parts[3]); G = float(parts[4])
            da_char = parts[5][0]
            da = ord(da_char) if da_char in ("D", "A") else ord("-")
        except (ValueError, IndexError):
            continue
        table[t] = (alpha, N, A, G, da)
    return table


def get_vdw_params(type_i: int) -> Optional[Tuple[float, float, float, float, int]]:
    """Look up MMFF94 VdW parameters. Returns (alpha, N, A, G, da_int) or None."""
    return _load()["vdw"].get(type_i)


# ---------------------------------------------------------------------------
# Bond-charge increment: MMFFCHG
# ---------------------------------------------------------------------------

def _parse_bci_table(lines: list) -> Dict[Tuple[int, int, int], float]:
    """Parse MMFFCHG. Format: bond_type i j bci. Stored asymmetric (i,j,bt)."""
    table = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            bt = int(parts[0]); i = int(parts[1]); j = int(parts[2])
            bci = float(parts[3])
        except ValueError:
            continue
        table[(i, j, bt)] = bci
    return table


def get_bci(type_i: int, type_j: int, bond_type: int = 0) -> Optional[float]:
    """Look up the MMFF94 BCI for type_i in an i-j bond.

    By convention bci(i,j) = -bci(j,i). The parsed table stores both
    directions if Jmol's file does; we try direct then negated swap.
    """
    table = _load()["bci"]
    direct = table.get((type_i, type_j, bond_type))
    if direct is not None:
        return direct
    swapped = table.get((type_j, type_i, bond_type))
    if swapped is not None:
        return -swapped
    return None


# ---------------------------------------------------------------------------
# Partial BCI: MMFFPBCI (fallback for missing MMFFCHG entries)
# ---------------------------------------------------------------------------

def _parse_pbci_table(lines: list) -> Dict[int, Tuple[float, float]]:
    """Parse MMFFPBCI. Format: '0 type pbci fcadj'.

    pbci is used as a fallback when MMFFCHG has no entry for an (i,j) pair:
    bci(i,j) ~ pbci(j) - pbci(i)
    fcadj is the formal-charge adjustment factor (used for charged systems).
    """
    table = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            t = int(parts[1]); pbci = float(parts[2]); fcadj = float(parts[3])
        except (ValueError, IndexError):
            continue
        table[t] = (pbci, fcadj)
    return table


def get_pbci(type_i: int) -> Optional[Tuple[float, float]]:
    """Look up MMFF94 partial BCI and formal-charge adjustment for an atom type.

    Returns (pbci, fcadj) or None.
    """
    return _load()["pbci"].get(type_i)


# ---------------------------------------------------------------------------
# Lazy loader with pickle cache
# ---------------------------------------------------------------------------

def _load() -> Dict[str, dict]:
    """Lazy load + cache the parsed parameter dicts.

    Cache behavior:
        - If cache file exists and is newer than the .txt source: deserialize.
        - Else: parse the .txt, write fresh cache, return.

    Thread-safe via _load_lock.
    """
    global _PARAMS
    if _PARAMS is not None:
        return _PARAMS
    with _load_lock:
        if _PARAMS is not None:
            return _PARAMS

        # Try to use the pickle cache
        if (_CACHE_FILE.exists()
                and _CACHE_FILE.stat().st_mtime >= _PARAM_FILE.stat().st_mtime):
            try:
                with _CACHE_FILE.open("rb") as f:
                    _PARAMS = pickle.load(f)
                return _PARAMS
            except (pickle.PickleError, EOFError, OSError):
                # Corrupted cache -- fall through to fresh parse
                pass

        # Fresh parse
        sections = _read_sections()
        _PARAMS = {
            "bond":    _parse_bond_table(sections.get("MMFFBOND", [])),
            "angle":   _parse_angle_table(sections.get("MMFFANG", [])),
            "sb":      _parse_sb_table(sections.get("MMFFSTBN", [])),
            "torsion": _parse_torsion_table(sections.get("MMFFTOR", [])),
            "oop":     _parse_oop_table(sections.get("MMFFOOP", [])),
            "vdw":     _parse_vdw_table(sections.get("MMFFVDW", [])),
            "bci":     _parse_bci_table(sections.get("MMFFCHG", [])),
            "pbci":    _parse_pbci_table(sections.get("MMFFPBCI", [])),
        }

        # Write fresh cache (best-effort; ignore write failures)
        try:
            _CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with _CACHE_FILE.open("wb") as f:
                pickle.dump(_PARAMS, f, protocol=pickle.HIGHEST_PROTOCOL)
        except OSError:
            pass

        return _PARAMS


# ---------------------------------------------------------------------------
# Atom-type properties: parsed from mmff94_atom_types.txt (separate file)
# ---------------------------------------------------------------------------

def _parse_atom_props_file() -> Dict[int, dict]:
    """Parse mmff94_atom_types.txt. Each non-comment line:
        AtSym  ElemNo  mmType  HType  formalCharge*12  valence  Desc  Smiles

    The MMFF94 'equivalence class' is the HType column in Jmol's file.
    When HType == 0 (no hydrogen-equivalence aggregation), we fall back to
    the type number itself.
    """
    table: Dict[int, dict] = {}
    with _ATOM_TYPES_FILE.open("r", encoding="ascii", errors="replace") as f:
        for raw in f:
            stripped = raw.strip()
            if not stripped or stripped.startswith(("#", "*")):
                continue
            parts = stripped.split(None, 7)  # split on whitespace, max 8 cols
            if len(parts) < 6:
                continue
            try:
                element_no = int(parts[1])
                mm_type = int(parts[2])
                h_type = int(parts[3])
                fc12 = int(parts[4])
                valence = int(parts[5])
            except ValueError:
                continue
            # Don't overwrite first occurrence (multiple SMARTS rows can share
            # the same mm_type -- the first entry wins per Jmol convention).
            if mm_type in table:
                continue
            table[mm_type] = {
                "symbol": parts[0],
                "element_no": element_no,
                "equiv_class": h_type if h_type != 0 else mm_type,
                "formal_charge_x12": fc12,
                "valence": valence,
                "description": parts[6] if len(parts) > 6 else "",
            }
    return table


def get_atom_props(mmff_type: int) -> Optional[dict]:
    """Look up element/class/formal-charge/valence for an MMFF94 atom type.

    Returns a dict with keys: symbol, element_no, equiv_class,
    formal_charge_x12, valence, description.  None if type unknown.
    """
    global _ATOM_PROPS
    if _ATOM_PROPS is None:
        with _atom_props_lock:
            if _ATOM_PROPS is None:
                _ATOM_PROPS = _parse_atom_props_file()
    return _ATOM_PROPS.get(mmff_type)
