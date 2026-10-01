"""
Atom correspondence — decide *which* atoms of two molecules pair up before a
rigid-body fit.

Three strategies are provided, matching the three practical alignment cases:

* :func:`pairs_by_index`     — conformers / docking poses of the *same* molecule
  (atom order preserved): pair filtered atom *k* of A with atom *k* of B.
* :func:`pairs_by_element`   — same as index but only keep pairs whose element
  symbols agree (a safety filter for slightly different atom orders).
* :func:`pairs_by_sequence`  — two protein chains: align their one-letter Cα
  sequences (Needleman-Wunsch) and pair the matched residues' Cα atoms.

An *atom selection* narrows which atoms participate: ``all`` / ``heavy``
(non-hydrogen) / ``ca`` (protein α-carbons) / ``backbone`` (protein N, CA, C, O).
"""
from __future__ import annotations

from src.services.alignment._engine.sequence_align import global_align

# Standard + a few non-standard residues → one-letter codes.
_THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'SEC': 'U', 'PYL': 'O',
}

_BACKBONE_NAMES = ('N', 'CA', 'C', 'O')


def is_protein(mol) -> bool:
    """True if the molecule looks like a protein.

    Prefers the ``is_protein`` flag set by the loader, but falls back to
    counting recognizable α-carbons — the parallel PDB path and some readers
    don't set the flag, so relying on it alone would misclassify large
    structures.
    """
    if mol.properties.get('is_protein'):
        return True
    count = 0
    for atom in mol.atoms:
        if _is_ca(atom):
            count += 1
            if count >= 3:
                return True
    return False


def has_coords(atom) -> bool:
    return atom.x is not None


def selection_indices(mol, selection='heavy'):
    """Return the atom indices participating in a fit for *selection*.

    Args:
        mol:       Molecule to select from.
        selection: ``'all'`` | ``'heavy'`` | ``'ca'`` | ``'backbone'``.

    Returns:
        List of atom indices (only atoms carrying 3-D coordinates).
    """
    sel = (selection or 'heavy').lower()
    out = []
    for atom in mol.atoms:
        if not has_coords(atom):
            continue
        if sel == 'all':
            out.append(atom.index)
        elif sel == 'heavy':
            if atom.symbol != 'H':
                out.append(atom.index)
        elif sel == 'ca':
            if _is_ca(atom):
                out.append(atom.index)
        elif sel == 'backbone':
            if _is_backbone(atom):
                out.append(atom.index)
        else:
            raise ValueError(f"unknown selection: {selection!r}")
    return out


def _is_ca(atom) -> bool:
    name = (getattr(atom, 'pdb_name', None) or '').strip()
    return name == 'CA' and (getattr(atom, 'res_name', None) in _THREE_TO_ONE)


def _is_backbone(atom) -> bool:
    name = (getattr(atom, 'pdb_name', None) or '').strip()
    return name in _BACKBONE_NAMES and (getattr(atom, 'res_name', None) in _THREE_TO_ONE)


def pairs_by_index(mol_a, mol_b, selection='heavy'):
    """Pair the k-th selected atom of A with the k-th selected atom of B."""
    ia = selection_indices(mol_a, selection)
    ib = selection_indices(mol_b, selection)
    n = min(len(ia), len(ib))
    return list(zip(ia[:n], ib[:n]))


def pairs_by_element(mol_a, mol_b, selection='heavy'):
    """Index pairing, but keep only pairs whose element symbols agree."""
    ia = selection_indices(mol_a, selection)
    ib = selection_indices(mol_b, selection)
    pairs = []
    for a_idx, b_idx in zip(ia, ib):
        if mol_a.atoms[a_idx].symbol == mol_b.atoms[b_idx].symbol:
            pairs.append((a_idx, b_idx))
    return pairs


def extract_ca_sequence(mol):
    """Return ``(ca_indices, one_letter_sequence)`` for a protein.

    One Cα per residue, in chain/sequence order.  ``ca_indices[k]`` is the atom
    index of the α-carbon whose residue contributes ``sequence[k]``.
    """
    ca_indices = []
    sequence = []
    seen = set()
    for atom in mol.atoms:
        if not has_coords(atom) or not _is_ca(atom):
            continue
        key = (getattr(atom, 'chain_id', None), getattr(atom, 'res_seq', None))
        if key in seen:
            continue  # only the first Cα of a residue (guards altlocs / duplicates)
        seen.add(key)
        ca_indices.append(atom.index)
        sequence.append(_THREE_TO_ONE.get(atom.res_name, 'X'))
    return ca_indices, ''.join(sequence)


def pairs_by_sequence(mol_a, mol_b, **scoring):
    """Pair Cα atoms of two proteins via Needleman-Wunsch sequence alignment.

    Extra keyword arguments (``match``, ``mismatch``, ``gap_open``,
    ``gap_extend``, ``substitution``) are forwarded to
    :func:`~src.services.alignment._engine.sequence_align.global_align`.
    """
    ca_a, seq_a = extract_ca_sequence(mol_a)
    ca_b, seq_b = extract_ca_sequence(mol_b)
    if not seq_a or not seq_b:
        return []
    alignment = global_align(seq_a, seq_b, **scoring)
    return [(ca_a[si], ca_b[ti]) for si, ti in alignment.matched_pairs]
