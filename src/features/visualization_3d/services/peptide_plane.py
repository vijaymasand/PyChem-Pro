"""
Peptide Plane — Faithful port of CPPCartoon's PeptidePlane.hpp

Each PeptidePlane spans three consecutive residues (r1, r2, r3) and defines
a local coordinate frame at the midpoint between CA1 and CA2:

    Forward = normalize(CA2 - CA1)
    b       = normalize(O1  - CA1)
    Normal  = normalize(cross(Forward, b))
    Side    = normalize(cross(Normal, Forward))

The Transition() function determines the effective SS type at segment
boundaries, and Flip() corrects orientation discontinuities.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Tuple

# SS type constants matching CPPCartoon
HELIX  = 'H'
STRAND = 'E'
COIL   = 'C'


@dataclass
class PeptidePlane:
    """Mathematical representation of a peptide bond plane (CPPCartoon-faithful)."""
    position: np.ndarray         # Midpoint of CA1-CA2
    normal: np.ndarray           # Plane normal (cross product direction)
    forward: np.ndarray          # Backbone tangent CA1 -> CA2
    side: np.ndarray             # Side vector (in-plane, perpendicular to forward)
    flipped: bool = False
    # SS types for the 3 residues that define this plane
    ss1: str = COIL              # SS of residue 1
    ss2: str = COIL              # SS of residue 2
    ss3: str = COIL              # SS of residue 3
    # Residue IDs for discontinuity detection
    res_id1: int = 0
    res_id2: int = 0
    res_id3: int = 0
    chain_id: str = ''


def _safe_normalize(v):
    """Normalize a vector, returning zero vector if length is near zero."""
    n = np.linalg.norm(v)
    if n < 1e-10:
        return np.zeros(3)
    return v / n


def new_peptide_plane(ca1_pos, o1_pos, ca2_pos,
                      ss1=COIL, ss2=COIL, ss3=COIL,
                      res_id1=0, res_id2=0, res_id3=0,
                      chain_id=''):
    """
    Construct a PeptidePlane from atom coordinates.
    Matches CPPCartoon's NewPeptidePlane exactly.
    
    Args:
        ca1_pos: CA position of residue 1 (np.ndarray shape (3,))
        o1_pos:  O  position of residue 1 (np.ndarray shape (3,))
        ca2_pos: CA position of residue 2 (np.ndarray shape (3,))
        ss1/ss2/ss3: secondary structure types for residues 1, 2, 3
    
    Returns:
        PeptidePlane or None if atoms are degenerate
    """
    # Forward: direction from CA1 to CA2
    a = _safe_normalize(ca2_pos - ca1_pos)
    if np.linalg.norm(a) < 1e-10:
        return None

    # b: direction from CA1 to O1
    b = _safe_normalize(o1_pos - ca1_pos)
    if np.linalg.norm(b) < 1e-10:
        return None

    # Normal = cross(forward, b) — perpendicular to the peptide plane
    c = _safe_normalize(np.cross(a, b))
    if np.linalg.norm(c) < 1e-10:
        return None

    # Side = cross(Normal, forward) — in the peptide plane, perpendicular to backbone
    d = _safe_normalize(np.cross(c, a))

    # Position = midpoint of CA1 and CA2
    p = (ca1_pos + ca2_pos) / 2.0

    return PeptidePlane(
        position=p,
        normal=c,
        forward=a,
        side=d,
        flipped=False,
        ss1=ss1, ss2=ss2, ss3=ss3,
        res_id1=res_id1, res_id2=res_id2, res_id3=res_id3,
        chain_id=chain_id,
    )


def transition(pp):
    """
    Determine the effective SS types at the start and end of this plane's segment.
    Matches CPPCartoon's Transition() function.
    
    Returns: (type1, type2)
    """
    t1 = pp.ss1
    t2 = pp.ss2
    t3 = pp.ss3

    type1 = t2
    type2 = t2

    # If SS changed going INTO this residue, use the previous type for type1
    if t2 > t1 and t2 == t3:
        type1 = t1
    # If SS changes going OUT of this residue, use the next type for type2
    if t2 > t3 and t1 == t2:
        type2 = t3

    return type1, type2


def flip(pp):
    """Flip the Side and Normal vectors. Matches CPPCartoon's Flip()."""
    pp.side = -pp.side
    pp.normal = -pp.normal
    pp.flipped = not pp.flipped


def compute_chain_planes(molecule):
    """
    Extract PeptidePlanes for each chain in the molecule.
    
    Matches CPPCartoon's createChainMesh plane-building loop:
    - Each plane uses 3 consecutive residues (r1, r2, r3)
    - Boundary residues are clamped (doubled)
    - Start plane position is pinned to first CA
    - End plane position is pinned to last CA
    - Side vectors are flipped for consistency
    
    Returns: dict of {chain_id: list[PeptidePlane]}
    """
    # 1. Group atoms by chain -> residue -> atom name
    chains_data = {}
    for atom in molecule.atoms:
        chain = getattr(atom, 'chain_id', '') or 'A'
        res_seq = getattr(atom, 'res_seq', 0)
        pdb_name_raw = getattr(atom, 'pdb_name', '')
        pdb_name = pdb_name_raw.strip() if pdb_name_raw else ''
        if not pdb_name:
            continue

        if chain not in chains_data:
            chains_data[chain] = {}
        if res_seq not in chains_data[chain]:
            chains_data[chain][res_seq] = {'atoms': {}, 'ss': COIL}
        chains_data[chain][res_seq]['atoms'][pdb_name] = atom

        # Record SS type
        ss = getattr(atom, 'ss_type', COIL) or COIL
        if pdb_name == 'CA':
            chains_data[chain][res_seq]['ss'] = ss

    result = {}

    for chain_id in sorted(chains_data.keys()):
        residues = chains_data[chain_id]
        sorted_seqs = sorted(residues.keys())
        n = len(sorted_seqs)

        if n < 2:
            continue

        # Build planes using 3-residue windows, matching CPPCartoon's loop
        # i goes from -1 to n-1 (inclusive), covering n+1 planes
        planes = []
        for i in range(-1, n):
            # Clamp indices at boundaries (exactly like CPPCartoon)
            if i == -1:
                id0, id1, id2 = 0, 0, 0 + 1 if n > 1 else 0
            elif i >= n - 2:
                if i == n - 2:
                    id0, id1, id2 = i, i + 1, i + 1
                else:  # i == n - 1
                    id0, id1, id2 = i - 1, i, i
            else:
                id0, id1, id2 = i, i + 1, i + 2

            # Clamp to valid range
            id0 = max(0, min(id0, n - 1))
            id1 = max(0, min(id1, n - 1))
            id2 = max(0, min(id2, n - 1))

            seq0 = sorted_seqs[id0]
            seq1 = sorted_seqs[id1]
            seq2 = sorted_seqs[id2]

            r1_atoms = residues[seq0]['atoms']
            r2_atoms = residues[seq1]['atoms']

            ca1 = r1_atoms.get('CA')
            o1 = r1_atoms.get('O')
            ca2 = r2_atoms.get('CA')

            if ca1 is None or o1 is None or ca2 is None:
                continue
            if not all(hasattr(a, 'x') for a in [ca1, o1, ca2]):
                continue

            ca1_pos = np.array([ca1.x, ca1.y, ca1.z], dtype=np.float64)
            o1_pos = np.array([o1.x, o1.y, o1.z], dtype=np.float64)
            ca2_pos = np.array([ca2.x, ca2.y, ca2.z], dtype=np.float64)

            ss1 = residues[seq0]['ss']
            ss2 = residues[seq1]['ss']
            ss3 = residues[seq2]['ss']

            pp = new_peptide_plane(
                ca1_pos, o1_pos, ca2_pos,
                ss1=ss1, ss2=ss2, ss3=ss3,
                res_id1=seq0, res_id2=seq1, res_id3=seq2,
                chain_id=chain_id,
            )
            if pp is None:
                continue

            # Pin start position to first CA
            if i <= 0:
                first_ca = residues[sorted_seqs[0]]['atoms'].get('CA')
                if first_ca:
                    pp.position = np.array([first_ca.x, first_ca.y, first_ca.z], dtype=np.float64)

            # Pin end position to last CA
            if i >= n - 2:
                last_ca = residues[sorted_seqs[-1]]['atoms'].get('CA')
                if last_ca:
                    pp.position = np.array([last_ca.x, last_ca.y, last_ca.z], dtype=np.float64)

            planes.append(pp)

        # Flip side vectors for consistency (matches CPPCartoon lines 513-519)
        if len(planes) > 1:
            for i in range(1, len(planes)):
                if np.dot(planes[i].side, planes[i - 1].side) < 0:
                    flip(planes[i])

        if len(planes) >= 4:
            result[chain_id] = planes

    return result
