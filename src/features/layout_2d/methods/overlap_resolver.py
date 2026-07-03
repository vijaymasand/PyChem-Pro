"""
Deterministic 2D overlap resolver.

OASA's coordinate generator occasionally produces a depiction where an acyclic
substituent is drawn on top of another part of the molecule (a classic case is
two ortho substituents on a ring whose atoms collide, e.g. aspirin).  This pass
repairs such cases *deterministically* and conservatively:

  For every acyclic single bond (a "rotatable" bond whose removal splits the
  molecule in two), reflect the smaller side across the bond axis — a
  chemically-equivalent 2-D flip that never changes connectivity or bond
  lengths — and keep the flip **only if it strictly reduces** the number of
  heavy-atom clashes plus heavy–heavy bond crossings.

Because it only ever accepts a strict improvement and processes bonds in a fixed
order, the result is stable and can never make a layout worse.  Rings are never
touched.  Operates in place on ``atom.x2d`` / ``atom.y2d``.
"""
from __future__ import annotations

from collections import deque

_CLASH_DIST = 0.5      # heavy atoms closer than this (bond length ~1.0) clash
_CLASH_DIST_SQ = _CLASH_DIST * _CLASH_DIST
_MAX_HEAVY = 80        # skip very large molecules (cost is ~O(bonds * heavy^2))


def resolve_overlaps_2d(mol) -> bool:
    """Reduce 2-D overlaps in *mol* by flipping acyclic branches.

    Returns True if any coordinates were changed.
    """
    heavy = [a for a in mol.atoms if a.symbol != 'H' and a.x2d is not None]
    if len(heavy) < 4 or len(heavy) > _MAX_HEAVY:
        return False

    coords = {a.index: [float(a.x2d), float(a.y2d)] for a in mol.atoms if a.x2d is not None}
    heavy_idx = [a.index for a in heavy]
    heavy_bonds = [(b.begin_atom_idx, b.end_atom_idx) for b in mol.bonds
                   if b.begin_atom_idx in coords and b.end_atom_idx in coords
                   and mol.atoms[b.begin_atom_idx].symbol != 'H'
                   and mol.atoms[b.end_atom_idx].symbol != 'H']

    base = _count_problems(coords, heavy_idx, heavy_bonds)
    if base == 0:
        return False

    changed = False
    # Deterministic order: by bond index.
    for bond in mol.bonds:
        if base == 0:
            break
        if bond.bond_type != 1 or bond.is_in_ring:   # single, non-ring only
            continue
        u, v = bond.begin_atom_idx, bond.end_atom_idx
        if u not in coords or v not in coords:
            continue
        branch = _branch_atoms(mol, keep=u, start=v, coords=coords)
        if branch is None:
            continue  # bond is inside a ring / not a clean split
        # Flip the smaller side for minimal disruption.
        if len(branch) * 2 > len(coords):
            branch = _branch_atoms(mol, keep=v, start=u, coords=coords)
            if branch is None:
                continue
            ax, ay = coords[v]
            bx, by = coords[u]
        else:
            ax, ay = coords[u]
            bx, by = coords[v]

        trial = {i: xy[:] for i, xy in coords.items()}
        for i in branch:
            trial[i] = _reflect(coords[i], ax, ay, bx, by)
        score = _count_problems(trial, heavy_idx, heavy_bonds)
        if score < base:
            coords = trial
            base = score
            changed = True

    if changed:
        for a in mol.atoms:
            if a.index in coords:
                a.x2d, a.y2d = coords[a.index][0], coords[a.index][1]
    return changed


def _branch_atoms(mol, keep, start, coords):
    """Atoms on *start*'s side when the keep–start bond is cut.

    Returns None if the traversal reaches *keep* by another route (the bond is
    part of a ring, so cutting it does not isolate a branch).
    """
    seen = {start}
    queue = deque([start])
    while queue:
        cur = queue.popleft()
        for nb in mol.get_neighbors(cur):
            if nb == keep and cur == start:
                continue  # the cut bond itself
            if nb == keep:
                return None  # cycle back to keep → ring bond, not a branch
            if nb not in seen and nb in coords:
                seen.add(nb)
                queue.append(nb)
    return seen


def _reflect(p, ax, ay, bx, by):
    """Reflect point *p* across the line through (ax,ay)-(bx,by)."""
    dx, dy = bx - ax, by - ay
    dd = dx * dx + dy * dy
    if dd < 1e-12:
        return p[:]
    t = ((p[0] - ax) * dx + (p[1] - ay) * dy) / dd
    projx, projy = ax + t * dx, ay + t * dy
    return [2.0 * projx - p[0], 2.0 * projy - p[1]]


def _count_problems(coords, heavy_idx, heavy_bonds):
    """Heavy-atom clashes + heavy–heavy bond crossings."""
    n = 0
    # Atom clashes (non-bonded pairs closer than the clash distance).
    for i in range(len(heavy_idx)):
        xi, yi = coords[heavy_idx[i]]
        for j in range(i + 1, len(heavy_idx)):
            xj, yj = coords[heavy_idx[j]]
            dx, dy = xi - xj, yi - yj
            if dx * dx + dy * dy < _CLASH_DIST_SQ:
                n += 1
    # Bond crossings (segments not sharing an endpoint).
    m = len(heavy_bonds)
    for i in range(m):
        a, b = heavy_bonds[i]
        pa, pb = coords[a], coords[b]
        for j in range(i + 1, m):
            c, d = heavy_bonds[j]
            if a == c or a == d or b == c or b == d:
                continue
            if _segments_cross(pa, pb, coords[c], coords[d]):
                n += 1
    return n


def _segments_cross(p1, p2, p3, p4):
    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])
    return (ccw(p1, p3, p4) != ccw(p2, p3, p4)) and (ccw(p1, p2, p3) != ccw(p1, p2, p4))
