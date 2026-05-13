"""
Pre-built NumPy arrays for the MMFF94 engine.

InteractionArrays is a frozen dataclass holding all 7 term's
parameter arrays. Built once per optimize_geometry call by ArraysBuilder
(implemented in Task 19) and consumed by every per-term calculator.

Design contract:
    - All arrays contiguous, dtype-locked.
    - Never resized between optimization steps.
    - Index arrays are int32 (matching numpy default for small ints).
    - Parameter arrays are float64.
"""
from __future__ import annotations
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class InteractionArrays:
    """All pre-cooked per-term data for one molecule, ready for vectorized ops."""

    n_atoms: int

    # ────────── Bond stretching ──────────
    bond_i: np.ndarray   # (Nb,) int32 — atom index i
    bond_j: np.ndarray   # (Nb,) int32 — atom index j
    bond_kb: np.ndarray  # (Nb,) float64 — force constant
    bond_r0: np.ndarray  # (Nb,) float64 — equilibrium length (Å)

    # ────────── Angle bending ──────────
    # theta0 in DEGREES (intentional — locks the units fix from spec §7.2)
    angle_i: np.ndarray; angle_j: np.ndarray; angle_k: np.ndarray   # int32
    angle_ka: np.ndarray; angle_theta0_deg: np.ndarray              # float64
    angle_is_linear: np.ndarray                                      # bool

    # ────────── Stretch-bend ──────────
    sb_i: np.ndarray; sb_j: np.ndarray; sb_k: np.ndarray             # int32
    sb_kbai: np.ndarray; sb_kbak: np.ndarray                         # float64
    sb_r0_ij: np.ndarray; sb_r0_jk: np.ndarray                       # float64
    sb_theta0_deg: np.ndarray                                        # float64

    # ────────── Torsion ──────────
    tor_i: np.ndarray; tor_j: np.ndarray
    tor_k: np.ndarray; tor_l: np.ndarray                             # int32
    tor_v1: np.ndarray; tor_v2: np.ndarray; tor_v3: np.ndarray       # float64

    # ────────── Out-of-plane bending ──────────
    oop_center: np.ndarray                                           # int32
    oop_i: np.ndarray; oop_j: np.ndarray; oop_k: np.ndarray          # int32
    oop_koop: np.ndarray                                             # float64

    # ────────── Van der Waals ──────────
    vdw_i: np.ndarray; vdw_j: np.ndarray                             # int32
    vdw_rs: np.ndarray; vdw_eps: np.ndarray                          # float64
                                                                     # combined R*, ε

    # ────────── Electrostatic ──────────
    es_i: np.ndarray; es_j: np.ndarray                               # int32
    es_qq: np.ndarray                                                # float64 — qi*qj precomputed
    es_factor: np.ndarray                                            # float64 — 249.0537 (1-4) or 332.0716


# ──────────────────────────── ArraysBuilder ────────────────────────────
import math
from src.core.domain.models.molecule import Molecule
from src.services.forcefield.parameters import (
    get_bond_params, get_angle_params, get_sb_params,
    get_torsion_params, get_oop_params, get_vdw_params,
)
from src.services.forcefield._engine.exclusions import build_exclusions

# MMFF94 OOP-eligible classes (Jmol's CalculationsMMFF.isInvertible whitelist).
_OOP_ELIGIBLE_CLASSES = frozenset({
    2, 3, 10, 30, 37, 39, 40, 41, 45, 49,
    54, 55, 56, 57, 58, 63, 64, 67, 69, 78, 80, 81,
})

# ES factors per MMFF94: 1-4 pairs scaled to 0.75 x 332.0716.
_ES_FACTOR_15PLUS = 332.0716
_ES_FACTOR_14 = 249.0537


class ArraysBuilder:
    """Build all 7 term arrays from a Molecule. Stateless — call build()."""

    @staticmethod
    def build(mol: Molecule) -> InteractionArrays:
        """Build complete InteractionArrays.

        Pre: mol.atoms have mmff_type and mmff_class set (run AtomTyper first).
        Pre: mol.atoms have partial_charge set (run assign_bci_charges first).
        """
        n_atoms = len(mol.atoms)
        excl = build_exclusions(mol)

        # ── Bonds ──
        bond_data = []
        for b in mol.bonds:
            ti = mol.atoms[b.begin_atom_idx].mmff_type
            tj = mol.atoms[b.end_atom_idx].mmff_type
            if ti == 0 or tj == 0:
                # Use class-level fallback parameters if type is missing
                params = (4.0, 1.50)  # default
            else:
                params = get_bond_params(ti, tj, 0)
                if params is None:
                    # Try class fallback
                    ci = mol.atoms[b.begin_atom_idx].mmff_class
                    cj = mol.atoms[b.end_atom_idx].mmff_class
                    params = get_bond_params(ci, cj, 0)
            kb, r0 = params if params is not None else (4.0, 1.50)
            bond_data.append((b.begin_atom_idx, b.end_atom_idx, kb, r0))

        if bond_data:
            bi, bj, bkb, br0 = zip(*bond_data)
        else:
            bi = bj = (); bkb = br0 = ()

        bond_i = np.asarray(bi, dtype=np.int32)
        bond_j = np.asarray(bj, dtype=np.int32)
        bond_kb = np.asarray(bkb, dtype=np.float64)
        bond_r0 = np.asarray(br0, dtype=np.float64)

        # ── Angles ──
        angle_data = []
        sb_data = []
        for atom in mol.atoms:
            neighbors = mol.get_neighbors(atom.index)
            if len(neighbors) < 2:
                continue
            for a in range(len(neighbors)):
                for c in range(a + 1, len(neighbors)):
                    i_idx = neighbors[a]
                    k_idx = neighbors[c]
                    ti = mol.atoms[i_idx].mmff_type
                    tj = atom.mmff_type
                    tk = mol.atoms[k_idx].mmff_type
                    if any(t == 0 for t in (ti, tj, tk)):
                        ka, t0 = 0.5, 109.5
                    else:
                        ap = get_angle_params(ti, tj, tk, 0)
                        if ap is None:
                            # class fallback
                            ci = mol.atoms[i_idx].mmff_class
                            cj = atom.mmff_class
                            ck = mol.atoms[k_idx].mmff_class
                            ap = get_angle_params(ci, cj, ck, 0)
                        ka, t0 = ap if ap is not None else (0.5, 109.5)
                    is_linear = abs(t0 - 180.0) < 1.0
                    angle_data.append((i_idx, atom.index, k_idx, ka, t0, is_linear))

                    # Stretch-bend (skip for linear angles).
                    if not is_linear and not (ti == 0 or tj == 0 or tk == 0):
                        sbp = get_sb_params(ti, tj, tk, 0)
                        if sbp is not None:
                            r0_ij_params = get_bond_params(ti, tj, 0) or (4.0, 1.50)
                            r0_jk_params = get_bond_params(tj, tk, 0) or (4.0, 1.50)
                            sb_data.append((
                                i_idx, atom.index, k_idx,
                                sbp[0], sbp[1],
                                r0_ij_params[1], r0_jk_params[1],
                                t0,
                            ))

        if angle_data:
            ai_l, aj_l, ak_l, aka_l, at0_l, alin_l = zip(*angle_data)
        else:
            ai_l = aj_l = ak_l = (); aka_l = at0_l = (); alin_l = ()

        angle_i = np.asarray(ai_l, dtype=np.int32)
        angle_j = np.asarray(aj_l, dtype=np.int32)
        angle_k = np.asarray(ak_l, dtype=np.int32)
        angle_ka = np.asarray(aka_l, dtype=np.float64)
        angle_theta0_deg = np.asarray(at0_l, dtype=np.float64)
        angle_is_linear = np.asarray(alin_l, dtype=bool)

        # ── SB ──
        if sb_data:
            si, sj, sk, skbai, skbak, sr0ij, sr0jk, st0 = zip(*sb_data)
        else:
            si = sj = sk = (); skbai = skbak = sr0ij = sr0jk = st0 = ()
        sb_i = np.asarray(si, dtype=np.int32)
        sb_j = np.asarray(sj, dtype=np.int32)
        sb_k = np.asarray(sk, dtype=np.int32)
        sb_kbai = np.asarray(skbai, dtype=np.float64)
        sb_kbak = np.asarray(skbak, dtype=np.float64)
        sb_r0_ij = np.asarray(sr0ij, dtype=np.float64)
        sb_r0_jk = np.asarray(sr0jk, dtype=np.float64)
        sb_theta0_deg = np.asarray(st0, dtype=np.float64)

        # ── Torsions ──
        tor_data = []
        for b in mol.bonds:
            j, k = b.begin_atom_idx, b.end_atom_idx
            for i_idx in mol.get_neighbors(j):
                if i_idx == k:
                    continue
                for l_idx in mol.get_neighbors(k):
                    if l_idx == j or l_idx == i_idx:
                        continue
                    ti = mol.atoms[i_idx].mmff_type
                    tj = mol.atoms[j].mmff_type
                    tk = mol.atoms[k].mmff_type
                    tl = mol.atoms[l_idx].mmff_type
                    if any(t == 0 for t in (ti, tj, tk, tl)):
                        v1, v2, v3 = 0.0, 0.0, 0.0
                    else:
                        tp = get_torsion_params(ti, tj, tk, tl, 0)
                        if tp is None:
                            ci = mol.atoms[i_idx].mmff_class
                            cj = mol.atoms[j].mmff_class
                            ck = mol.atoms[k].mmff_class
                            cl = mol.atoms[l_idx].mmff_class
                            tp = get_torsion_params(ci, cj, ck, cl, 0)
                        v1, v2, v3 = tp if tp is not None else (0.0, 0.0, 0.0)
                    tor_data.append((i_idx, j, k, l_idx, v1, v2, v3))

        if tor_data:
            ti_l, tj_l, tk_l, tl_l, tv1, tv2, tv3 = zip(*tor_data)
        else:
            ti_l = tj_l = tk_l = tl_l = (); tv1 = tv2 = tv3 = ()
        tor_i = np.asarray(ti_l, dtype=np.int32)
        tor_j = np.asarray(tj_l, dtype=np.int32)
        tor_k = np.asarray(tk_l, dtype=np.int32)
        tor_l = np.asarray(tl_l, dtype=np.int32)
        tor_v1 = np.asarray(tv1, dtype=np.float64)
        tor_v2 = np.asarray(tv2, dtype=np.float64)
        tor_v3 = np.asarray(tv3, dtype=np.float64)

        # ── OOP ──
        oop_data = []
        for atom in mol.atoms:
            if len(mol.get_neighbors(atom.index)) != 3:
                continue
            if atom.mmff_class not in _OOP_ELIGIBLE_CLASSES:
                continue
            neighbors = list(mol.get_neighbors(atom.index))
            # Three Wilson-angle entries per center (one per outer atom).
            triples = [
                (neighbors[0], atom.index, neighbors[1], neighbors[2]),
                (neighbors[1], atom.index, neighbors[2], neighbors[0]),
                (neighbors[2], atom.index, neighbors[0], neighbors[1]),
            ]
            tc = atom.mmff_type
            ts = [mol.atoms[n].mmff_type for n in neighbors]
            koop = get_oop_params(tc, ts[0], ts[1], ts[2]) or 0.0
            for tup in triples:
                oop_data.append((tup[1], tup[0], tup[2], tup[3], koop))

        if oop_data:
            oc, oi, oj, ok, okk = zip(*oop_data)
        else:
            oc = oi = oj = ok = (); okk = ()
        oop_center = np.asarray(oc, dtype=np.int32)
        oop_i = np.asarray(oi, dtype=np.int32)
        oop_j = np.asarray(oj, dtype=np.int32)
        oop_k = np.asarray(ok, dtype=np.int32)
        oop_koop = np.asarray(okk, dtype=np.float64)

        # ── VdW + ES — small-N path or cell-list ──
        N_CELL_THRESHOLD = 200
        if n_atoms < N_CELL_THRESHOLD:
            candidate_pairs = [
                (i, j)
                for i in range(n_atoms)
                for j in range(i + 1, n_atoms)
            ]
        else:
            from src.services.forcefield._engine.neighbor_list import (
                build_neighbor_pairs, VDW_CUTOFF, ES_CUTOFF,
            )
            coords = np.array([[a.x, a.y, a.z] for a in mol.atoms],
                              dtype=np.float64)
            # Use the larger cutoff (ES) for the candidate list; filter
            # by VdW cutoff within the inner branch.
            candidate_arr = build_neighbor_pairs(coords, ES_CUTOFF)
            candidate_pairs = [tuple(p) for p in candidate_arr]

        vdw_data = []
        es_data = []
        for i, j in candidate_pairs:
            pair = frozenset({i, j})
            if pair in excl.pairs_12 or pair in excl.pairs_13:
                continue
            ti = mol.atoms[i].mmff_type
            tj = mol.atoms[j].mmff_type
            if ti != 0 and tj != 0:
                vp_i = get_vdw_params(ti)
                vp_j = get_vdw_params(tj)
                if vp_i is not None and vp_j is not None:
                    rs, eps = _combine_vdw(vp_i, vp_j)
                    if eps > 0:
                        vdw_data.append((i, j, rs, eps))
            qi = mol.atoms[i].partial_charge
            qj = mol.atoms[j].partial_charge
            if abs(qi) > 1e-6 and abs(qj) > 1e-6:
                f = _ES_FACTOR_14 if pair in excl.pairs_14 else _ES_FACTOR_15PLUS
                es_data.append((i, j, qi * qj, f))

        if vdw_data:
            vi, vj, vrs, veps = zip(*vdw_data)
        else:
            vi = vj = (); vrs = veps = ()
        vdw_i = np.asarray(vi, dtype=np.int32)
        vdw_j = np.asarray(vj, dtype=np.int32)
        vdw_rs = np.asarray(vrs, dtype=np.float64)
        vdw_eps = np.asarray(veps, dtype=np.float64)

        if es_data:
            ei, ej, eqq, ef = zip(*es_data)
        else:
            ei = ej = (); eqq = ef = ()
        es_i = np.asarray(ei, dtype=np.int32)
        es_j = np.asarray(ej, dtype=np.int32)
        es_qq = np.asarray(eqq, dtype=np.float64)
        es_factor = np.asarray(ef, dtype=np.float64)

        return InteractionArrays(
            n_atoms=n_atoms,
            bond_i=bond_i, bond_j=bond_j, bond_kb=bond_kb, bond_r0=bond_r0,
            angle_i=angle_i, angle_j=angle_j, angle_k=angle_k,
            angle_ka=angle_ka, angle_theta0_deg=angle_theta0_deg,
            angle_is_linear=angle_is_linear,
            sb_i=sb_i, sb_j=sb_j, sb_k=sb_k,
            sb_kbai=sb_kbai, sb_kbak=sb_kbak,
            sb_r0_ij=sb_r0_ij, sb_r0_jk=sb_r0_jk, sb_theta0_deg=sb_theta0_deg,
            tor_i=tor_i, tor_j=tor_j, tor_k=tor_k, tor_l=tor_l,
            tor_v1=tor_v1, tor_v2=tor_v2, tor_v3=tor_v3,
            oop_center=oop_center, oop_i=oop_i, oop_j=oop_j, oop_k=oop_k,
            oop_koop=oop_koop,
            vdw_i=vdw_i, vdw_j=vdw_j, vdw_rs=vdw_rs, vdw_eps=vdw_eps,
            es_i=es_i, es_j=es_j, es_qq=es_qq, es_factor=es_factor,
        )


# Donor/Acceptor codes (matches CalculationsMMFF.DA_D / DA_A in Jmol).
_DA_D = ord("D")
_DA_A = ord("A")
_DA_DA = _DA_D + _DA_A


def _combine_vdw(vp_i, vp_j):
    """MMFF94 combined VdW radius and well depth.

    vp_x = (alpha, N, A, G, da_flag) where da_flag is ord('-'/'D'/'A').

    Returns (rs_combined, eps_combined).
    """
    alpha_i, N_i, A_i, G_i, da_i = vp_i
    alpha_j, N_j, A_j, G_j, da_j = vp_j

    rs_ii = A_i * (alpha_i ** 0.25)
    rs_jj = A_j * (alpha_j ** 0.25)
    if rs_ii + rs_jj == 0.0:
        return (0.0, 0.0)

    gamma = (rs_ii - rs_jj) / (rs_ii + rs_jj)
    rs = 0.5 * (rs_ii + rs_jj)

    # Donor exception: no rs inflation if either atom is a donor.
    if da_i != _DA_D and da_j != _DA_D:
        rs *= (1.0 + 0.2 * (1.0 - math.exp(-12.0 * gamma * gamma)))

    eps_denom = math.sqrt(alpha_i / N_i) + math.sqrt(alpha_j / N_j)
    if eps_denom == 0.0:
        return (rs, 0.0)
    eps = ((181.16 * G_i * G_j * alpha_i * alpha_j) / eps_denom) * (rs ** -6.0)

    # Donor-acceptor exception
    if da_i + da_j == _DA_DA:
        rs *= 0.8
        eps *= 0.5

    return (rs, eps)
