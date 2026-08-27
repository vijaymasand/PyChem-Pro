"""
Hand-coded MMFF94 atom typer (Phase 1 implementation).

Covers ~43 numerical MMFF94 atom types via per-element decision-tree rules.
Designed to be replaced wholesale in Phase 3 by SmartsAtomTyper behind
the same AtomTyper interface — engine and parameter code unaffected.

Two-pass typing:
    1. Heavy atoms first (rules read only neighbors' element/hybridization,
       not their type — so order within heavies doesn't matter).
    2. Hydrogens second (rules read the parent heavy atom's mmff_type).
"""
from __future__ import annotations
from typing import Tuple, Callable, Dict

from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom
from src.services.forcefield.atom_typing.types import class_for_type
from src.services.forcefield.atom_typing._ring_info import build_ring_info, RingInfo


class HandCodedAtomTyper:
    """Phase 1 MMFF94 atom typer. ~43-type coverage of organic chemistry."""

    def type_atoms(self, mol: Molecule) -> None:
        """Populate atom.mmff_type and atom.mmff_class on every atom in mol."""
        # Setup: ring perception, hybridization, aromaticity.
        mol.assign_hybridization()
        mol.perceive_aromaticity()
        ring_info = build_ring_info(mol)

        # Pass 1: heavy atoms.
        for atom in mol.atoms:
            if atom.symbol == "H":
                continue
            atom.mmff_type, atom.mmff_class = self._type_heavy(atom, mol, ring_info)

        # Pass 2: hydrogens (depend on parent heavy's mmff_type).
        for atom in mol.atoms:
            if atom.symbol != "H":
                continue
            atom.mmff_type, atom.mmff_class = self._type_hydrogen(atom, mol)

    def coverage(self, mol: Molecule) -> float:
        if not mol.atoms:
            return 1.0
        typed = sum(1 for a in mol.atoms if a.mmff_type > 0)
        return typed / len(mol.atoms)

    # ────────────────────────── Heavy dispatch ──────────────────────────

    def _type_heavy(self, atom: Atom, mol: Molecule,
                    ring_info: Dict[int, RingInfo]) -> Tuple[int, int]:
        """Return (mmff_type, mmff_class) for a non-hydrogen atom."""
        ri = ring_info[atom.index]
        sym = atom.symbol
        dispatch: Dict[str, Callable[[Atom, Molecule, RingInfo], Tuple[int, int]]] = {
            "C":  self._type_carbon,
            "N":  self._type_nitrogen,
            "O":  self._type_oxygen,
            "S":  self._type_sulfur,
            "P":  self._type_phosphorus,
            "F":  lambda a, m, r: (11, 11),
            "Cl": lambda a, m, r: (12, 12),
            "Br": lambda a, m, r: (13, 13),
            "I":  lambda a, m, r: (14, 14),
        }
        rule = dispatch.get(sym)
        if rule is None:
            # Unsupported heavy element (e.g., metals): mark unmatched.
            return (0, 0)
        return rule(atom, mol, ri)

    # ───────────────────────────── Carbon ─────────────────────────────

    def _type_carbon(self, atom: Atom, mol: Molecule, ri: RingInfo) -> Tuple[int, int]:
        """Carbon types: 1, 2, 3, 4, 20, 22, 30, 37, 41, 57."""
        neighbors = mol.get_neighbors(atom.index)
        neighbor_atoms = [mol.atoms[i] for i in neighbors]

        # Type 37: aromatic carbon (in aromatic ring).
        if ri.in_aromatic_ring:
            return (37, 37)

        # Type 22: C in cyclopropyl (3-ring).
        if 3 in ri.ring_sizes:
            # Type 22 is sp3 cyclopropyl. If double-bonded, may be type 30.
            hyb = atom.hybridization or "sp3"
            if hyb == "sp3":
                return (22, 22)

        # Determine double/triple bond environment.
        has_double = False
        has_triple = False
        double_neighbor_syms = []
        for n_idx in neighbors:
            for nb_other, bond in mol.get_neighbor_bonds(atom.index):
                if nb_other == n_idx:
                    if bond.order == 2:
                        has_double = True
                        double_neighbor_syms.append(mol.atoms[n_idx].symbol)
                    elif bond.order == 3:
                        has_triple = True

        # Type 4: alkyne (sp C with triple bond).
        if has_triple:
            return (4, 4)

        # Type 41: carboxylate C — C with two O neighbors, both with no H
        # and at least one bearing a formal negative charge.
        if atom.hybridization == "sp2" and has_double:
            o_neighbors = [na for na in neighbor_atoms if na.symbol == "O"]
            if len(o_neighbors) == 2:
                # Anion: at least one O has formal_charge < 0.
                has_anion = any(o.formal_charge < 0 for o in o_neighbors)
                # Both O have no H (i.e., no -OH; this distinguishes from -COOH).
                both_bare = all(
                    (o.num_explicit_h or 0) + o.num_implicit_h == 0
                    for o in o_neighbors
                )
                if has_anion or both_bare:
                    return (41, 3)

        # Type 57: guanidinium C — C bonded to 3 N's with at least one double bond.
        n_neighbors = [na for na in neighbor_atoms if na.symbol == "N"]
        if len(n_neighbors) == 3 and atom.hybridization == "sp2":
            return (57, 3)

        # Type 3: general carbonyl C — sp2 C with =O, =N, =S, or =P.
        if atom.hybridization == "sp2" and has_double:
            if any(sym in ("O", "N", "S", "P") for sym in double_neighbor_syms):
                return (3, 3)

        # Type 30: C=C in 4-ring.
        if 4 in ri.ring_sizes and has_double and "C" in double_neighbor_syms:
            return (30, 20)

        # Type 2: vinylic C (sp2 with C=C).
        if atom.hybridization == "sp2" and has_double and "C" in double_neighbor_syms:
            return (2, 2)

        # Type 20: C in 4-ring (sp3).
        if 4 in ri.ring_sizes:
            return (20, 20)

        # Type 1: alkyl C (default sp3).
        return (1, 1)

    # Other element rules — placeholders, filled by Tasks 11-13.

    def _type_nitrogen(self, atom: Atom, mol: Molecule, ri: RingInfo) -> Tuple[int, int]:
        """Nitrogen types: 8, 9, 10, 34, 38, 39, 40, 45, 67, 68."""
        neighbors = mol.get_neighbors(atom.index)
        neighbor_atoms = [mol.atoms[i] for i in neighbors]

        # Type 38: pyridine N (aromatic N with 2 neighbors, no H, formal_charge 0).
        if ri.in_aromatic_ring:
            n_h_explicit_neighbors = sum(1 for na in neighbor_atoms if na.symbol == "H")
            n_h = (atom.num_explicit_h or 0) + atom.num_implicit_h + n_h_explicit_neighbors
            if n_h == 0 and atom.formal_charge == 0:
                return (38, 38)
            # Type 39: pyrrole N (aromatic N with H).
            if n_h > 0:
                return (39, 39)
            # Aromatic N+ -> still type 38 class
            return (38, 38)

        # Determine double bonds + bond-order details.
        has_double_to = []  # list of (symbol) per double-bond
        for n_idx in neighbors:
            for nb_other, bond in mol.get_neighbor_bonds(atom.index):
                if nb_other == n_idx and bond.order == 2:
                    has_double_to.append(mol.atoms[n_idx].symbol)
                    break

        # Type 45: nitro N (N bonded to 2 O's, at least one double-bond O).
        o_neighbors = [na for na in neighbor_atoms if na.symbol == "O"]
        if len(o_neighbors) >= 2 and "O" in has_double_to:
            return (45, 45)

        # Type 67: N-oxide N (N with one =O and one or more C, formal_charge +0 or +1).
        if "O" in has_double_to and any(na.symbol == "C" for na in neighbor_atoms) \
                and len(o_neighbors) == 1:
            # Distinguish nitroso (N=O alone) - type 67 covers N-oxide of amines
            return (67, 8)

        # Type 10: amide N (sp3 N bonded to a carbonyl C - a C that has =O).
        for na in neighbor_atoms:
            if na.symbol == "C":
                # Check if that C has a =O
                for c_nb, c_bond in mol.get_neighbor_bonds(na.index):
                    if (mol.atoms[c_nb].symbol == "O"
                            and c_bond.order == 2):
                        return (10, 10)

        # Type 9: imine N (sp2 N with =C).
        if atom.hybridization == "sp2" and "C" in has_double_to:
            return (9, 9)

        # Type 40: aniline N (sp3 N bonded to an aromatic C).
        if any(na.symbol == "C" and na.is_aromatic for na in neighbor_atoms):
            return (40, 8)

        # Type 34 / 68: ammonium N+ / amine N+ (formal_charge >= 1, sp3).
        if atom.formal_charge >= 1:
            # Type 34: ammonium (4 single-bonded substituents, no double bonds)
            if len(neighbors) == 4 and not has_double_to:
                return (34, 8)
            return (68, 8)

        # Type 8: default amine N (sp3 amine).
        return (8, 8)

    def _type_oxygen(self, atom: Atom, mol: Molecule, ri: RingInfo) -> Tuple[int, int]:
        """Oxygen types: 6, 7, 32, 35, 49."""
        neighbors = mol.get_neighbors(atom.index)
        neighbor_atoms = [mol.atoms[i] for i in neighbors]

        # Type 32: anion O (formal_charge < 0).
        if atom.formal_charge < 0:
            return (32, 32)

        # Determine if double-bonded.
        has_double = False
        for n_idx in neighbors:
            for nb_other, bond in mol.get_neighbor_bonds(atom.index):
                if nb_other == n_idx and bond.order == 2:
                    has_double = True
                    break
            if has_double:
                break

        # Type 7: carbonyl O (double-bonded to C, N, S, or P).
        if has_double:
            return (7, 7)

        # Type 49: oxonium O (formal_charge > 0).
        if atom.formal_charge > 0:
            return (49, 6)

        # Count H attached: explicit H atom neighbors + num_explicit_h + implicit.
        n_h_explicit_neighbors = sum(1 for na in neighbor_atoms if na.symbol == "H")
        n_h = (atom.num_explicit_h or 0) + atom.num_implicit_h + n_h_explicit_neighbors
        # Type 35: hydroxyl O (H on O - alcohol, water, phenol).
        if n_h > 0:
            return (35, 6)

        # Type 6: ether/ester linker O (no double bonds, no H, no formal charge).
        return (6, 6)

    def _type_sulfur(self, atom: Atom, mol: Molecule, ri: RingInfo) -> Tuple[int, int]:
        """Sulfur types: 15, 16, 17, 18."""
        neighbors = mol.get_neighbors(atom.index)

        # Count =O double bonds.
        n_double_o = 0
        n_double_other = 0
        for n_idx in neighbors:
            for nb_other, bond in mol.get_neighbor_bonds(atom.index):
                if nb_other == n_idx and bond.order == 2:
                    if mol.atoms[n_idx].symbol == "O":
                        n_double_o += 1
                    else:
                        n_double_other += 1

        # Type 18: sulfone (two S=O double bonds).
        if n_double_o >= 2:
            return (18, 16)

        # Type 17: sulfoxide (one S=O double bond).
        if n_double_o == 1:
            return (17, 17)

        # Type 16: thiocarbonyl (S=C/N/P).
        if n_double_other > 0:
            return (16, 16)

        # Type 15: sulfide / thiol (no double bonds).
        return (15, 15)

    def _type_phosphorus(self, atom: Atom, mol: Molecule, ri: RingInfo) -> Tuple[int, int]:
        """Phosphorus types: 25, 26."""
        # Count =O double bonds.
        n_double_o = 0
        for n_idx in mol.get_neighbors(atom.index):
            for nb_other, bond in mol.get_neighbor_bonds(atom.index):
                if nb_other == n_idx and bond.order == 2:
                    if mol.atoms[n_idx].symbol == "O":
                        n_double_o += 1

        # Type 25: phosphate / phosphine oxide (P=O present).
        if n_double_o > 0:
            return (25, 25)

        # Type 26: phosphite (P with 3 single bonds, no =O).
        return (26, 26)

    def _type_hydrogen(self, atom: Atom, mol: Molecule) -> Tuple[int, int]:
        """Hydrogen types by parent heavy atom's mmff_type.

        H is always bonded to exactly one heavy atom in valid chemistry.
        Returns (type, class). class for hydrogens is mostly the same as
        the type number (H types collapse to fewer classes).
        """
        neighbors = mol.get_neighbors(atom.index)
        if not neighbors:
            return (0, 0)  # orphan H

        parent = mol.atoms[neighbors[0]]
        parent_type = parent.mmff_type

        # Map parent_type -> H type
        # See spec section 5.2 hydrogen types:
        # 5 H-on-C, 21 H-on-OH, 23 H-on-amine, 24 H-on-NH+,
        # 28 H-on-amide, 31 H-on-water, 36 H-on-guanidinium, 71 H-on-SH
        if parent_type in (1, 2, 3, 4, 20, 22, 30, 37, 41, 57):
            # H on C - type 5
            return (5, 5)
        if parent_type == 35:
            # H on alcohol-OH or water (we use type 21 for OH;
            # MMFF distinguishes water as type 31, but the typing
            # for water requires checking parent has 2 H neighbors).
            n_h_on_parent = sum(
                1 for n_idx in mol.get_neighbors(parent.index)
                if mol.atoms[n_idx].symbol == "H"
            )
            if n_h_on_parent >= 2:  # water-like (parent has >=2 H)
                return (31, 31)
            return (21, 21)
        if parent_type == 6:
            # H on ether O? Rare but possible (e.g., phosphate -OH).
            # Default to type 21 (alcohol H).
            return (21, 21)
        if parent_type == 32:
            # H on anion O - shouldn't exist, but fall back to type 21.
            return (21, 21)
        if parent_type == 10:
            # H on amide N - type 28
            return (28, 23)
        if parent_type in (8, 9, 38, 39, 40, 67):
            # H on amine / imine / aromatic / aniline / N-oxide - type 23
            return (23, 23)
        if parent_type in (34, 68):
            # H on N+ ammonium - type 24
            return (24, 23)
        if parent_type == 57:
            # H on guanidinium C - already covered above by 'H on C'
            # but H on guanidinium N is type 36 (rare case)
            return (5, 5)
        if parent_type in (15, 16, 17, 18):
            # H on S (thiol -SH) - type 71
            return (71, 23)
        if parent_type in (25, 26):
            # H on P (rare; phosphine) - fall back to type 5 (no proper P-H type
            # in our Phase 1 set; phosphate -OH is on O not P)
            return (5, 5)

        # Unknown parent type - leave H unmatched.
        return (0, 0)
