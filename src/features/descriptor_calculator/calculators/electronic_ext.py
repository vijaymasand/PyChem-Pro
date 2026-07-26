"""
Extended electronic descriptors.

Charged partial surface area (CPSA) family of Stanton & Jurs plus a set of
charge and electronegativity statistics. The CPSA descriptors need both
partial charges and 3D coordinates; charges are computed on the fly with the
in-house Gasteiger implementation when the molecule does not carry them yet.
"""
import numpy as np

from .electronic import ElectronicCalculator


class ExtendedElectronicCalculator(ElectronicCalculator):
    """Electronic descriptors, basic set plus the CPSA family."""

    # ------------------------------------------------------------ helpers
    def charges(self, molecule, selection):
        """Partial charges of the selection as an array (computed if missing)."""
        return np.array(self._atom_charges(molecule, selection))

    def _atom_sasa(self, molecule, selection):
        """Per-atom solvent accessible surface, empty when there are no coordinates."""
        atoms = [molecule.atoms[i] for i in self.get_selected_atoms(molecule, selection)]
        if not atoms or not all(getattr(a, 'has_coords', False) and a.x is not None
                                for a in atoms):
            return None
        from src.features.cheminformatics.services.geometry_utils import calculate_sasa_per_atom
        coords = np.array([[a.x, a.y, a.z] for a in atoms])
        symbols = [a.symbol for a in atoms]
        try:
            return np.asarray(calculate_sasa_per_atom(coords, symbols), dtype=float)
        except Exception:
            return None

    def _cpsa_parts(self, molecule, selection):
        """(charges, per atom surface) or (None, None) when unavailable."""
        surface = self._atom_sasa(molecule, selection)
        if surface is None:
            return None, None
        charges = self.charges(molecule, selection)
        if charges.size != surface.size:
            return None, None
        return charges, surface

    # -------------------------------------------------- charge statistics
    def calc_total_absolute_charge(self, molecule, selection) -> float:
        """Sum of absolute partial charges, the electronic 'size' of a molecule."""
        return float(np.sum(np.abs(self.charges(molecule, selection))))

    def calc_charge_variance(self, molecule, selection) -> float:
        charges = self.charges(molecule, selection)
        return float(np.var(charges)) if charges.size else 0.0

    def calc_sum_squared_charges(self, molecule, selection) -> float:
        charges = self.charges(molecule, selection)
        return float(np.sum(charges ** 2))

    def calc_submolecular_polarity(self, molecule, selection) -> float:
        """Largest difference between any two atomic charges (Todeschini SPP)."""
        charges = self.charges(molecule, selection)
        if charges.size == 0:
            return 0.0
        return float(np.max(charges) - np.min(charges))

    def calc_relative_positive_charge(self, molecule, selection) -> float:
        """Most positive charge divided by the total positive charge."""
        charges = self.charges(molecule, selection)
        positive = charges[charges > 0]
        total = float(np.sum(positive))
        if total <= 0:
            return 0.0
        return float(np.max(positive) / total)

    def calc_relative_negative_charge(self, molecule, selection) -> float:
        charges = self.charges(molecule, selection)
        negative = charges[charges < 0]
        total = float(np.sum(np.abs(negative)))
        if total <= 0:
            return 0.0
        return float(np.max(np.abs(negative)) / total)

    def calc_mean_electronegativity(self, molecule, selection) -> float:
        values = [getattr(molecule.atoms[i].element, 'electronegativity', 0.0) or 0.0
                  for i in self.get_selected_atoms(molecule, selection)]
        return float(np.mean(values)) if values else 0.0

    def calc_electronegativity_variance(self, molecule, selection) -> float:
        values = [getattr(molecule.atoms[i].element, 'electronegativity', 0.0) or 0.0
                  for i in self.get_selected_atoms(molecule, selection)]
        return float(np.var(values)) if values else 0.0

    def calc_max_bond_polarity(self, molecule, selection) -> float:
        """Largest electronegativity difference across a bond."""
        chosen = self.get_selected_set(selection)
        best = 0.0
        for bond in molecule.bonds:
            i, j = bond.begin_atom_idx, bond.end_atom_idx
            if i not in chosen or j not in chosen:
                continue
            ei = getattr(molecule.atoms[i].element, 'electronegativity', 0.0) or 0.0
            ej = getattr(molecule.atoms[j].element, 'electronegativity', 0.0) or 0.0
            best = max(best, abs(ei - ej))
        return float(best)

    def calc_topological_polar_charge(self, molecule, selection) -> float:
        """Absolute charge carried by N, O, S and P atoms."""
        polar = {'N', 'O', 'S', 'P'}
        self.ensure_charges(molecule)
        total = 0.0
        for idx in self.get_selected_atoms(molecule, selection):
            atom = molecule.atoms[idx]
            if atom.symbol in polar:
                total += abs(getattr(atom, 'partial_charge', 0.0) or 0.0)
        return float(total)

    # ------------------------------------------------------------- CPSA
    def calc_ppsa1(self, molecule, selection) -> float:
        """Partial positive surface area: surface of the positively charged atoms."""
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        return float(np.sum(surface[charges > 0]))

    def calc_pnsa1(self, molecule, selection) -> float:
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        return float(np.sum(surface[charges < 0]))

    def calc_ppsa2(self, molecule, selection) -> float:
        """Total charge weighted positive surface area."""
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        positive = charges > 0
        return float(np.sum(charges[positive]) * np.sum(surface[positive]))

    def calc_pnsa2(self, molecule, selection) -> float:
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        negative = charges < 0
        return float(np.sum(charges[negative]) * np.sum(surface[negative]))

    def calc_ppsa3(self, molecule, selection) -> float:
        """Atomic charge weighted positive surface area."""
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        positive = charges > 0
        return float(np.sum(charges[positive] * surface[positive]))

    def calc_pnsa3(self, molecule, selection) -> float:
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        negative = charges < 0
        return float(np.sum(charges[negative] * surface[negative]))

    def calc_dpsa1(self, molecule, selection) -> float:
        """Difference between the positive and the negative surface area."""
        return self.calc_ppsa1(molecule, selection) - self.calc_pnsa1(molecule, selection)

    def calc_fpsa1(self, molecule, selection) -> float:
        """Fractional positive surface area, PPSA1 over the total surface."""
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        total = float(np.sum(surface))
        if total <= 0:
            return 0.0
        return float(np.sum(surface[charges > 0]) / total)

    def calc_fnsa1(self, molecule, selection) -> float:
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        total = float(np.sum(surface))
        if total <= 0:
            return 0.0
        return float(np.sum(surface[charges < 0]) / total)

    def calc_wpsa1(self, molecule, selection) -> float:
        """Surface weighted positive surface area, PPSA1 * total surface / 1000."""
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        return float(np.sum(surface[charges > 0]) * np.sum(surface) / 1000.0)

    def calc_wnsa1(self, molecule, selection) -> float:
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None:
            return 0.0
        return float(np.sum(surface[charges < 0]) * np.sum(surface) / 1000.0)

    def calc_rpcs(self, molecule, selection) -> float:
        """Relative positive charged surface area (Jurs RPCS)."""
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None or charges.size == 0:
            return 0.0
        positive = charges > 0
        if not positive.any():
            return 0.0
        total_positive = float(np.sum(charges[positive]))
        if total_positive <= 0:
            return 0.0
        most_positive = int(np.argmax(charges))
        return float(surface[most_positive] * np.max(charges) / total_positive)

    def calc_rncs(self, molecule, selection) -> float:
        charges, surface = self._cpsa_parts(molecule, selection)
        if charges is None or charges.size == 0:
            return 0.0
        negative = charges < 0
        if not negative.any():
            return 0.0
        total_negative = float(np.sum(np.abs(charges[negative])))
        if total_negative <= 0:
            return 0.0
        most_negative = int(np.argmin(charges))
        return float(surface[most_negative] * abs(np.min(charges)) / total_negative)
