"""
Hybrid/drug-like descriptor calculations.
Includes: lipophilicity, polarizability, molar refractivity, Lipinski rules, drug-likeness scores.
"""
from . import BaseCalculator


class HybridCalculator(BaseCalculator):
    """Calculator for hybrid/drug-like molecular descriptors."""

    def calc_lipophilicity(self, molecule, selection) -> float:
        """Calculate molecular lipophilicity (logP) using MLP method."""
        from ...cheminformatics.services.lipophilicity_service import calculate_logp

        # The MLP method integrates over the solvent accessible surface, so the
        # molecule needs a geometry; one is generated when it has none.
        if not self.ensure_coordinates(molecule):
            return 0.0
        # If selection is not the whole molecule, we could subset,
        # but logP is generally a molecular property.
        # For consistency with other calculators, we'll calculate it for the full molecule
        # unless it's a very specific fragment selection.
        return calculate_logp(molecule)

    def calc_polarizability(self, molecule, selection) -> float:
        """Calculate molecular polarizability."""
        atomic_polarizabilities = {'H': 0.667, 'C': 1.76, 'N': 1.10, 'O': 0.802,
                                  'F': 0.557, 'Cl': 2.18, 'Br': 3.05, 'I': 5.35}

        return sum(atomic_polarizabilities.get(molecule.atoms[idx].symbol, 1.0)
                  for idx in self.get_selected_atoms(molecule, selection))

    def calc_molar_refractivity(self, molecule, selection) -> float:
        """Calculate molar refractivity."""
        atomic_refractivities = {'H': 1.463, 'C': 2.491, 'N': 2.743, 'O': 1.529,
                               'F': 0.923, 'Cl': 5.607, 'Br': 8.865, 'I': 13.900}

        return sum(atomic_refractivities.get(molecule.atoms[idx].symbol, 2.0)
                  for idx in self.get_selected_atoms(molecule, selection))


    def calc_lipinski_hba(self, molecule, selection) -> int:
        """Calculate Lipinski H-bond acceptors (N + O count)."""
        return (self.count_atoms_by_symbol(molecule, selection, 'N') +
                self.count_atoms_by_symbol(molecule, selection, 'O'))

    def calc_lipinski_hbd(self, molecule, selection) -> int:
        """Calculate Lipinski H-bond donors (OH + NH count)."""
        from ...cheminformatics.services.atom_properties import AtomPropertyAnalyzer
        analyzer = AtomPropertyAnalyzer(molecule)
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if analyzer.is_hbond_donor(idx))

    def calc_lipinski_violations(self, molecule, selection) -> int:
        """Calculate number of Lipinski Rule of 5 violations."""
        violations = 0

        # Molecular weight > 500
        mw = self.calc_molecular_weight(molecule, selection)
        if mw > 500:
            violations += 1

        # LogP > 5
        logp = self.calc_lipophilicity(molecule, selection)
        if logp > 5:
            violations += 1

        # HBA > 10
        hba = self.calc_lipinski_hba(molecule, selection)
        if hba > 10:
            violations += 1

        # HBD > 5
        hbd = self.calc_lipinski_hbd(molecule, selection)
        if hbd > 5:
            violations += 1

        return violations

    def calc_veber_tpsa(self, molecule, selection) -> float:
        """Calculate Veber polar surface area."""
        from ...cheminformatics.services.atom_properties import AtomPropertyAnalyzer
        analyzer = AtomPropertyAnalyzer(molecule)
        polar_set = set(analyzer.POLAR_ATOMS)
        atomic_areas = {'O': 17.0, 'N': 12.0, 'S': 25.0, 'P': 25.0}

        return sum(atomic_areas.get(molecule.atoms[idx].symbol, 15.0)
                  for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol in polar_set)

    def calc_drug_likeness_score(self, molecule, selection) -> float:
        """Calculate composite drug-likeness score (0-1 scale)."""
        score = 1.0

        # Lipinski violations
        lipinski_violations = self.calc_lipinski_violations(molecule, selection)
        if lipinski_violations > 0:
            score -= lipinski_violations * 0.15

        # Molecular weight preference (150-500)
        mw = self.calc_molecular_weight(molecule, selection)
        if mw < 150 or mw > 500:
            score -= 0.1

        # LogP preference (0-5)
        logp = self.calc_lipophilicity(molecule, selection)
        if logp < 0 or logp > 5:
            score -= 0.1

        # Rotatable bonds (<10 preferred)
        rotb = self.calc_rotatable_bonds(molecule, selection)
        if rotb > 10:
            score -= 0.05

        return max(0.0, min(1.0, score))

    def calc_synthetic_accessibility(self, molecule, selection) -> float:
        """Calculate synthetic accessibility score (1-10 scale, lower is easier)."""
        score = 5.0  # Base score

        # Complexity factors
        rings = molecule.find_rings()
        selected_set = self.get_selected_set(selection)
        ring_count = sum(1 for ring in rings if all(idx in selected_set for idx in ring))
        score += ring_count * 0.5

        rotb = self.calc_rotatable_bonds(molecule, selection)
        score += min(rotb * 0.1, 1.0)

        # Heteroatoms increase complexity
        hetero = sum(1 for idx in self.get_selected_atoms(molecule, selection)
                    if molecule.atoms[idx].symbol not in ['C', 'H'])
        score += min(hetero * 0.05, 1.0)

        return max(1.0, min(10.0, score))

    def calc_fsp3(self, molecule, selection) -> float:
        """Calculate fraction of sp3 hybridized carbons (Fsp3)."""
        total_c = self.count_atoms_by_symbol(molecule, selection, 'C')
        if total_c == 0:
            return 0.0

        molecule.assign_hybridization()
        sp3_c = sum(1 for idx in self.get_selected_atoms(molecule, selection)
                   if molecule.atoms[idx].symbol == 'C' and molecule.atoms[idx].hybridization == 'sp3')
        return sp3_c / total_c

    def calc_molecular_weight(self, molecule, selection) -> float:
        """Calculate molecular weight."""
        atomic_weights = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
                        'F': 18.998, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904,
                        'S': 32.066, 'P': 30.974, 'Si': 28.086}

        return sum(atomic_weights.get(molecule.atoms[idx].symbol, 0.0)
                  for idx in self.get_selected_atoms(molecule, selection))

    def calc_rotatable_bonds(self, molecule, selection) -> int:
        """Calculate number of rotatable bonds."""
        selected_set = self.get_selected_set(selection)
        count = 0
        for bond in molecule.bonds:
            if (bond.begin_atom_idx in selected_set and
                bond.end_atom_idx in selected_set and
                bond.is_single and not bond.is_in_ring and
                not (molecule.atoms[bond.begin_atom_idx].symbol in ['H', 'F', 'Cl', 'Br', 'I'] or
                     molecule.atoms[bond.end_atom_idx].symbol in ['H', 'F', 'Cl', 'Br', 'I'])):
                count += 1
        return count

    def calc_hetero_atom_count(self, molecule, selection) -> int:
        """Calculate hetero atom count (non-C, non-H)."""
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol not in ['C', 'H'])
