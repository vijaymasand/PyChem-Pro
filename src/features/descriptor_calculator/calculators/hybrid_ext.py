"""
Extended hybrid / drug-likeness descriptors.

Topological polar surface area from Ertl's fragment contributions, the Delaney
aqueous solubility estimate, McGowan volume, and the medicinal chemistry filters
that are normally reported next to them (Ghose, Veber, Egan, Muegge, rule of
three, lead-likeness). None of these need 3D coordinates except the ones that
build on the MLP logP.
"""
import numpy as np

from .hybrid import HybridCalculator
from .graph_utils import get_subgraph

# Ertl fragment contributions to the topological polar surface area (A^2).
# Key: (symbol, aromatic, n_neighbours, n_hydrogens, charge, max_bond_order)
TPSA_NITROGEN = {
    (3, 0, 0, 1): 3.24,    # N(-*)(-*)-*
    (2, 0, 0, 2): 12.36,   # N(-*)=*
    (1, 0, 0, 3): 23.79,   # N#*
    (3, 0, 0, 2): 11.68,   # N(-*)(=*)=* , nitro
    (2, 0, 0, 3): 13.60,   # N(=*)#*
    (2, 1, 0, 1): 12.03,   # NH(-*)-*
    (1, 1, 0, 2): 23.85,   # NH=*
    (1, 2, 0, 1): 26.02,   # NH2-*
    (4, 0, 1, 1): 0.00,    # N+(-*)(-*)(-*)-*
    (3, 0, 1, 2): 3.01,    # N+(-*)(-*)=*
    (2, 0, 1, 3): 4.36,    # N+(-*)#*
    (3, 1, 1, 1): 4.44,    # NH+(-*)(-*)-*
    (2, 1, 1, 2): 13.97,   # NH+(-*)=*
    (2, 2, 1, 1): 16.61,   # NH2+(-*)-*
    (1, 2, 1, 2): 25.59,   # NH2+=*
    (1, 3, 1, 1): 27.64,   # NH3+-*
}
TPSA_AROMATIC_NITROGEN = {
    (2, 0, 0): 12.89,      # n(:*):*
    (3, 0, 0): 4.41,       # n(:*)(:*):*
    (3, 0, 0, 'sub'): 4.93,  # n(-*)(:*):*
    (2, 1, 0): 15.79,      # nH(:*):*
    (3, 0, 1): 4.10,       # n+(:*)(:*):*
    (2, 1, 1): 14.14,      # nH+(:*):*
}
TPSA_OXYGEN = {
    (2, 0, 0, 1): 9.23,    # O(-*)-*
    (1, 0, 0, 2): 17.07,   # O=*
    (1, 1, 0, 1): 20.23,   # OH-*
    (1, 0, -1, 1): 23.06,  # O- -*
}
TPSA_SULFUR = {
    (2, 0, 1): 25.30,      # S(-*)-*
    (1, 0, 2): 32.09,      # S=*
    (3, 0, 2): 19.21,      # S(-*)(-*)=*
    (4, 0, 2): 8.38,       # S(-*)(-*)(=*)=*
    (1, 1, 1): 38.80,      # SH-*
}
TPSA_PHOSPHORUS = {
    (3, 0, 1): 13.59,      # P(-*)(-*)-*
    (2, 0, 2): 34.14,      # P(-*)=*
    (4, 0, 2): 9.81,       # P(-*)(-*)(-*)=*
    (3, 1, 2): 23.47,      # PH(-*)(-*)=*
}

# McGowan characteristic atomic volumes (cm3/mol / 100)
MCGOWAN_VOLUME = {
    'H': 8.71, 'C': 16.35, 'N': 14.39, 'O': 12.43, 'F': 10.48, 'Si': 26.83,
    'P': 24.87, 'S': 22.91, 'Cl': 20.95, 'Br': 26.21, 'I': 34.53, 'B': 18.32,
    'Se': 27.44,
}


class ExtendedHybridCalculator(HybridCalculator):
    """Drug-likeness descriptors, basic set plus filters and property models."""

    # ------------------------------------------------------------ helpers
    def _heavy_neighbors(self, molecule, idx, chosen):
        return [n for n in molecule.get_neighbors(idx)
                if n in chosen and molecule.atoms[n].symbol != 'H']

    def _h_count(self, molecule, idx, chosen):
        atom = molecule.atoms[idx]
        explicit = sum(1 for n in molecule.get_neighbors(idx)
                       if molecule.atoms[n].symbol == 'H')
        return getattr(atom, 'total_h', 0) + explicit

    def _max_bond_order(self, molecule, idx, chosen):
        order = 1
        for neighbor in molecule.get_neighbors(idx):
            if neighbor not in chosen:
                continue
            bond = molecule.get_bond_between(idx, neighbor)
            if bond is None:
                continue
            if bond.is_triple:
                order = max(order, 3)
            elif bond.is_double:
                order = max(order, 2)
        return order

    # --------------------------------------------------------------- TPSA
    def calc_tpsa(self, molecule, selection) -> float:
        """Topological polar surface area, Ertl fragment contributions (N and O)."""
        return self._tpsa(molecule, selection, include_sp=False)

    def calc_tpsa_with_sp(self, molecule, selection) -> float:
        """TPSA counting sulfur and phosphorus as polar as well."""
        return self._tpsa(molecule, selection, include_sp=True)

    def _tpsa(self, molecule, selection, include_sp=False) -> float:
        chosen = set(self.get_selected_atoms(molecule, selection))
        total = 0.0
        for idx in chosen:
            atom = molecule.atoms[idx]
            symbol = atom.symbol
            if symbol not in ('N', 'O') and not (include_sp and symbol in ('S', 'P')):
                continue
            neighbors = len(self._heavy_neighbors(molecule, idx, chosen))
            hydrogens = int(self._h_count(molecule, idx, chosen))
            charge = int(getattr(atom, 'formal_charge', 0) or 0)
            order = self._max_bond_order(molecule, idx, chosen)

            if symbol == 'N':
                if atom.is_aromatic:
                    total += self._aromatic_nitrogen_contribution(neighbors, hydrogens, charge)
                else:
                    total += TPSA_NITROGEN.get((neighbors, hydrogens, charge, order), 3.24)
            elif symbol == 'O':
                if atom.is_aromatic:
                    total += 13.14  # o(:*):*
                else:
                    total += TPSA_OXYGEN.get((neighbors, hydrogens, charge, order), 9.23)
            elif symbol == 'S':
                if atom.is_aromatic:
                    total += 28.24
                else:
                    total += TPSA_SULFUR.get((neighbors, hydrogens, order), 25.30)
            elif symbol == 'P':
                total += TPSA_PHOSPHORUS.get((neighbors, hydrogens, order), 13.59)
        return float(total)

    def _aromatic_nitrogen_contribution(self, neighbors, hydrogens, charge):
        if charge > 0:
            return 4.10 if neighbors >= 3 else 14.14
        if hydrogens:
            return 15.79
        if neighbors >= 3:
            return 4.93   # substituted, e.g. N-methyl imidazole
        return 12.89

    # ------------------------------------------------------ property models
    def calc_esol_logs(self, molecule, selection) -> float:
        """Delaney ESOL aqueous solubility, log10 of mol/L.

        LogS = 0.16 - 0.63*logP - 0.0062*MW + 0.066*RotB - 0.74*AromaticProportion
        """
        graph = get_subgraph(molecule, selection)
        if graph.n == 0:
            return 0.0
        logp = self.calc_lipophilicity(molecule, selection)
        mw = self.calc_molecular_weight(molecule, selection)
        rotatable = self.calc_rotatable_bonds(molecule, selection)
        aromatic = sum(1 for a in graph.atoms if a.is_aromatic) / graph.n
        return float(0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotatable - 0.74 * aromatic)

    def calc_mcgowan_volume(self, molecule, selection) -> float:
        """McGowan characteristic volume (cm3/mol/100), V = sum(atoms) - 6.56*bonds."""
        chosen = set(self.get_selected_atoms(molecule, selection))
        total = 0.0
        n_atoms = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            total += MCGOWAN_VOLUME.get(atom.symbol, 16.35)
            n_atoms += 1
            if atom.symbol != 'H':
                hydrogens = self._h_count(molecule, idx, chosen)
                total += hydrogens * MCGOWAN_VOLUME['H']
                n_atoms += hydrogens
        bonds = self.get_bond_count_in_selection(molecule, chosen)
        # implicit hydrogens add one bond each
        bonds += sum(self._h_count(molecule, i, chosen) for i in chosen
                     if molecule.atoms[i].symbol != 'H')
        return float(total - 6.56 * bonds)

    def calc_hydrophilic_factor(self, molecule, selection) -> float:
        """Todeschini hydrophilic factor, high for polar hydrogen rich molecules."""
        chosen = set(self.get_selected_atoms(molecule, selection))
        heavy = [i for i in chosen if molecule.atoms[i].symbol != 'H']
        n_sk = len(heavy)
        if n_sk == 0:
            return 0.0
        n_c = sum(1 for i in heavy if molecule.atoms[i].symbol == 'C')
        n_hy = 0
        for idx in heavy:
            atom = molecule.atoms[idx]
            if atom.symbol in ('N', 'O', 'S') and self._h_count(molecule, idx, chosen) > 0:
                n_hy += 1
        numerator = (1 + n_hy) * np.log2(1 + n_hy)
        numerator += n_c * (1.0 / n_sk) * np.log2(1.0 / n_sk)
        numerator += np.sqrt(n_hy / (n_sk ** 2))
        return float(numerator / np.log2(1.0 + n_sk))

    def calc_molecular_flexibility(self, molecule, selection) -> float:
        """Rotatable bonds per heavy atom."""
        graph = get_subgraph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return float(self.calc_rotatable_bonds(molecule, selection) / graph.n)

    def calc_ligand_efficiency_scale(self, molecule, selection) -> float:
        """LogP per heavy atom, the lipophilic efficiency scale used in lead work."""
        graph = get_subgraph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return float(self.calc_lipophilicity(molecule, selection) / graph.n)

    # ------------------------------------------------------------- filters
    def _profile(self, molecule, selection):
        """The property set every drug-likeness filter is expressed in."""
        graph = get_subgraph(molecule, selection)
        return {
            'mw': self.calc_molecular_weight(molecule, selection),
            'logp': self.calc_lipophilicity(molecule, selection),
            'hba': self.calc_lipinski_hba(molecule, selection),
            'hbd': self.calc_lipinski_hbd(molecule, selection),
            'tpsa': self.calc_tpsa(molecule, selection),
            'rotb': self.calc_rotatable_bonds(molecule, selection),
            'mr': self.calc_molar_refractivity(molecule, selection),
            'heavy': graph.n,
            'rings': len(graph.rings()),
            'carbons': sum(1 for a in graph.atoms if a.symbol == 'C'),
            'hetero': sum(1 for a in graph.atoms if a.symbol not in ('C', 'H')),
        }

    def calc_ghose_violations(self, molecule, selection) -> int:
        """Ghose filter: MW 160-480, logP -0.4 to 5.6, 20-70 atoms, MR 40-130."""
        p = self._profile(molecule, selection)
        violations = 0
        if not 160 <= p['mw'] <= 480:
            violations += 1
        if not -0.4 <= p['logp'] <= 5.6:
            violations += 1
        if not 20 <= p['heavy'] <= 70:
            violations += 1
        if not 40 <= p['mr'] <= 130:
            violations += 1
        return violations

    def calc_veber_violations(self, molecule, selection) -> int:
        """Veber oral bioavailability: rotatable bonds <= 10 and TPSA <= 140."""
        p = self._profile(molecule, selection)
        return int(p['rotb'] > 10) + int(p['tpsa'] > 140)

    def calc_egan_violations(self, molecule, selection) -> int:
        """Egan absorption egg: TPSA <= 131.6 and logP <= 5.88."""
        p = self._profile(molecule, selection)
        return int(p['tpsa'] > 131.6) + int(p['logp'] > 5.88)

    def calc_muegge_violations(self, molecule, selection) -> int:
        """Muegge pharmacophore filter, seven simultaneous criteria."""
        p = self._profile(molecule, selection)
        violations = 0
        if not 200 <= p['mw'] <= 600:
            violations += 1
        if not -2 <= p['logp'] <= 5:
            violations += 1
        if p['tpsa'] > 150:
            violations += 1
        if p['rings'] > 7:
            violations += 1
        if p['carbons'] <= 4:
            violations += 1
        if p['hetero'] <= 1:
            violations += 1
        if p['rotb'] > 15:
            violations += 1
        if p['hba'] > 10:
            violations += 1
        if p['hbd'] > 5:
            violations += 1
        return violations

    def calc_rule_of_three_violations(self, molecule, selection) -> int:
        """Congreve rule of three for fragment screening libraries."""
        p = self._profile(molecule, selection)
        violations = 0
        if p['mw'] > 300:
            violations += 1
        if p['logp'] > 3:
            violations += 1
        if p['hbd'] > 3:
            violations += 1
        if p['hba'] > 3:
            violations += 1
        if p['rotb'] > 3:
            violations += 1
        return violations

    def calc_lead_likeness_violations(self, molecule, selection) -> int:
        """Lead-likeness: MW 250-350, logP <= 3.5, rotatable bonds <= 7."""
        p = self._profile(molecule, selection)
        violations = 0
        if not 250 <= p['mw'] <= 350:
            violations += 1
        if p['logp'] > 3.5:
            violations += 1
        if p['rotb'] > 7:
            violations += 1
        return violations

    def calc_bbb_score(self, molecule, selection) -> float:
        """Blood-brain barrier likelihood on a 0-1 scale.

        Built from the properties that separate CNS from non-CNS drugs:
        TPSA below 90 A^2, logP between 1 and 4, MW below 400, at most 3
        H-bond donors and few rotatable bonds. """
        p = self._profile(molecule, selection)
        score = 0.0
        score += 0.3 if p['tpsa'] <= 90 else (0.15 if p['tpsa'] <= 120 else 0.0)
        score += 0.25 if 1.0 <= p['logp'] <= 4.0 else (0.1 if 0 <= p['logp'] <= 5 else 0.0)
        score += 0.2 if p['mw'] <= 400 else (0.1 if p['mw'] <= 500 else 0.0)
        score += 0.15 if p['hbd'] <= 3 else 0.0
        score += 0.1 if p['rotb'] <= 8 else 0.0
        return float(round(score, 3))

    def calc_oral_bioavailability_score(self, molecule, selection) -> float:
        """Fraction of the six common oral filters the molecule passes."""
        checks = [
            self.calc_lipinski_violations(molecule, selection) == 0,
            self.calc_ghose_violations(molecule, selection) == 0,
            self.calc_veber_violations(molecule, selection) == 0,
            self.calc_egan_violations(molecule, selection) == 0,
            self.calc_muegge_violations(molecule, selection) == 0,
            self.calc_lead_likeness_violations(molecule, selection) == 0,
        ]
        return float(sum(checks) / len(checks))

    def calc_polar_atom_fraction(self, molecule, selection) -> float:
        """N, O, S and P atoms as a fraction of the heavy atoms."""
        graph = get_subgraph(molecule, selection)
        if graph.n == 0:
            return 0.0
        polar = sum(1 for a in graph.atoms if a.symbol in ('N', 'O', 'S', 'P'))
        return float(polar / graph.n)

    def calc_tpsa_per_heavy_atom(self, molecule, selection) -> float:
        """TPSA normalised by size, comparable across molecular weights."""
        graph = get_subgraph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return float(self.calc_tpsa(molecule, selection) / graph.n)
