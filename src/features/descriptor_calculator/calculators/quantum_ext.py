"""
Quantum chemical descriptors from simple Huckel molecular orbital theory.

The base QuantumCalculator returned the same constants for every molecule,
which made every quantum descriptor useless. This module solves the Huckel
secular determinant of the pi system instead, which is cheap (a symmetric
eigenvalue problem the size of the pi system) and gives frontier orbital
energies that actually respond to conjugation, ring size and heteroatoms.

Orbital energies follow E = alpha + x * beta with the parametrisation fitted to
ionisation potentials (alpha = -6.15 eV, beta = -3.32 eV), which puts the HOMO
of benzene at -9.47 eV against an experimental IP of 9.24 eV and naphthalene at
-8.20 eV against 8.14 eV. Heteroatoms use the h_X / k_XY parameters of
Streitwieser. Pi electron and delocalisation energies are reported in units of
beta, the scale free convention, because no single beta reproduces both
spectroscopic and thermochemical quantities.

Molecules without a pi system have no Huckel frontier orbitals. For those the
frontier energies are estimated from the valence orbital ionisation potentials
of the atoms present, which keeps the descriptors molecule dependent instead of
constant. Values are labelled as estimates in the descriptor documentation.
"""
import numpy as np

from .quantum import QuantumCalculator

ALPHA = -6.15   # eV, Coulomb integral of a carbon 2p orbital (IP fitted)
BETA = -3.32    # eV, carbon-carbon resonance integral (IP fitted)

# Streitwieser heteroatom parameters: alpha_X = alpha + h_X * beta
H_PARAMS = {
    'C': 0.0, 'N': 0.5, 'N+': 2.0, 'O': 1.0, 'O2': 2.0, 'S': 1.0,
    'F': 3.0, 'Cl': 2.0, 'Br': 1.5, 'I': 1.4, 'P': 0.4, 'B': -0.45,
}
# beta_CX = k_CX * beta
K_PARAMS = {
    'C': 1.0, 'N': 0.8, 'O': 0.8, 'S': 0.7, 'F': 0.7,
    'Cl': 0.4, 'Br': 0.3, 'I': 0.3, 'P': 0.6, 'B': 0.7,
}

# Saturated molecules have no Huckel pi system. Their HOMO is the highest lone
# pair (or the C-C/C-H sigma manifold when there is none) and their LUMO the
# corresponding sigma*. These are the typical ionisation energies (eV) of those
# localised orbitals in organic molecules; the smallest one present wins,
# because that is the electron a molecule gives up first.
SIGMA_IP = {
    'C': 10.9, 'H': 10.9, 'N': 9.0, 'O': 10.3, 'S': 8.7, 'P': 9.5,
    'F': 12.5, 'Cl': 11.0, 'Br': 10.5, 'I': 9.5, 'Si': 10.0, 'B': 10.5,
    'Se': 8.4, 'As': 9.2,
}
# electron affinity (eV) of the matching sigma* orbital; the largest one present
# wins, since that is the orbital an electron enters first
SIGMA_EA = {
    'C': -1.3, 'H': -1.3, 'N': -1.0, 'O': -0.9, 'S': -0.4, 'P': -0.6,
    'F': -0.5, 'Cl': 0.2, 'Br': 0.4, 'I': 0.6, 'Si': -0.5, 'B': -0.2,
    'Se': 0.5, 'As': -0.3,
}


class HuckelResult:
    """Eigenvalues and occupation of a Huckel pi system.

    `x` holds the eigenvalues of the topological matrix in descending order
    (most bonding first); `energies` the same levels as eV, ascending. """

    __slots__ = ('x', 'energies', 'n_electrons', 'n_orbitals', 'systems')

    def __init__(self, x, n_electrons, systems=0):
        self.x = np.sort(np.asarray(x, dtype=float))[::-1] if len(x) else np.array([])
        self.energies = ALPHA + self.x * BETA   # beta < 0, so this is ascending
        self.n_electrons = n_electrons
        self.n_orbitals = len(self.x)
        self.systems = systems

    @property
    def homo(self):
        occupied = self.n_electrons // 2
        if occupied <= 0 or occupied > self.n_orbitals:
            return None
        return float(self.energies[occupied - 1])

    @property
    def lumo(self):
        occupied = self.n_electrons // 2
        if occupied < 0 or occupied >= self.n_orbitals:
            return None
        return float(self.energies[occupied])

    @property
    def pi_energy_beta(self):
        """Total pi electron energy above n*alpha, in units of beta."""
        occupied = min(self.n_electrons // 2, self.n_orbitals)
        if occupied <= 0:
            return 0.0
        energy = 2.0 * float(np.sum(self.x[:occupied]))
        if self.n_electrons % 2 and occupied < self.n_orbitals:
            energy += float(self.x[occupied])
        return energy


class ExtendedQuantumCalculator(QuantumCalculator):
    """Quantum descriptors backed by a Huckel MO calculation."""

    # ---------------------------------------------------------- pi system
    def _pi_atoms(self, molecule, selection):
        """Atoms that contribute a p orbital to the pi system.

        Included are atoms carrying a double or triple bond or flagged as
        aromatic, plus N/O/S neighbours of such atoms, which donate a lone
        pair (this is what makes pyrrole and furan come out right). """
        chosen = {i for i in selection.atom_indices if i < len(molecule.atoms)}
        unsaturated = set()
        for bond in molecule.bonds:
            i, j = bond.begin_atom_idx, bond.end_atom_idx
            if i not in chosen or j not in chosen:
                continue
            if bond.is_double or bond.is_triple or bond.is_aromatic:
                unsaturated.add(i)
                unsaturated.add(j)
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.is_aromatic:
                unsaturated.add(idx)

        donors = set()
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.symbol not in ('N', 'O', 'S'):
                continue
            if idx in unsaturated:
                continue
            if any(n in unsaturated for n in molecule.get_neighbors(idx)):
                donors.add(idx)
        return sorted(unsaturated | donors), unsaturated, donors

    def _pi_electron_count(self, molecule, indices, unsaturated, donors):
        """Electrons the pi system holds: one per unsaturated atom, two per donor."""
        electrons = 0
        for idx in indices:
            atom = molecule.atoms[idx]
            if idx in donors:
                electrons += 2
            elif atom.symbol in ('N', 'O', 'S') and atom.is_aromatic:
                # aromatic heteroatoms with three connections donate a lone pair
                if len(molecule.get_neighbors(idx)) >= 3:
                    electrons += 2
                else:
                    electrons += 1
            else:
                electrons += 1
            electrons -= getattr(atom, 'formal_charge', 0) or 0
        return max(electrons, 0)

    def huckel(self, molecule, selection) -> HuckelResult:
        """Solve the Huckel secular problem for the pi system of a selection."""
        indices, unsaturated, donors = self._pi_atoms(molecule, selection)
        if not indices:
            return HuckelResult(np.array([]), 0, 0)

        position = {a: k for k, a in enumerate(indices)}
        n = len(indices)
        matrix = np.zeros((n, n))
        for k, idx in enumerate(indices):
            atom = molecule.atoms[idx]
            key = atom.symbol
            if key == 'O' and len(molecule.get_neighbors(idx)) >= 2:
                key = 'O2'   # ether/furan type oxygen holds its pair more tightly
            matrix[k, k] = H_PARAMS.get(key, H_PARAMS.get(atom.symbol, 0.0))

        edges = 0
        for bond in molecule.bonds:
            i, j = bond.begin_atom_idx, bond.end_atom_idx
            if i not in position or j not in position:
                continue
            si = molecule.atoms[i].symbol
            sj = molecule.atoms[j].symbol
            k_ij = K_PARAMS.get(si, 1.0) if si != 'C' else K_PARAMS.get(sj, 1.0)
            matrix[position[i], position[j]] = k_ij
            matrix[position[j], position[i]] = k_ij
            edges += 1

        if edges == 0:
            return HuckelResult(np.array([]), 0, 0)

        # eigenvalues of the topological matrix, E = alpha + x*beta
        x_values = np.linalg.eigvalsh(matrix)
        electrons = self._pi_electron_count(molecule, indices, unsaturated, donors)
        systems = self._count_pi_systems(molecule, indices)
        return HuckelResult(x_values, electrons, systems)

    def _count_pi_systems(self, molecule, indices):
        """Number of separate conjugated systems."""
        remaining = set(indices)
        systems = 0
        while remaining:
            systems += 1
            stack = [remaining.pop()]
            while stack:
                current = stack.pop()
                for neighbor in molecule.get_neighbors(current):
                    if neighbor in remaining:
                        remaining.discard(neighbor)
                        stack.append(neighbor)
        return systems

    # ------------------------------------------------------- frontier MOs
    def _fallback_frontier(self, molecule, selection):
        """Localised orbital estimate for molecules without a pi system.

        The HOMO is taken as the most easily ionised lone pair or sigma bond
        present and the LUMO as the most accessible sigma*, which keeps both
        molecule dependent (propanol ionises more easily than butane). """
        symbols = [molecule.atoms[i].symbol for i in selection.atom_indices
                   if i < len(molecule.atoms)]
        heavy = [s for s in symbols if s != 'H'] or symbols
        if not heavy:
            return -10.9, 1.3
        homo = -min(SIGMA_IP.get(s, 10.9) for s in heavy)
        lumo = -max(SIGMA_EA.get(s, -1.3) for s in heavy)
        return homo, lumo

    def calc_homo_energy(self, molecule, selection) -> float:
        result = self.huckel(molecule, selection)
        homo = result.homo
        if homo is None:
            return float(self._fallback_frontier(molecule, selection)[0])
        return float(homo)

    def calc_lumo_energy(self, molecule, selection) -> float:
        result = self.huckel(molecule, selection)
        lumo = result.lumo
        if lumo is None:
            return float(self._fallback_frontier(molecule, selection)[1])
        return float(lumo)

    # --------------------------------------------------- derived indices
    def calc_ionization_potential(self, molecule, selection) -> float:
        """Koopmans ionisation potential, IP = -E(HOMO)."""
        return -self.calc_homo_energy(molecule, selection)

    def calc_electron_affinity(self, molecule, selection) -> float:
        """Koopmans electron affinity, EA = -E(LUMO)."""
        return -self.calc_lumo_energy(molecule, selection)

    def calc_mulliken_electronegativity(self, molecule, selection) -> float:
        """chi = (IP + EA) / 2, the negative of the chemical potential."""
        return (self.calc_ionization_potential(molecule, selection) +
                self.calc_electron_affinity(molecule, selection)) / 2.0

    def calc_electrophilicity_index(self, molecule, selection) -> float:
        """Parr electrophilicity omega = mu^2 / (2 * eta)."""
        hardness = self.calc_chemical_hardness(molecule, selection)
        if abs(hardness) < 1e-6:
            return 0.0
        potential = self.calc_chemical_potential(molecule, selection)
        return float(potential ** 2 / (2.0 * hardness))

    def calc_nucleophilicity_index(self, molecule, selection) -> float:
        """Reciprocal of the electrophilicity index."""
        omega = self.calc_electrophilicity_index(molecule, selection)
        return float(1.0 / omega) if abs(omega) > 1e-9 else 0.0

    def calc_electroaccepting_power(self, molecule, selection) -> float:
        """omega+ = (IP + 3*EA)^2 / (16*(IP - EA)), Gazquez electron accepting power."""
        ip = self.calc_ionization_potential(molecule, selection)
        ea = self.calc_electron_affinity(molecule, selection)
        denominator = 16.0 * (ip - ea)
        if abs(denominator) < 1e-6:
            return 0.0
        return float((ip + 3 * ea) ** 2 / denominator)

    def calc_electrodonating_power(self, molecule, selection) -> float:
        """omega- = (3*IP + EA)^2 / (16*(IP - EA))."""
        ip = self.calc_ionization_potential(molecule, selection)
        ea = self.calc_electron_affinity(molecule, selection)
        denominator = 16.0 * (ip - ea)
        if abs(denominator) < 1e-6:
            return 0.0
        return float((3 * ip + ea) ** 2 / denominator)

    def calc_pi_electron_count(self, molecule, selection) -> int:
        return int(self.huckel(molecule, selection).n_electrons)

    def calc_pi_system_size(self, molecule, selection) -> int:
        """Number of atoms taking part in the pi system."""
        return int(self.huckel(molecule, selection).n_orbitals)

    def calc_pi_system_count(self, molecule, selection) -> int:
        """Number of separate conjugated fragments."""
        return int(self.huckel(molecule, selection).systems)

    def calc_total_pi_energy(self, molecule, selection) -> float:
        """Total Huckel pi electron energy in units of beta (8.0 for benzene)."""
        return float(self.huckel(molecule, selection).pi_energy_beta)

    def calc_delocalization_energy(self, molecule, selection) -> float:
        """Resonance energy in units of beta: the pi energy an equivalent set of
        isolated double bonds would not have. Benzene gives 2.000, butadiene
        0.472 and naphthalene 3.683, the textbook Huckel values. """
        result = self.huckel(molecule, selection)
        if result.n_orbitals == 0:
            return 0.0
        pairs = result.n_electrons // 2
        return float(result.pi_energy_beta - 2.0 * pairs)

    def calc_homo_lumo_gap_per_atom(self, molecule, selection) -> float:
        """Gap normalised by the size of the pi system, tracks conjugation length."""
        result = self.huckel(molecule, selection)
        gap = self.calc_homo_lumo_gap(molecule, selection)
        if result.n_orbitals == 0:
            return gap
        return float(gap / result.n_orbitals)

    def calc_absolute_hardness_ratio(self, molecule, selection) -> float:
        """Hardness over electronegativity, a scale free reactivity ratio."""
        electronegativity = self.calc_mulliken_electronegativity(molecule, selection)
        if abs(electronegativity) < 1e-6:
            return 0.0
        return float(self.calc_chemical_hardness(molecule, selection) / electronegativity)
