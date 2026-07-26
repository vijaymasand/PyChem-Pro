"""
Constitutional descriptor calculations.
Includes: elemental counts, bond counts, ring counts, functional groups, ratios, proportions.
"""
from typing import List, Set
from . import BaseCalculator


class ConstitutionalCalculator(BaseCalculator):
    """Calculator for constitutional (1D) molecular descriptors."""

    def calc_molecular_weight(self, molecule, selection) -> float:
        """Calculate molecular weight.

        Uses the element masses from the periodic table (so every element
        contributes, not just a hard coded handful) and adds the hydrogens that
        a structure from SMILES only carries implicitly. """
        total_weight = 0.0
        for idx in self.get_selected_atoms(molecule, selection):
            atom = molecule.atoms[idx]
            total_weight += getattr(atom, 'mass', 0.0)
            if atom.symbol != 'H':
                total_weight += self.count_hydrogens(molecule, idx) * 1.008
        return total_weight

    def count_hydrogens(self, molecule, idx) -> int:
        """Hydrogens on an atom, implicit and explicit alike."""
        atom = molecule.atoms[idx]
        explicit = sum(1 for n in molecule.get_neighbors(idx)
                       if n < len(molecule.atoms) and molecule.atoms[n].symbol == 'H')
        return int(getattr(atom, 'total_h', 0) or 0) + explicit

    def calc_atom_count(self, molecule, selection) -> int:
        """Calculate total atom count."""
        return len(selection.atom_indices)

    def calc_heavy_atom_count(self, molecule, selection) -> int:
        """Calculate heavy atom count."""
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol != 'H')

    def calc_bond_count(self, molecule, selection) -> int:
        """Calculate bond count in selection."""
        return self.get_bond_count_in_selection(molecule, self.get_selected_set(selection))

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

    def calc_ring_count(self, molecule, selection) -> int:
        """Calculate number of rings in selection."""
        rings = molecule.find_rings()
        selected_set = self.get_selected_set(selection)
        return sum(1 for ring in rings if all(idx in selected_set for idx in ring))

    def calc_aromatic_ring_count(self, molecule, selection) -> int:
        """Calculate number of aromatic rings in selection."""
        rings = molecule.find_rings()
        selected_set = self.get_selected_set(selection)
        count = 0
        for ring in rings:
            if all(idx in selected_set for idx in ring):
                if all(molecule.atoms[idx].is_aromatic for idx in ring):
                    count += 1
        return count

    def calc_h_donor_count(self, molecule, selection) -> int:
        """Calculate hydrogen bond donor count."""
        from ...cheminformatics.services.atom_properties import AtomPropertyAnalyzer
        analyzer = AtomPropertyAnalyzer(molecule)
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if analyzer.is_hbond_donor(idx))

    def calc_h_acceptor_count(self, molecule, selection) -> int:
        """Calculate hydrogen bond acceptor count."""
        from ...cheminformatics.services.atom_properties import AtomPropertyAnalyzer
        analyzer = AtomPropertyAnalyzer(molecule)
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if analyzer.is_hbond_acceptor(idx))

    def calc_lipophilic_count(self, molecule, selection) -> int:
        """Calculate lipophilic atom count."""
        from ...cheminformatics.services.atom_properties import AtomPropertyAnalyzer
        analyzer = AtomPropertyAnalyzer(molecule)
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if analyzer.is_lipophilic(idx))

    # Elemental counts
    def calc_carbon_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'C')

    def calc_nitrogen_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'N')

    def calc_oxygen_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'O')

    def calc_sulfur_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'S')

    def calc_phosphorus_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'P')

    def calc_halogen_count(self, molecule, selection) -> int:
        """Calculate halogen atom count."""
        halogens = {'F', 'Cl', 'Br', 'I'}
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol in halogens)

    def calc_hydrogen_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'H')

    # Individual halogen counts
    def calc_fluorine_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'F')

    def calc_chlorine_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'Cl')

    def calc_bromine_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'Br')

    def calc_iodine_count(self, molecule, selection) -> int:
        return self.count_atoms_by_symbol(molecule, selection, 'I')

    # Bond counts
    def calc_single_bond_count(self, molecule, selection) -> int:
        self.ensure_perception(molecule)
        selected_set = self.get_selected_set(selection)
        return sum(1 for bond in molecule.bonds
                  if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set and bond.is_single)

    def calc_double_bond_count(self, molecule, selection) -> int:
        self.ensure_perception(molecule)
        selected_set = self.get_selected_set(selection)
        return sum(1 for bond in molecule.bonds
                  if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set and bond.is_double)

    def calc_triple_bond_count(self, molecule, selection) -> int:
        self.ensure_perception(molecule)
        selected_set = self.get_selected_set(selection)
        return sum(1 for bond in molecule.bonds
                  if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set and bond.is_triple)

    def calc_aromatic_bond_count(self, molecule, selection) -> int:
        """Calculate aromatic bond count by identifying bonds within aromatic rings."""
        selected_set = self.get_selected_set(selection)
        aromatic_bonds = set()
        
        # 1. First, count bonds explicitly marked as aromatic
        for bond in molecule.bonds:
            if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set:
                if bond.is_aromatic:
                    aromatic_bonds.add(tuple(sorted((bond.begin_atom_idx, bond.end_atom_idx))))
        
        # 2. Fallback: Identify bonds between aromatic atoms in the same ring
        # This is very robust for SDF/MOL files where bond orders might be Kekulized
        rings = molecule.find_rings()
        for ring in rings:
            # If all atoms in the ring are marked aromatic, count all its bonds as aromatic
            if all(idx < len(molecule.atoms) and molecule.atoms[idx].is_aromatic for idx in ring):
                for i in range(len(ring)):
                    a1 = ring[i]
                    a2 = ring[(i+1) % len(ring)]
                    if a1 in selected_set and a2 in selected_set:
                         aromatic_bonds.add(tuple(sorted((a1, a2))))
                             
        return len(aromatic_bonds)

    # Hybridization-dependent carbon counts
    def calc_sp3_carbon_count(self, molecule, selection) -> int:
        """Calculate sp3 hybridized carbon atom count."""
        molecule.assign_hybridization()
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'C' and molecule.atoms[idx].hybridization == 'sp3')

    def calc_sp2_carbon_count(self, molecule, selection) -> int:
        """Calculate sp2 hybridized carbon atom count."""
        molecule.assign_hybridization()
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'C' and molecule.atoms[idx].hybridization == 'sp2')

    def calc_sp_carbon_count(self, molecule, selection) -> int:
        """Calculate sp hybridized carbon atom count."""
        molecule.assign_hybridization()
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'C' and molecule.atoms[idx].hybridization == 'sp')

    # Hybridization-dependent heteroatom counts
    def calc_sp3_nitrogen_count(self, molecule, selection) -> int:
        molecule.assign_hybridization()
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'N' and molecule.atoms[idx].hybridization == 'sp3')

    def calc_sp2_nitrogen_count(self, molecule, selection) -> int:
        molecule.assign_hybridization()
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'N' and molecule.atoms[idx].hybridization == 'sp2')

    def calc_sp3_oxygen_count(self, molecule, selection) -> int:
        molecule.assign_hybridization()
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'O' and molecule.atoms[idx].hybridization == 'sp3')

    def calc_sp2_oxygen_count(self, molecule, selection) -> int:
        molecule.assign_hybridization()
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'O' and molecule.atoms[idx].hybridization == 'sp2')

    # Ring-related descriptors
    def calc_hetero_ring_count(self, molecule, selection) -> int:
        """Calculate number of rings containing heteroatoms."""
        rings = molecule.find_rings()
        selected_set = self.get_selected_set(selection)
        count = 0
        for ring in rings:
            if all(idx in selected_set for idx in ring):
                if any(molecule.atoms[idx].symbol not in ['C', 'H'] for idx in ring):
                    count += 1
        return count

    def calc_ring5_count(self, molecule, selection) -> int:
        return sum(1 for size in self.get_ring_sizes(molecule, selection) if size == 5)

    def calc_ring6_count(self, molecule, selection) -> int:
        return sum(1 for size in self.get_ring_sizes(molecule, selection) if size == 6)

    def calc_largest_ring_size(self, molecule, selection) -> int:
        sizes = self.get_ring_sizes(molecule, selection)
        return max(sizes) if sizes else 0

    # Functional group counts
    def calc_amide_bond_count(self, molecule, selection) -> int:
        """Calculate amide bond count (-C(=O)-N-)."""
        selected_set = self.get_selected_set(selection)
        count = 0
        for bond in molecule.bonds:
            if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set:
                atom1 = molecule.atoms[bond.begin_atom_idx]
                atom2 = molecule.atoms[bond.end_atom_idx]
                if (atom1.symbol == 'C' and atom2.symbol == 'N') or (atom1.symbol == 'N' and atom2.symbol == 'C'):
                    c_idx = bond.begin_atom_idx if atom1.symbol == 'C' else bond.end_atom_idx
                    for nb_idx in molecule.get_neighbors(c_idx):
                        nb_bond = molecule.get_bond_between(c_idx, nb_idx)
                        if nb_bond and nb_bond.is_double and molecule.atoms[nb_idx].symbol == 'O':
                            count += 1
                            break
        return count

    def calc_ester_bond_count(self, molecule, selection) -> int:
        """Calculate ester group count (R-COO-R', where R' != H)."""
        selected_set = self.get_selected_set(selection)
        count = 0
        for idx in self.get_selected_atoms(molecule, selection):
            atom = molecule.atoms[idx]
            if atom.symbol == 'O':
                neighbors = molecule.get_neighbors(idx)
                if len(neighbors) == 2:
                    # Check if neighbors are two carbons (one carbonyl, one alkyl/aryl)
                    n_atoms = [molecule.atoms[n] for n in neighbors if n < len(molecule.atoms)]
                    if all(a.symbol == 'C' for a in n_atoms):
                        # At least one Carbon must be a carbonyl (C=O)
                        is_ester = False
                        for n_idx in neighbors:
                            for nn_idx in molecule.get_neighbors(n_idx):
                                if nn_idx != idx and nn_idx < len(molecule.atoms) and molecule.atoms[nn_idx].symbol == 'O':
                                    bond = molecule.get_bond_between(n_idx, nn_idx)
                                    if bond and bond.is_double:
                                        is_ester = True
                                        break
                            if is_ester: break
                        if is_ester:
                            count += 1
        return count

    def calc_carbonyl_bond_count(self, molecule, selection) -> int:
        """Calculate carbonyl group count (C=O)."""
        selected_set = self.get_selected_set(selection)
        count = 0
        for bond in molecule.bonds:
            if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set and bond.is_double:
                atom1 = molecule.atoms[bond.begin_atom_idx]
                atom2 = molecule.atoms[bond.end_atom_idx]
                if (atom1.symbol == 'C' and atom2.symbol == 'O') or (atom1.symbol == 'O' and atom2.symbol == 'C'):
                    count += 1
        return count

    def calc_hydroxyl_group_count(self, molecule, selection) -> int:
        """Calculate hydroxyl group (-OH) count.

        Counts implicit hydrogens as well, otherwise every structure coming
        from SMILES would report zero hydroxyls. """
        count = 0
        for idx in self.get_selected_atoms(molecule, selection):
            if molecule.atoms[idx].symbol == 'O' and self.count_hydrogens(molecule, idx) > 0:
                count += 1
        return count

    def calc_carboxyl_group_count(self, molecule, selection) -> int:
        """Calculate carboxyl group (-COOH) count.

        A carboxyl needs a carbonyl oxygen *and* a hydroxyl oxygen on the same
        carbon; without the double bond check a sugar's anomeric centre
        (ring O plus OH) would be counted as an acid. """
        count = 0
        for idx in self.get_selected_atoms(molecule, selection):
            if molecule.atoms[idx].symbol != 'C':
                continue
            carbonyl = False
            hydroxyl = False
            for n in molecule.get_neighbors(idx):
                if n >= len(molecule.atoms) or molecule.atoms[n].symbol != 'O':
                    continue
                bond = molecule.get_bond_between(idx, n)
                if bond is not None and bond.is_double:
                    carbonyl = True
                elif self.count_hydrogens(molecule, n) > 0:
                    hydroxyl = True
            if carbonyl and hydroxyl:
                count += 1
        return count

    def calc_amine_group_count(self, molecule, selection) -> int:
        """Calculate amine group count."""
        count = 0
        for idx in self.get_selected_atoms(molecule, selection):
            if molecule.atoms[idx].symbol == 'N':
                neighbors = molecule.get_neighbors(idx)
                is_amide = False
                for n in neighbors:
                    if n < len(molecule.atoms) and molecule.atoms[n].symbol == 'C':
                        for nn in molecule.get_neighbors(n):
                            if nn < len(molecule.atoms) and molecule.atoms[nn].symbol == 'O':
                                bond = molecule.get_bond_between(n, nn)
                                if bond and bond.is_double:
                                    is_amide = True
                                    break
                if not is_amide:
                    count += 1
        return count

    def calc_methyl_group_count(self, molecule, selection) -> int:
        """Calculate methyl group (-CH3) count, implicit hydrogens included."""
        count = 0
        for idx in self.get_selected_atoms(molecule, selection):
            if molecule.atoms[idx].symbol == 'C' and self.count_hydrogens(molecule, idx) == 3:
                count += 1
        return count

    # Chain and branch descriptors
    def calc_longest_chain(self, molecule, selection) -> int:
        """Calculate longest carbon chain length."""
        selected_set = self.get_selected_set(selection)
        carbon_indices = [idx for idx in selection.atom_indices
                         if idx < len(molecule.atoms) and molecule.atoms[idx].symbol == 'C']
        if not carbon_indices:
            return 0

        max_chain = 0
        visited = set()

        def dfs(idx, depth):
            nonlocal max_chain
            max_chain = max(max_chain, depth)
            visited.add(idx)
            for neighbor in molecule.get_neighbors(idx):
                if neighbor in selected_set and neighbor not in visited:
                    if neighbor < len(molecule.atoms) and molecule.atoms[neighbor].symbol == 'C':
                        dfs(neighbor, depth + 1)
            visited.remove(idx)

        for start in carbon_indices:
            dfs(start, 1)

        return max_chain

    def calc_branch_count(self, molecule, selection) -> int:
        """Calculate number of branch points (carbons with >2 carbon neighbors)."""
        selected_set = self.get_selected_set(selection)
        count = 0
        for idx in self.get_selected_atoms(molecule, selection):
            if molecule.atoms[idx].symbol == 'C':
                neighbors = molecule.get_neighbors(idx)
                carbon_neighbors = sum(1 for n in neighbors
                                      if n in selected_set and n < len(molecule.atoms)
                                      and molecule.atoms[n].symbol == 'C')
                if carbon_neighbors > 2:
                    count += 1
        return count

    def calc_fragment_count(self, molecule, selection) -> int:
        """Calculate number of disconnected fragments."""
        selected_set = self.get_selected_set(selection)
        if not selected_set:
            return 0

        visited = set()
        fragments = 0

        for start in selection.atom_indices:
            if start in visited:
                continue
            fragments += 1
            queue = [start]
            while queue:
                current = queue.pop(0)
                if current in visited:
                    continue
                visited.add(current)
                for neighbor in molecule.get_neighbors(current):
                    if neighbor in selected_set and neighbor not in visited:
                        queue.append(neighbor)

        return fragments

    # Ratios and proportions
    def calc_hetero_atom_count(self, molecule, selection) -> int:
        """Calculate hetero atom count (non-C, non-H)."""
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol not in ['C', 'H'])

    def calc_aliphatic_carbon_count(self, molecule, selection) -> int:
        """Calculate aliphatic carbon atom count."""
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'C' and not molecule.atoms[idx].is_aromatic)

    def calc_aromatic_carbon_count(self, molecule, selection) -> int:
        """Calculate aromatic carbon atom count."""
        return sum(1 for idx in self.get_selected_atoms(molecule, selection)
                  if molecule.atoms[idx].symbol == 'C' and molecule.atoms[idx].is_aromatic)

    def calc_nonring_carbon_count(self, molecule, selection) -> int:
        """Calculate non-ring carbon atoms."""
        rings = molecule.find_rings()
        selected_set = self.get_selected_set(selection)
        ring_carbons = set()
        for ring in rings:
            if all(idx in selected_set for idx in ring):
                for idx in ring:
                    if molecule.atoms[idx].symbol == 'C':
                        ring_carbons.add(idx)
        return len([idx for idx in selection.atom_indices
                   if idx < len(molecule.atoms) and molecule.atoms[idx].symbol == 'C'
                   and idx not in ring_carbons])

    def calc_formal_charge(self, molecule, selection) -> int:
        """Calculate total formal charge."""
        return sum(getattr(molecule.atoms[idx], 'formal_charge', 0)
                  for idx in self.get_selected_atoms(molecule, selection))

    def calc_unsaturation_count(self, molecule, selection) -> float:
        """Calculate degree of unsaturation (double bond equivalents)."""
        c_count = self.calc_carbon_count(molecule, selection)
        h_count = self.calc_hydrogen_count(molecule, selection)
        n_count = self.calc_nitrogen_count(molecule, selection)
        x_count = self.calc_halogen_count(molecule, selection)
        return c_count - (h_count + x_count)/2.0 + n_count/2.0 + 1.0

    def calc_aromatic_proportion(self, molecule, selection) -> float:
        """Calculate proportion of aromatic atoms (relative to heavy atoms)."""
        heavy_count = self.calc_heavy_atom_count(molecule, selection)
        if heavy_count == 0:
            return 0.0
        aromatic_count = sum(1 for idx in self.get_selected_atoms(molecule, selection)
                            if molecule.atoms[idx].symbol != 'H' and molecule.atoms[idx].is_aromatic)
        return aromatic_count / heavy_count

    def calc_aliphatic_proportion(self, molecule, selection) -> float:
        """Calculate proportion of aliphatic atoms (relative to heavy atoms)."""
        heavy_count = self.calc_heavy_atom_count(molecule, selection)
        if heavy_count == 0:
            return 0.0
        aliphatic_count = sum(1 for idx in self.get_selected_atoms(molecule, selection)
                             if molecule.atoms[idx].symbol != 'H' and not molecule.atoms[idx].is_aromatic)
        return aliphatic_count / heavy_count

    def calc_sp3_proportion(self, molecule, selection) -> float:
        """Calculate proportion of sp3 hybridized carbons."""
        total_c = self.calc_carbon_count(molecule, selection)
        if total_c == 0:
            return 0.0
        return self.calc_sp3_carbon_count(molecule, selection) / total_c

    def calc_hetero_proportion(self, molecule, selection) -> float:
        """Calculate proportion of hetero atoms (relative to heavy atoms)."""
        heavy_count = self.calc_heavy_atom_count(molecule, selection)
        if heavy_count == 0:
            return 0.0
        return self.calc_hetero_atom_count(molecule, selection) / heavy_count

    def calc_hc_ratio(self, molecule, selection) -> float:
        """Calculate hydrogen to carbon ratio."""
        c_count = self.calc_carbon_count(molecule, selection)
        if c_count == 0:
            return 0.0
        return self.calc_hydrogen_count(molecule, selection) / c_count

    def calc_nc_ratio(self, molecule, selection) -> float:
        """Calculate nitrogen to carbon ratio."""
        c_count = self.calc_carbon_count(molecule, selection)
        if c_count == 0:
            return 0.0
        return self.calc_nitrogen_count(molecule, selection) / c_count

    def calc_oc_ratio(self, molecule, selection) -> float:
        """Calculate oxygen to carbon ratio."""
        c_count = self.calc_carbon_count(molecule, selection)
        if c_count == 0:
            return 0.0
        return self.calc_oxygen_count(molecule, selection) / c_count
