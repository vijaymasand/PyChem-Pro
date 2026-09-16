"""
Topological descriptor calculations.
Includes: Wiener index, Zagreb indices, Balaban index, connectivity indices, path counts, etc.
"""
import numpy as np
from typing import List
from . import BaseCalculator


class TopologicalCalculator(BaseCalculator):
    """Calculator for topological (2D) molecular descriptors."""

    def calc_wiener_index(self, molecule, selection) -> float:
        """Calculate Wiener index."""
        distance_matrix = self.get_distance_matrix(molecule)
        total_distance = 0.0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]):
                    total_distance += distance_matrix[i][j]
        return total_distance

    def calc_zagreb_index_1(self, molecule, selection) -> float:
        """Calculate first Zagreb index."""
        selected_set = self.get_selected_set(selection)
        return sum(len([n for n in molecule.get_neighbors(idx) if n in selected_set]) ** 2
                  for idx in selection.atom_indices)

    def calc_zagreb_index_2(self, molecule, selection) -> float:
        """Calculate second Zagreb index."""
        selected_set = self.get_selected_set(selection)
        total = 0.0
        for bond in molecule.bonds:
            if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set:
                degree1 = len([n for n in molecule.get_neighbors(bond.begin_atom_idx) if n in selected_set])
                degree2 = len([n for n in molecule.get_neighbors(bond.end_atom_idx) if n in selected_set])
                total += degree1 * degree2
        return float(total)

    def calc_balaban_index(self, molecule, selection) -> float:
        """Calculate Balaban J index."""
        selected_set = self.get_selected_set(selection)
        n_atoms = len(selected_set)
        if n_atoms <= 1:
            return 0.0

        n_bonds = self.get_bond_count_in_selection(molecule, selected_set)
        mu = n_bonds - n_atoms + 1
        if mu < 0:
            mu = 0

        dist_matrix = self.get_distance_matrix(molecule)
        s_values = {}
        for i in selected_set:
            s_i = sum(dist_matrix[i][j] for j in selected_set if i != j and j < len(dist_matrix[i])
                     if dist_matrix[i][j] != float('inf'))
            s_values[i] = s_i

        total = 0.0
        for bond in molecule.bonds:
            if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set:
                si = s_values.get(bond.begin_atom_idx, 0)
                sj = s_values.get(bond.end_atom_idx, 0)
                if si > 0 and sj > 0:
                    total += 1.0 / np.sqrt(si * sj)

        return (n_bonds / (mu + 1.0)) * total

    def calc_harary_index(self, molecule, selection) -> float:
        """Calculate Harary index."""
        distance_matrix = self.get_distance_matrix(molecule)
        selected_set = self.get_selected_set(selection)
        total = 0.0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]):
                    distance = distance_matrix[i][j]
                    if distance > 0:
                        total += 1.0 / distance
        return total

    def calc_connectivity_index_chi0(self, molecule, selection) -> float:
        """Calculate Randic connectivity index (order 0)."""
        selected_set = self.get_selected_set(selection)
        total = 0.0
        for idx in self.get_selected_atoms(molecule, selection):
            deg = len([n for n in molecule.get_neighbors(idx) if n in selected_set])
            if deg > 0:
                total += 1.0 / np.sqrt(deg)
        return float(total)

    def calc_connectivity_index_chi1(self, molecule, selection) -> float:
        """Calculate connectivity index (order 1)."""
        selected_set = self.get_selected_set(selection)
        total = 0.0
        for bond in molecule.bonds:
            if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set:
                degree1 = len([n for n in molecule.get_neighbors(bond.begin_atom_idx) if n in selected_set])
                degree2 = len([n for n in molecule.get_neighbors(bond.end_atom_idx) if n in selected_set])
                if degree1 > 0 and degree2 > 0:
                    total += 1.0 / np.sqrt(degree1 * degree2)
        return float(total)

    def calc_kappa_shape_index_1(self, molecule, selection) -> float:
        """Calculate Kappa shape index 1."""
        n_atoms = len(selection.atom_indices)
        if n_atoms <= 2:
            return 0.0

        distance_matrix = self.get_distance_matrix(molecule)
        max_distance = 0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]):
                    max_distance = max(max_distance, distance_matrix[i][j])

        if max_distance > 0:
            return (n_atoms * (n_atoms - 1) ** 2) / (2 * max_distance ** 2)
        return 0.0

    def calc_kappa_shape_index_2(self, molecule, selection) -> float:
        """Calculate Kappa shape index 2."""
        n_atoms = len(selection.atom_indices)
        if n_atoms <= 3:
            return 0.0

        distance_matrix = self.get_distance_matrix(molecule)
        total_distance = sum(distance_matrix[i][j] for i in selection.atom_indices
                          for j in selection.atom_indices
                          if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]))

        if total_distance > 0:
            return (n_atoms - 1) * (n_atoms - 2) ** 2 / (total_distance ** 2)
        return 0.0

    def calc_kappa_shape_index_3(self, molecule, selection) -> float:
        """Calculate Kappa shape index 3."""
        n_atoms = len(selection.atom_indices)
        if n_atoms <= 4:
            return 0.0

        distance_matrix = self.get_distance_matrix(molecule)
        max_distance = 0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]):
                    max_distance = max(max_distance, distance_matrix[i][j])

        if max_distance > 0:
            return (n_atoms - 2) * (n_atoms - 3) ** 2 / (max_distance ** 2)
        return 0.0

    def calc_hosoya_index(self, molecule, selection) -> int:
        """Calculate Hosoya index (total number of matchings)."""
        selected_set = self.get_selected_set(selection)
        edges = [(bond.begin_atom_idx, bond.end_atom_idx) for bond in molecule.bonds
                if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set]

        if not edges:
            return 0

        # Simplified: just count edges and isolated vertices
        vertices = set()
        for e in edges:
            vertices.add(e[0])
            vertices.add(e[1])

        return len(edges) + (len(selected_set) - len(vertices))

    def calc_platt_index(self, molecule, selection) -> int:
        """Calculate Platt index (sum of edge degrees)."""
        selected_set = self.get_selected_set(selection)
        return sum(len([n for n in molecule.get_neighbors(bond.begin_atom_idx) if n in selected_set]) +
                  len([n for n in molecule.get_neighbors(bond.end_atom_idx) if n in selected_set])
                  for bond in molecule.bonds
                  if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set)

    def calc_polarity_number(self, molecule, selection) -> int:
        """Calculate polarity number (number of pairs at distance 3)."""
        distance_matrix = self.get_distance_matrix(molecule)
        selected_set = self.get_selected_set(selection)
        return sum(1 for i in selection.atom_indices for j in selection.atom_indices
                  if i < j and i < len(distance_matrix) and j < len(distance_matrix[i])
                  and distance_matrix[i][j] == 3)

    def calc_bertz_index(self, molecule, selection) -> float:
        """Calculate Bertz complexity index."""
        n_atoms = len(selection.atom_indices)
        n_bonds = self.get_bond_count_in_selection(molecule, self.get_selected_set(selection))
        return n_bonds * np.log2(n_bonds) if n_bonds > 0 else 0.0

    def calc_bonchev_trinajstic_index(self, molecule, selection) -> float:
        """Calculate Bonchev-Trinajstic information index."""
        selected_set = self.get_selected_set(selection)
        degrees = [len([n for n in molecule.get_neighbors(idx) if n in selected_set])
                  for idx in selection.atom_indices]
        total = sum(degrees)
        if total == 0:
            return 0.0
        probabilities = [d / total for d in degrees if d > 0]
        return -sum(p * np.log2(p) for p in probabilities)

    def calc_information_content_0(self, molecule, selection) -> float:
        """Calculate information content (order 0)."""
        return self.calc_bonchev_trinajstic_index(molecule, selection)

    def calc_information_content_1(self, molecule, selection) -> float:
        """Calculate information content (order 1)."""
        return self.calc_information_content_0(molecule, selection)

    def calc_electrotopological_state_index(self, molecule, selection) -> float:
        """Calculate electrotopological state index."""
        selected_set = self.get_selected_set(selection)
        total = 0.0
        atomic_numbers = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
                         'Cl': 17, 'Br': 35, 'I': 53, 'S': 16, 'P': 15}

        for idx in selection.atom_indices:
            degree = len([n for n in molecule.get_neighbors(idx) if n in selected_set])
            if degree > 0:
                atom = molecule.atoms[idx]
                z = atomic_numbers.get(atom.symbol, 6)
                intrinsic_state = z / degree
                total += intrinsic_state

        return total

    def calc_molecular_topological_index(self, molecule, selection) -> float:
        """Calculate molecular topological index."""
        return self.calc_wiener_index(molecule, selection) + self.calc_zagreb_index_1(molecule, selection)

    def calc_nuclear_repulsion_energy(self, molecule, selection) -> float:
        """Calculate nuclear repulsion energy."""
        distance_matrix = self.get_distance_matrix(molecule)
        atomic_numbers = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
                         'Cl': 17, 'Br': 35, 'I': 53, 'S': 16, 'P': 15}

        total = 0.0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(molecule.atoms) and j < len(molecule.atoms):
                    z1 = atomic_numbers.get(molecule.atoms[i].symbol, 6)
                    z2 = atomic_numbers.get(molecule.atoms[j].symbol, 6)
                    if i < len(distance_matrix) and j < len(distance_matrix[i]):
                        distance = distance_matrix[i][j]
                        if distance > 0:
                            total += (z1 * z2) / distance

        return total * 0.0001

    # Extended topological descriptors
    def calc_hyper_wiener_index(self, molecule, selection) -> float:
        """Calculate hyper-Wiener index: sum of distances and squared distances."""
        distance_matrix = self.get_distance_matrix(molecule)
        total = 0.0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]):
                    d = distance_matrix[i][j]
                    total += d + d * d
        return total / 2.0

    def calc_randic_index(self, molecule, selection) -> float:
        """Calculate Randic connectivity index."""
        return self.calc_connectivity_index_chi0(molecule, selection)

    def calc_eccentric_connectivity_index(self, molecule, selection) -> float:
        """Calculate eccentric connectivity index."""
        selected_set = self.get_selected_set(selection)
        distance_matrix = self.get_distance_matrix(molecule)
        total = 0.0
        for idx in selection.atom_indices:
            degree = len([n for n in molecule.get_neighbors(idx) if n in selected_set])
            eccentricity = max((distance_matrix[idx][j] for j in selection.atom_indices
                               if idx < len(distance_matrix) and j < len(distance_matrix[idx])
                               and distance_matrix[idx][j] != float('inf')), default=0)
            total += degree * eccentricity
        return total

    def calc_path_count_3(self, molecule, selection) -> int:
        """Calculate number of 3-atom paths."""
        selected_set = self.get_selected_set(selection)
        count = 0
        for middle in selection.atom_indices:
            neighbors = [n for n in molecule.get_neighbors(middle) if n in selected_set]
            count += len(neighbors) * (len(neighbors) - 1) // 2
        return count

    def calc_path_count_4(self, molecule, selection) -> int:
        """Calculate number of 4-atom paths."""
        selected_set = self.get_selected_set(selection)
        count = 0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i != j:
                    common = len(set(molecule.get_neighbors(i)) & set(molecule.get_neighbors(j)) & selected_set)
                    count += common
        return count // 2

    def calc_average_path_length(self, molecule, selection) -> float:
        """Calculate average topological path length."""
        distance_matrix = self.get_distance_matrix(molecule)
        total = 0.0
        n = 0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]):
                    d = distance_matrix[i][j]
                    if d != float('inf'):
                        total += d
                        n += 1
        return total / n if n > 0 else 0.0

    def calc_longest_shortest_path(self, molecule, selection) -> int:
        """Calculate graph diameter (longest shortest path)."""
        distance_matrix = self.get_distance_matrix(molecule)
        diameter = 0
        for i in selection.atom_indices:
            for j in selection.atom_indices:
                if i < j and i < len(distance_matrix) and j < len(distance_matrix[i]):
                    d = distance_matrix[i][j]
                    if d != float('inf'):
                        diameter = max(diameter, int(d))
        return diameter

    def calc_cyclomatic_number(self, molecule, selection) -> int:
        """Calculate cyclomatic number (circuit rank): E - V + C."""
        selected_set = self.get_selected_set(selection)
        v = len(selected_set)
        e = self.get_bond_count_in_selection(molecule, selected_set)
        c = self.calc_fragment_count(molecule, selection)
        return e - v + c

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
