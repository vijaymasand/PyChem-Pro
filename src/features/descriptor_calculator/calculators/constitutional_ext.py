"""
Extended constitutional (1D) descriptors.

Functional group counts, ring system topology (spiro/bridgehead/fused),
carbon substitution classes and mass related quantities. All of them work on
the molecular graph alone, so they are available even without 3D coordinates.
"""
from .constitutional import ConstitutionalCalculator
from .graph_utils import get_subgraph, VDW_VOLUME

HALOGENS = {'F', 'Cl', 'Br', 'I'}


class ExtendedConstitutionalCalculator(ConstitutionalCalculator):
    """Constitutional descriptors, basic set plus the extended families."""

    # ------------------------------------------------------------ helpers
    def _selected(self, molecule, selection):
        return set(self.get_selected_atoms(molecule, selection))

    def _heavy_neighbors(self, molecule, idx, chosen):
        return [n for n in molecule.get_neighbors(idx)
                if n in chosen and molecule.atoms[n].symbol != 'H']

    def _bond_between(self, molecule, i, j):
        return molecule.get_bond_between(i, j)

    def _double_bonded_to(self, molecule, idx, symbol, chosen):
        """Neighbours of idx connected by a double bond and matching symbol."""
        found = []
        for n in molecule.get_neighbors(idx):
            if n not in chosen:
                continue
            if molecule.atoms[n].symbol != symbol:
                continue
            bond = self._bond_between(molecule, idx, n)
            if bond is not None and bond.is_double:
                found.append(n)
        return found

    def _total_h(self, molecule, idx, chosen):
        atom = molecule.atoms[idx]
        explicit = sum(1 for n in molecule.get_neighbors(idx)
                       if molecule.atoms[n].symbol == 'H')
        return getattr(atom, 'total_h', 0) + explicit

    # -------------------------------------------------------- mass / size
    def calc_exact_mass(self, molecule, selection) -> float:
        """Monoisotopic-style mass from the element masses, implicit H included."""
        chosen = self._selected(molecule, selection)
        total = 0.0
        for idx in chosen:
            atom = molecule.atoms[idx]
            total += getattr(atom, 'mass', 0.0)
            total += self._total_h(molecule, idx, chosen) * 1.008 if atom.symbol != 'H' else 0.0
        return total

    def calc_heavy_atom_mol_weight(self, molecule, selection) -> float:
        """Molecular weight counting heavy atoms only."""
        chosen = self._selected(molecule, selection)
        return sum(getattr(molecule.atoms[idx], 'mass', 0.0)
                   for idx in chosen if molecule.atoms[idx].symbol != 'H')

    def calc_average_atomic_mass(self, molecule, selection) -> float:
        """Mean mass per heavy atom, a crude measure of how heavy the elements are."""
        chosen = [i for i in self._selected(molecule, selection)
                  if molecule.atoms[i].symbol != 'H']
        if not chosen:
            return 0.0
        return sum(getattr(molecule.atoms[i], 'mass', 0.0) for i in chosen) / len(chosen)

    def calc_valence_electron_count(self, molecule, selection) -> int:
        """Total number of valence electrons, correlates with polarizability."""
        from .graph_utils import valence_electrons
        chosen = self._selected(molecule, selection)
        total = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            total += valence_electrons(atom.symbol)
            if atom.symbol != 'H':
                total += self._total_h(molecule, idx, chosen)
        return total

    def calc_vdw_volume_sum(self, molecule, selection) -> float:
        """Sum of atomic van der Waals volumes (Bondi), no packing correction."""
        chosen = self._selected(molecule, selection)
        return sum(VDW_VOLUME.get(molecule.atoms[i].symbol, 20.58) for i in chosen)

    # --------------------------------------------------- carbon classes
    def _carbon_substitution(self, molecule, selection, degree):
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.symbol != 'C':
                continue
            carbons = sum(1 for n in self._heavy_neighbors(molecule, idx, chosen)
                          if molecule.atoms[n].symbol == 'C')
            if carbons == degree:
                count += 1
        return count

    def calc_primary_carbon_count(self, molecule, selection) -> int:
        return self._carbon_substitution(molecule, selection, 1)

    def calc_secondary_carbon_count(self, molecule, selection) -> int:
        return self._carbon_substitution(molecule, selection, 2)

    def calc_tertiary_carbon_count(self, molecule, selection) -> int:
        return self._carbon_substitution(molecule, selection, 3)

    def calc_quaternary_carbon_count(self, molecule, selection) -> int:
        return self._carbon_substitution(molecule, selection, 4)

    def calc_terminal_atom_count(self, molecule, selection) -> int:
        """Heavy atoms with a single heavy neighbour."""
        chosen = self._selected(molecule, selection)
        return sum(1 for idx in chosen
                   if molecule.atoms[idx].symbol != 'H'
                   and len(self._heavy_neighbors(molecule, idx, chosen)) == 1)

    def calc_max_atom_degree(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        degrees = [len(self._heavy_neighbors(molecule, i, chosen)) for i in chosen]
        return max(degrees) if degrees else 0

    def calc_mean_atom_degree(self, molecule, selection) -> float:
        chosen = [i for i in self._selected(molecule, selection)
                  if molecule.atoms[i].symbol != 'H']
        if not chosen:
            return 0.0
        selected_set = self._selected(molecule, selection)
        return sum(len(self._heavy_neighbors(molecule, i, selected_set))
                   for i in chosen) / len(chosen)

    # -------------------------------------------------------- ring system
    def _selected_rings(self, molecule, selection):
        chosen = self._selected(molecule, selection)
        return [r for r in molecule.find_rings() if all(i in chosen for i in r)]

    def calc_ring_atom_count(self, molecule, selection) -> int:
        atoms = set()
        for ring in self._selected_rings(molecule, selection):
            atoms.update(ring)
        return len(atoms)

    def calc_ring_bond_count(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        return sum(1 for bond in molecule.bonds
                   if bond.is_in_ring
                   and bond.begin_atom_idx in chosen and bond.end_atom_idx in chosen)

    def calc_aliphatic_ring_count(self, molecule, selection) -> int:
        """Rings with at least one non-aromatic atom."""
        rings = self._selected_rings(molecule, selection)
        return sum(1 for r in rings
                   if not all(molecule.atoms[i].is_aromatic for i in r))

    def calc_saturated_ring_count(self, molecule, selection) -> int:
        """Rings whose bonds are all single."""
        count = 0
        for ring in self._selected_rings(molecule, selection):
            saturated = True
            size = len(ring)
            for k in range(size):
                bond = self._bond_between(molecule, ring[k], ring[(k + 1) % size])
                if bond is None or not bond.is_single or bond.is_aromatic:
                    saturated = False
                    break
            if saturated:
                count += 1
        return count

    def calc_aromatic_heterocycle_count(self, molecule, selection) -> int:
        rings = self._selected_rings(molecule, selection)
        return sum(1 for r in rings
                   if all(molecule.atoms[i].is_aromatic for i in r)
                   and any(molecule.atoms[i].symbol != 'C' for i in r))

    def calc_aromatic_carbocycle_count(self, molecule, selection) -> int:
        rings = self._selected_rings(molecule, selection)
        return sum(1 for r in rings
                   if all(molecule.atoms[i].is_aromatic for i in r)
                   and all(molecule.atoms[i].symbol == 'C' for i in r))

    def calc_saturated_heterocycle_count(self, molecule, selection) -> int:
        count = 0
        for ring in self._selected_rings(molecule, selection):
            if not any(molecule.atoms[i].symbol != 'C' for i in ring):
                continue
            size = len(ring)
            saturated = True
            for k in range(size):
                bond = self._bond_between(molecule, ring[k], ring[(k + 1) % size])
                if bond is None or not bond.is_single or bond.is_aromatic:
                    saturated = False
                    break
            if saturated:
                count += 1
        return count

    def calc_spiro_atom_count(self, molecule, selection) -> int:
        """Atoms shared by two rings that have no bond in common."""
        rings = [set(r) for r in self._selected_rings(molecule, selection)]
        spiro = set()
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                shared = rings[i] & rings[j]
                if len(shared) == 1:
                    spiro.update(shared)
        return len(spiro)

    def calc_bridgehead_atom_count(self, molecule, selection) -> int:
        """Atoms shared by two rings that also share at least one bond."""
        rings = [set(r) for r in self._selected_rings(molecule, selection)]
        bridgeheads = set()
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                shared = rings[i] & rings[j]
                if len(shared) >= 2:
                    for idx in shared:
                        neighbors = set(molecule.get_neighbors(idx)) & shared
                        # a bridgehead sits at the end of the shared path
                        if len(neighbors) <= 1 or len(shared) > 2:
                            bridgeheads.add(idx)
        return len(bridgeheads)

    def calc_fused_ring_count(self, molecule, selection) -> int:
        """Number of ring pairs sharing a bond."""
        rings = [set(r) for r in self._selected_rings(molecule, selection)]
        fused = 0
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                if len(rings[i] & rings[j]) >= 2:
                    fused += 1
        return fused

    def calc_macrocycle_count(self, molecule, selection) -> int:
        """Rings with more than 12 atoms."""
        return sum(1 for r in self._selected_rings(molecule, selection) if len(r) > 12)

    def calc_ring_complexity(self, molecule, selection) -> float:
        """Ring atoms per ring, a simple measure of ring system condensation."""
        rings = self._selected_rings(molecule, selection)
        if not rings:
            return 0.0
        atoms = set()
        for ring in rings:
            atoms.update(ring)
        return len(atoms) / len(rings)

    # ---------------------------------------------------- functional groups
    def calc_nitro_group_count(self, molecule, selection) -> int:
        """R-NO2, both the charge separated and the pentavalent notation."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'N':
                continue
            oxygens = [n for n in molecule.get_neighbors(idx)
                       if n in chosen and molecule.atoms[n].symbol == 'O'
                       and len(self._heavy_neighbors(molecule, n, chosen)) == 1]
            if len(oxygens) >= 2:
                count += 1
        return count

    def calc_nitrile_count(self, molecule, selection) -> int:
        """C#N triple bonds."""
        chosen = self._selected(molecule, selection)
        count = 0
        for bond in molecule.bonds:
            if not bond.is_triple:
                continue
            i, j = bond.begin_atom_idx, bond.end_atom_idx
            if i not in chosen or j not in chosen:
                continue
            symbols = {molecule.atoms[i].symbol, molecule.atoms[j].symbol}
            if symbols == {'C', 'N'}:
                count += 1
        return count

    def calc_aldehyde_count(self, molecule, selection) -> int:
        """C(=O)H, a carbonyl carbon with a hydrogen and one carbon neighbour."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'C':
                continue
            if not self._double_bonded_to(molecule, idx, 'O', chosen):
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            carbons = [n for n in heavy if molecule.atoms[n].symbol == 'C']
            oxygens = [n for n in heavy if molecule.atoms[n].symbol == 'O']
            if len(oxygens) == 1 and len(carbons) <= 1 and self._total_h(molecule, idx, chosen) >= 1:
                count += 1
        return count

    def calc_ketone_count(self, molecule, selection) -> int:
        """C(=O) flanked by two carbons and no other heteroatom."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'C':
                continue
            if not self._double_bonded_to(molecule, idx, 'O', chosen):
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            carbons = [n for n in heavy if molecule.atoms[n].symbol == 'C']
            oxygens = [n for n in heavy if molecule.atoms[n].symbol == 'O']
            others = [n for n in heavy if molecule.atoms[n].symbol not in ('C', 'O')]
            if len(carbons) == 2 and len(oxygens) == 1 and not others:
                count += 1
        return count

    def calc_ether_count(self, molecule, selection) -> int:
        """C-O-C oxygens that are not part of an ester or acid."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.symbol != 'O' or atom.is_aromatic:
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            if len(heavy) != 2 or self._total_h(molecule, idx, chosen):
                continue
            if all(molecule.atoms[n].symbol == 'C' for n in heavy):
                carbonyl = any(self._double_bonded_to(molecule, n, 'O', chosen) for n in heavy)
                if not carbonyl:
                    count += 1
        return count

    def calc_thioether_count(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.symbol != 'S' or atom.is_aromatic:
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            if len(heavy) == 2 and not self._total_h(molecule, idx, chosen):
                if all(molecule.atoms[n].symbol == 'C' for n in heavy):
                    count += 1
        return count

    def calc_thiol_count(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        return sum(1 for idx in chosen
                   if molecule.atoms[idx].symbol == 'S'
                   and len(self._heavy_neighbors(molecule, idx, chosen)) == 1
                   and self._total_h(molecule, idx, chosen) >= 1)

    def calc_sulfonamide_count(self, molecule, selection) -> int:
        """S(=O)(=O)N."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'S':
                continue
            if len(self._double_bonded_to(molecule, idx, 'O', chosen)) >= 2:
                if any(molecule.atoms[n].symbol == 'N'
                       for n in self._heavy_neighbors(molecule, idx, chosen)):
                    count += 1
        return count

    def calc_sulfone_count(self, molecule, selection) -> int:
        """C-S(=O)(=O)-C."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'S':
                continue
            if len(self._double_bonded_to(molecule, idx, 'O', chosen)) >= 2:
                carbons = [n for n in self._heavy_neighbors(molecule, idx, chosen)
                           if molecule.atoms[n].symbol == 'C']
                if len(carbons) >= 2:
                    count += 1
        return count

    def calc_sulfoxide_count(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'S':
                continue
            if len(self._double_bonded_to(molecule, idx, 'O', chosen)) == 1:
                carbons = [n for n in self._heavy_neighbors(molecule, idx, chosen)
                           if molecule.atoms[n].symbol == 'C']
                if len(carbons) >= 2:
                    count += 1
        return count

    def calc_phosphate_count(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'P':
                continue
            oxygens = [n for n in self._heavy_neighbors(molecule, idx, chosen)
                       if molecule.atoms[n].symbol == 'O']
            if len(oxygens) >= 3:
                count += 1
        return count

    def calc_imine_count(self, molecule, selection) -> int:
        """C=N double bonds outside aromatic rings."""
        chosen = self._selected(molecule, selection)
        count = 0
        for bond in molecule.bonds:
            if not bond.is_double or bond.is_aromatic:
                continue
            i, j = bond.begin_atom_idx, bond.end_atom_idx
            if i not in chosen or j not in chosen:
                continue
            if {molecule.atoms[i].symbol, molecule.atoms[j].symbol} == {'C', 'N'}:
                count += 1
        return count

    def calc_azo_count(self, molecule, selection) -> int:
        """N=N double bonds."""
        chosen = self._selected(molecule, selection)
        count = 0
        for bond in molecule.bonds:
            if not bond.is_double or bond.is_aromatic:
                continue
            i, j = bond.begin_atom_idx, bond.end_atom_idx
            if i in chosen and j in chosen:
                if molecule.atoms[i].symbol == 'N' and molecule.atoms[j].symbol == 'N':
                    count += 1
        return count

    def calc_urea_count(self, molecule, selection) -> int:
        """N-C(=O)-N."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'C':
                continue
            if not self._double_bonded_to(molecule, idx, 'O', chosen):
                continue
            nitrogens = [n for n in self._heavy_neighbors(molecule, idx, chosen)
                         if molecule.atoms[n].symbol == 'N']
            if len(nitrogens) >= 2:
                count += 1
        return count

    def calc_guanidine_count(self, molecule, selection) -> int:
        """C bonded to three nitrogens, one of them doubly."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'C':
                continue
            nitrogens = [n for n in self._heavy_neighbors(molecule, idx, chosen)
                         if molecule.atoms[n].symbol == 'N']
            if len(nitrogens) >= 3 and self._double_bonded_to(molecule, idx, 'N', chosen):
                count += 1
        return count

    def calc_phenol_count(self, molecule, selection) -> int:
        """OH attached to an aromatic carbon."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.symbol != 'O' or self._total_h(molecule, idx, chosen) < 1:
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            if len(heavy) == 1 and molecule.atoms[heavy[0]].is_aromatic:
                count += 1
        return count

    def calc_aliphatic_hydroxyl_count(self, molecule, selection) -> int:
        """Aliphatic alcohols: hydroxyls that are neither phenolic nor acidic."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.symbol != 'O' or self._total_h(molecule, idx, chosen) < 1:
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            if len(heavy) != 1:
                continue
            neighbor = heavy[0]
            if molecule.atoms[neighbor].is_aromatic:
                continue
            if self._double_bonded_to(molecule, neighbor, 'O', chosen):
                continue  # carboxylic acid
            count += 1
        return count

    def calc_halide_on_aromatic_count(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol not in HALOGENS:
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            if heavy and molecule.atoms[heavy[0]].is_aromatic:
                count += 1
        return count

    def calc_acidic_group_count(self, molecule, selection) -> int:
        """Carboxylic acids, sulfonic acids, tetrazoles and phosphonic acids."""
        return (self.calc_carboxyl_group_count(molecule, selection) +
                self._sulfonic_acid_count(molecule, selection) +
                self.calc_phosphate_count(molecule, selection))

    def _sulfonic_acid_count(self, molecule, selection) -> int:
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            if molecule.atoms[idx].symbol != 'S':
                continue
            if len(self._double_bonded_to(molecule, idx, 'O', chosen)) >= 2:
                hydroxyls = [n for n in self._heavy_neighbors(molecule, idx, chosen)
                             if molecule.atoms[n].symbol == 'O'
                             and self._total_h(molecule, n, chosen) >= 1]
                if hydroxyls:
                    count += 1
        return count

    def calc_basic_group_count(self, molecule, selection) -> int:
        """Basic nitrogens: aliphatic amines, amidines and guanidines."""
        chosen = self._selected(molecule, selection)
        count = 0
        for idx in chosen:
            atom = molecule.atoms[idx]
            if atom.symbol != 'N' or atom.is_aromatic:
                continue
            heavy = self._heavy_neighbors(molecule, idx, chosen)
            # an amide nitrogen is not basic
            amide = any(self._double_bonded_to(molecule, n, 'O', chosen)
                        for n in heavy if molecule.atoms[n].symbol == 'C')
            if not amide and len(heavy) <= 3:
                count += 1
        return count

    def calc_rotatable_bond_fraction(self, molecule, selection) -> float:
        """Rotatable bonds per heavy atom, a size independent flexibility measure."""
        graph = get_subgraph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return self.calc_rotatable_bonds(molecule, selection) / graph.n

    def calc_heteroatom_ratio_heavy(self, molecule, selection) -> float:
        """Heteroatoms per heavy atom."""
        graph = get_subgraph(molecule, selection)
        if graph.n == 0:
            return 0.0
        hetero = sum(1 for a in graph.atoms if a.symbol not in ('C', 'H'))
        return hetero / graph.n
