"""
Enhanced Molecular Descriptor Engine for PyChem.
Calculates 300+ molecular descriptors across 8 categories.
Uses modular calculator classes for better organization and performance.

Each category has a base calculator with the classical descriptors and an
extended calculator that adds the families QSAR work normally expects
(connectivity and shape indices, graph spectra, autocorrelations, CPSA,
Huckel frontier orbitals, drug-likeness filters). Only numpy is used.
"""
import hashlib
import numpy as np
from typing import Dict, List, Optional, Callable, Any
from enum import Enum, auto

from ...core.domain.models.molecule import Molecule
from .descriptor_types import (
    DescriptorCategory,
    DescriptorInfo,
    DescriptorResult,
    CalculationProgress,
    AtomSelection,
    SelectionType
)

# Import calculators. The extended classes subclass the basic ones, so a single
# instance per category serves both the classical and the extended descriptors.
from .calculators import BaseCalculator
from .calculators.fingerprints import FingerprintCalculator
from .calculators.constitutional_ext import ExtendedConstitutionalCalculator
from .calculators.topological_ext import ExtendedTopologicalCalculator
from .calculators.geometric_ext import ExtendedGeometricCalculator
from .calculators.electronic_ext import ExtendedElectronicCalculator
from .calculators.quantum_ext import ExtendedQuantumCalculator
from .calculators.hybrid_ext import ExtendedHybridCalculator


class DescriptorCache:
    """Cache for descriptor calculation results."""

    def __init__(self, max_size: int = 100):
        self.cache = {}
        self.max_size = max_size

    def _get_molecule_hash(self, molecule: Molecule) -> str:
        """Generate hash for molecule state."""
        atoms_str = "".join(f"{a.symbol}{a.x}{a.y}{a.z}" for a in molecule.atoms)
        bonds_str = "".join(f"{b.begin_atom_idx}{b.end_atom_idx}{b.order}" for b in molecule.bonds)
        return hashlib.md5(f"{atoms_str}{bonds_str}".encode()).hexdigest()

    def _get_selection_hash(self, selection: AtomSelection) -> str:
        """Generate hash for atom selection."""
        return hashlib.md5(str(sorted(selection.atom_indices)).encode()).hexdigest()

    def get_cache_key(self, molecule: Molecule, selection: AtomSelection, descriptor_name: str) -> str:
        """Generate unique cache key."""
        mol_hash = self._get_molecule_hash(molecule)
        sel_hash = self._get_selection_hash(selection)
        return f"{mol_hash}:{sel_hash}:{descriptor_name}"

    def get(self, key: str) -> Optional[Any]:
        """Get cached value."""
        return self.cache.get(key)

    def set(self, key: str, value: Any):
        """Set cached value with LRU eviction."""
        if len(self.cache) >= self.max_size:
            oldest_key = next(iter(self.cache))
            del self.cache[oldest_key]
        self.cache[key] = value

    def clear(self):
        """Clear cache."""
        self.cache.clear()


class DescriptorEngine:
    """
    Enhanced molecular descriptor calculation engine.
    Calculates 90+ descriptors across 8 categories.
    """

    def __init__(self, enable_cache: bool = True):
        self.enable_cache = enable_cache
        self.cache = DescriptorCache() if enable_cache else None
        self.progress_callback = None  # For GUI progress updates

        # Initialize calculator instances
        self.constitutional_calc = ExtendedConstitutionalCalculator()
        self.topological_calc = ExtendedTopologicalCalculator()
        self.geometric_calc = ExtendedGeometricCalculator()
        self.electronic_calc = ExtendedElectronicCalculator()
        self.quantum_calc = ExtendedQuantumCalculator()
        self.fingerprint_calc = FingerprintCalculator()
        self.hybrid_calc = ExtendedHybridCalculator()

        # Build descriptor registry
        self.descriptors = self._initialize_descriptors()

    def _initialize_descriptors(self) -> Dict[DescriptorCategory, List[DescriptorInfo]]:
        """Initialize all 90+ descriptors organized by category."""
        return {
            DescriptorCategory.CONSTITUTIONAL: self._get_constitutional_descriptors(),
            DescriptorCategory.TOPOLOGICAL: self._get_topological_descriptors(),
            DescriptorCategory.GEOMETRIC: self._get_geometric_descriptors(),
            DescriptorCategory.ELECTRONIC: self._get_electronic_descriptors(),
            DescriptorCategory.QUANTUM: self._get_quantum_descriptors(),
            DescriptorCategory.FINGERPRINTS: self._get_fingerprint_descriptors(),
            DescriptorCategory.HYBRID: self._get_hybrid_descriptors(),
            DescriptorCategory.CUSTOM: self._get_custom_descriptors(),
        }

    def _get_constitutional_descriptors(self) -> List[DescriptorInfo]:
        """Get constitutional (1D) descriptors."""
        calc = self.constitutional_calc
        return [
            DescriptorInfo("MolecularWeight", "Average molecular weight from the standard atomic weights, including the hydrogens a structure from SMILES only carries implicitly (aspirin 180.16 g/mol)", DescriptorCategory.CONSTITUTIONAL, "sum(m_i) + n_H * 1.008", "g/mol", calculation_function=calc.calc_molecular_weight),
            DescriptorInfo("AtomCount", "Total number of atoms", DescriptorCategory.CONSTITUTIONAL, "len(atoms)", "count", calculation_function=calc.calc_atom_count),
            DescriptorInfo("HeavyAtomCount", "Number of non-hydrogen atoms", DescriptorCategory.CONSTITUTIONAL, "count(non-H)", "count", calculation_function=calc.calc_heavy_atom_count),
            DescriptorInfo("BondCount", "Total number of bonds", DescriptorCategory.CONSTITUTIONAL, "len(bonds)", "count", calculation_function=calc.calc_bond_count),
            DescriptorInfo("RotatableBonds", "Number of rotatable bonds: acyclic single bonds between heavy atoms, excluding terminal ones. The standard conformational flexibility measure used by the Veber rules", DescriptorCategory.CONSTITUTIONAL, "count(rotatable)", "count", calculation_function=calc.calc_rotatable_bonds),
            DescriptorInfo("RingCount", "Number of rings in the smallest set of smallest rings (SSSR); the ring count of the molecular skeleton", DescriptorCategory.CONSTITUTIONAL, "count(rings)", "count", calculation_function=calc.calc_ring_count),
            DescriptorInfo("AromaticRingCount", "Number of aromatic rings", DescriptorCategory.CONSTITUTIONAL, "count(aromatic_rings)", "count", calculation_function=calc.calc_aromatic_ring_count),
            DescriptorInfo("HeteroRingCount", "Number of heteroatom-containing rings", DescriptorCategory.CONSTITUTIONAL, "count(hetero_rings)", "count", calculation_function=calc.calc_hetero_ring_count),
            DescriptorInfo("Ring5Count", "Number of 5-membered rings", DescriptorCategory.CONSTITUTIONAL, "count(5-rings)", "count", calculation_function=calc.calc_ring5_count),
            DescriptorInfo("Ring6Count", "Number of 6-membered rings", DescriptorCategory.CONSTITUTIONAL, "count(6-rings)", "count", calculation_function=calc.calc_ring6_count),
            DescriptorInfo("LargestRingSize", "Size of the largest ring", DescriptorCategory.CONSTITUTIONAL, "max(ring_sizes)", "count", calculation_function=calc.calc_largest_ring_size),
            DescriptorInfo("HDonorCount", "Number of hydrogen bond donors: N-H and O-H groups", DescriptorCategory.CONSTITUTIONAL, "count(H-donors)", "count", calculation_function=calc.calc_h_donor_count),
            DescriptorInfo("HAcceptorCount", "Number of hydrogen bond acceptors: nitrogen and oxygen atoms with a free lone pair", DescriptorCategory.CONSTITUTIONAL, "count(H-acceptors)", "count", calculation_function=calc.calc_h_acceptor_count),
            DescriptorInfo("LipophilicAtomCount", "Number of atoms classified as lipophilic (carbons and halogens with no polar neighbour)", DescriptorCategory.CONSTITUTIONAL, "count(lipophilic)", "count", calculation_function=calc.calc_lipophilic_count),
            DescriptorInfo("CarbonCount", "Number of carbon atoms", DescriptorCategory.CONSTITUTIONAL, "count(C)", "count", calculation_function=calc.calc_carbon_count),
            DescriptorInfo("NitrogenCount", "Number of nitrogen atoms", DescriptorCategory.CONSTITUTIONAL, "count(N)", "count", calculation_function=calc.calc_nitrogen_count),
            DescriptorInfo("OxygenCount", "Number of oxygen atoms", DescriptorCategory.CONSTITUTIONAL, "count(O)", "count", calculation_function=calc.calc_oxygen_count),
            DescriptorInfo("SulfurCount", "Number of sulfur atoms", DescriptorCategory.CONSTITUTIONAL, "count(S)", "count", calculation_function=calc.calc_sulfur_count),
            DescriptorInfo("PhosphorusCount", "Number of phosphorus atoms", DescriptorCategory.CONSTITUTIONAL, "count(P)", "count", calculation_function=calc.calc_phosphorus_count),
            DescriptorInfo("HalogenCount", "Number of halogen atoms", DescriptorCategory.CONSTITUTIONAL, "count(F,Cl,Br,I)", "count", calculation_function=calc.calc_halogen_count),
            DescriptorInfo("FluorineCount", "Number of fluorine atoms", DescriptorCategory.CONSTITUTIONAL, "count(F)", "count", calculation_function=calc.calc_fluorine_count),
            DescriptorInfo("ChlorineCount", "Number of chlorine atoms", DescriptorCategory.CONSTITUTIONAL, "count(Cl)", "count", calculation_function=calc.calc_chlorine_count),
            DescriptorInfo("BromineCount", "Number of bromine atoms", DescriptorCategory.CONSTITUTIONAL, "count(Br)", "count", calculation_function=calc.calc_bromine_count),
            DescriptorInfo("IodineCount", "Number of iodine atoms", DescriptorCategory.CONSTITUTIONAL, "count(I)", "count", calculation_function=calc.calc_iodine_count),
            DescriptorInfo("HydrogenCount", "Number of hydrogen atoms", DescriptorCategory.CONSTITUTIONAL, "count(H)", "count", calculation_function=calc.calc_hydrogen_count),
            DescriptorInfo("SingleBondCount", "Number of single bonds", DescriptorCategory.CONSTITUTIONAL, "count(single_bonds)", "count", calculation_function=calc.calc_single_bond_count),
            DescriptorInfo("DoubleBondCount", "Number of double bonds", DescriptorCategory.CONSTITUTIONAL, "count(double_bonds)", "count", calculation_function=calc.calc_double_bond_count),
            DescriptorInfo("TripleBondCount", "Number of triple bonds", DescriptorCategory.CONSTITUTIONAL, "count(triple_bonds)", "count", calculation_function=calc.calc_triple_bond_count),
            DescriptorInfo("AromaticBondCount", "Number of aromatic bonds", DescriptorCategory.CONSTITUTIONAL, "count(aromatic_bonds)", "count", calculation_function=calc.calc_aromatic_bond_count),
            DescriptorInfo("HeteroAtomCount", "Number of hetero atoms", DescriptorCategory.CONSTITUTIONAL, "count(hetero)", "count", calculation_function=calc.calc_hetero_atom_count),
            DescriptorInfo("AliphaticCarbonCount", "Number of aliphatic carbons", DescriptorCategory.CONSTITUTIONAL, "count(aliphatic_C)", "count", calculation_function=calc.calc_aliphatic_carbon_count),
            DescriptorInfo("AromaticCarbonCount", "Number of aromatic carbons", DescriptorCategory.CONSTITUTIONAL, "count(aromatic_C)", "count", calculation_function=calc.calc_aromatic_carbon_count),
            DescriptorInfo("NonRingCarbonCount", "Number of non-ring carbon atoms", DescriptorCategory.CONSTITUTIONAL, "count(C_non-ring)", "count", calculation_function=calc.calc_nonring_carbon_count),
            DescriptorInfo("SP3CarbonCount", "Number of sp3 hybridized carbon atoms", DescriptorCategory.CONSTITUTIONAL, "count(sp3_C)", "count", calculation_function=calc.calc_sp3_carbon_count),
            DescriptorInfo("SP2CarbonCount", "Number of sp2 hybridized carbon atoms", DescriptorCategory.CONSTITUTIONAL, "count(sp2_C)", "count", calculation_function=calc.calc_sp2_carbon_count),
            DescriptorInfo("SPCarbonCount", "Number of sp hybridized carbon atoms", DescriptorCategory.CONSTITUTIONAL, "count(sp_C)", "count", calculation_function=calc.calc_sp_carbon_count),
            DescriptorInfo("SP3NitrogenCount", "Number of sp3 hybridized nitrogen atoms", DescriptorCategory.CONSTITUTIONAL, "count(sp3_N)", "count", calculation_function=calc.calc_sp3_nitrogen_count),
            DescriptorInfo("SP2NitrogenCount", "Number of sp2 hybridized nitrogen atoms", DescriptorCategory.CONSTITUTIONAL, "count(sp2_N)", "count", calculation_function=calc.calc_sp2_nitrogen_count),
            DescriptorInfo("SP3OxygenCount", "Number of sp3 hybridized oxygen atoms", DescriptorCategory.CONSTITUTIONAL, "count(sp3_O)", "count", calculation_function=calc.calc_sp3_oxygen_count),
            DescriptorInfo("SP2OxygenCount", "Number of sp2 hybridized oxygen atoms", DescriptorCategory.CONSTITUTIONAL, "count(sp2_O)", "count", calculation_function=calc.calc_sp2_oxygen_count),
            DescriptorInfo("FormalCharge", "Sum of the formal charges of the selected atoms, the net charge of the species", DescriptorCategory.CONSTITUTIONAL, "sum(formal_charges)", "e", calculation_function=calc.calc_formal_charge),
            DescriptorInfo("UnsaturationCount", "Degree of unsaturation (double bond equivalents): rings plus pi bonds, computed from the molecular formula", DescriptorCategory.CONSTITUTIONAL, "DBE", "count", calculation_function=calc.calc_unsaturation_count),
            DescriptorInfo("AromaticProportion", "Fraction of atoms belonging to an aromatic ring", DescriptorCategory.CONSTITUTIONAL, "aromatic/total", "fraction", calculation_function=calc.calc_aromatic_proportion),
            DescriptorInfo("AliphaticProportion", "Fraction of atoms that are not aromatic", DescriptorCategory.CONSTITUTIONAL, "aliphatic/total", "fraction", calculation_function=calc.calc_aliphatic_proportion),
            DescriptorInfo("SP3Proportion", "Fraction of carbons that are sp3 hybridised", DescriptorCategory.CONSTITUTIONAL, "sp3_C/total_C", "fraction", calculation_function=calc.calc_sp3_proportion),
            DescriptorInfo("HeteroProportion", "Fraction of atoms that are neither carbon nor hydrogen", DescriptorCategory.CONSTITUTIONAL, "hetero/total", "fraction", calculation_function=calc.calc_hetero_proportion),
            DescriptorInfo("HCRatio", "Hydrogen to carbon ratio", DescriptorCategory.CONSTITUTIONAL, "H/C", "ratio", calculation_function=calc.calc_hc_ratio),
            DescriptorInfo("NCRatio", "Nitrogen to carbon ratio", DescriptorCategory.CONSTITUTIONAL, "N/C", "ratio", calculation_function=calc.calc_nc_ratio),
            DescriptorInfo("OCRatio", "Oxygen to carbon ratio", DescriptorCategory.CONSTITUTIONAL, "O/C", "ratio", calculation_function=calc.calc_oc_ratio),
            DescriptorInfo("AmideBondCount", "Number of amide bonds (C(=O)-N), the linkage of peptides and a very common medicinal chemistry motif", DescriptorCategory.CONSTITUTIONAL, "count(amide_bonds)", "count", calculation_function=calc.calc_amide_bond_count),
            DescriptorInfo("EsterBondCount", "Number of ester linkages (C(=O)-O-C), frequent metabolic soft spots", DescriptorCategory.CONSTITUTIONAL, "count(ester_bonds)", "count", calculation_function=calc.calc_ester_bond_count),
            DescriptorInfo("CarbonylBondCount", "Number of C=O double bonds of any kind (aldehyde, ketone, acid, ester, amide)", DescriptorCategory.CONSTITUTIONAL, "count(carbonyl_bonds)", "count", calculation_function=calc.calc_carbonyl_bond_count),
            DescriptorInfo("HydroxylGroupCount", "Number of hydroxyl (-OH) groups", DescriptorCategory.CONSTITUTIONAL, "count(OH_groups)", "count", calculation_function=calc.calc_hydroxyl_group_count),
            DescriptorInfo("CarboxylGroupCount", "Number of carboxyl (-COOH) groups", DescriptorCategory.CONSTITUTIONAL, "count(COOH_groups)", "count", calculation_function=calc.calc_carboxyl_group_count),
            DescriptorInfo("AmineGroupCount", "Number of nitrogen atoms that are not part of an amide, i.e. amine-like nitrogens", DescriptorCategory.CONSTITUTIONAL, "count(amine_groups)", "count", calculation_function=calc.calc_amine_group_count),
            DescriptorInfo("MethylGroupCount", "Number of methyl (-CH3) groups", DescriptorCategory.CONSTITUTIONAL, "count(methyl_groups)", "count", calculation_function=calc.calc_methyl_group_count),
            DescriptorInfo("LongestChain", "Number of atoms in the longest carbon chain, the classic measure of molecular length", DescriptorCategory.CONSTITUTIONAL, "max(chain_length)", "count", calculation_function=calc.calc_longest_chain),
            DescriptorInfo("BranchCount", "Number of atoms carrying three or more heavy neighbours, i.e. branch points", DescriptorCategory.CONSTITUTIONAL, "count(branches)", "count", calculation_function=calc.calc_branch_count),
            DescriptorInfo("FragmentCount", "Number of disconnected fragments", DescriptorCategory.CONSTITUTIONAL, "count(fragments)", "count", calculation_function=calc.calc_fragment_count),
        ] + self._get_constitutional_extended()

    def _get_constitutional_extended(self) -> List[DescriptorInfo]:
        """Mass, atom environment, ring system and functional group descriptors."""
        calc = self.constitutional_calc
        C = DescriptorCategory.CONSTITUTIONAL
        return [
            # --- mass and size
            DescriptorInfo("ExactMass", "Mass computed from the element masses with implicit hydrogens added; differs from MolecularWeight when the structure carries implicit H", C, "sum(m_i) + n_H * 1.008", "g/mol", calculation_function=calc.calc_exact_mass),
            DescriptorInfo("HeavyAtomMolWeight", "Molecular weight of the non-hydrogen atoms only; useful for normalising properties by molecular size", C, "sum(m_i, i != H)", "g/mol", calculation_function=calc.calc_heavy_atom_mol_weight),
            DescriptorInfo("AverageAtomicMass", "Mean atomic mass per heavy atom; rises with heavy halogens and metals, a proxy for how 'heavy' the elements are", C, "MW_heavy / n_heavy", "g/mol", calculation_function=calc.calc_average_atomic_mass),
            DescriptorInfo("ValenceElectronCount", "Total number of valence electrons; correlates strongly with polarizability and dispersion interactions", C, "sum(Zv_i)", "count", calculation_function=calc.calc_valence_electron_count),
            DescriptorInfo("VdWVolumeSum", "Sum of Bondi van der Waals atomic volumes, without any packing correction; a coordinate-free size measure", C, "sum(V_vdW,i)", "Å³", calculation_function=calc.calc_vdw_volume_sum),
            # --- carbon substitution classes
            DescriptorInfo("PrimaryCarbonCount", "Carbons bonded to exactly one other carbon (methyl and chain ends)", C, "count(C with 1 C neighbour)", "count", calculation_function=calc.calc_primary_carbon_count),
            DescriptorInfo("SecondaryCarbonCount", "Carbons bonded to two other carbons (chain interior)", C, "count(C with 2 C neighbours)", "count", calculation_function=calc.calc_secondary_carbon_count),
            DescriptorInfo("TertiaryCarbonCount", "Carbons bonded to three other carbons (branch points)", C, "count(C with 3 C neighbours)", "count", calculation_function=calc.calc_tertiary_carbon_count),
            DescriptorInfo("QuaternaryCarbonCount", "Carbons bonded to four other carbons; sterically congested centres that slow metabolism", C, "count(C with 4 C neighbours)", "count", calculation_function=calc.calc_quaternary_carbon_count),
            DescriptorInfo("TerminalAtomCount", "Heavy atoms with a single heavy neighbour, i.e. the leaves of the molecular graph", C, "count(deg = 1)", "count", calculation_function=calc.calc_terminal_atom_count),
            DescriptorInfo("MaxAtomDegree", "Largest number of heavy neighbours on any atom; measures the most branched centre", C, "max(deg_i)", "count", calculation_function=calc.calc_max_atom_degree),
            DescriptorInfo("MeanAtomDegree", "Average heavy atom degree, equal to 2*bonds/atoms; grows with ring fusion and branching", C, "2 * E / V", "count", calculation_function=calc.calc_mean_atom_degree),
            # --- ring system
            DescriptorInfo("RingAtomCount", "Number of atoms that belong to at least one ring", C, "|union(rings)|", "count", calculation_function=calc.calc_ring_atom_count),
            DescriptorInfo("RingBondCount", "Number of bonds that belong to at least one ring", C, "count(ring bonds)", "count", calculation_function=calc.calc_ring_bond_count),
            DescriptorInfo("AliphaticRingCount", "Rings containing at least one non-aromatic atom", C, "count(non-aromatic rings)", "count", calculation_function=calc.calc_aliphatic_ring_count),
            DescriptorInfo("SaturatedRingCount", "Rings in which every bond is single; a marker of three-dimensional character", C, "count(all-single rings)", "count", calculation_function=calc.calc_saturated_ring_count),
            DescriptorInfo("AromaticHeterocycleCount", "Aromatic rings that contain at least one heteroatom (pyridine, imidazole, ...)", C, "count(aromatic hetero rings)", "count", calculation_function=calc.calc_aromatic_heterocycle_count),
            DescriptorInfo("AromaticCarbocycleCount", "Aromatic rings made only of carbon (benzene rings)", C, "count(aromatic all-C rings)", "count", calculation_function=calc.calc_aromatic_carbocycle_count),
            DescriptorInfo("SaturatedHeterocycleCount", "Fully saturated rings containing a heteroatom (piperidine, morpholine, ...)", C, "count(saturated hetero rings)", "count", calculation_function=calc.calc_saturated_heterocycle_count),
            DescriptorInfo("SpiroAtomCount", "Atoms shared by two rings that have no bond in common; spiro centres enforce perpendicular ring planes", C, "count(|R_i ∩ R_j| = 1)", "count", calculation_function=calc.calc_spiro_atom_count),
            DescriptorInfo("BridgeheadAtomCount", "Atoms at the junction of two rings sharing a bond or a bridge; characteristic of rigid cage systems", C, "count(shared ring junctions)", "count", calculation_function=calc.calc_bridgehead_atom_count),
            DescriptorInfo("FusedRingCount", "Number of ring pairs that share at least one bond, i.e. the number of fusion sites", C, "count(|R_i ∩ R_j| >= 2)", "count", calculation_function=calc.calc_fused_ring_count),
            DescriptorInfo("MacrocycleCount", "Rings larger than 12 atoms, which behave very differently from ordinary rings in ADME terms", C, "count(ring size > 12)", "count", calculation_function=calc.calc_macrocycle_count),
            DescriptorInfo("RingComplexity", "Ring atoms per ring; approaches the ring size for isolated rings and drops as rings become fused", C, "|union(rings)| / n_rings", "ratio", calculation_function=calc.calc_ring_complexity),
            # --- functional groups
            DescriptorInfo("NitroGroupCount", "Number of nitro groups (R-NO2), recognised in both the charge separated and the pentavalent notation", C, "count(N with 2 terminal O)", "count", calculation_function=calc.calc_nitro_group_count),
            DescriptorInfo("NitrileCount", "Number of nitrile (C#N) groups", C, "count(C#N)", "count", calculation_function=calc.calc_nitrile_count),
            DescriptorInfo("AldehydeCount", "Number of aldehyde groups: a carbonyl carbon carrying a hydrogen", C, "count(HC=O)", "count", calculation_function=calc.calc_aldehyde_count),
            DescriptorInfo("KetoneCount", "Number of ketone groups: a carbonyl carbon flanked by two carbons", C, "count(C-CO-C)", "count", calculation_function=calc.calc_ketone_count),
            DescriptorInfo("EtherCount", "Number of ether oxygens (C-O-C), excluding esters and acids", C, "count(C-O-C)", "count", calculation_function=calc.calc_ether_count),
            DescriptorInfo("ThioetherCount", "Number of thioether sulfurs (C-S-C)", C, "count(C-S-C)", "count", calculation_function=calc.calc_thioether_count),
            DescriptorInfo("ThiolCount", "Number of thiol (-SH) groups; frequent covalent-binding and metabolic liabilities", C, "count(SH)", "count", calculation_function=calc.calc_thiol_count),
            DescriptorInfo("SulfonamideCount", "Number of sulfonamide groups S(=O)(=O)N", C, "count(SO2N)", "count", calculation_function=calc.calc_sulfonamide_count),
            DescriptorInfo("SulfoneCount", "Number of sulfone groups C-S(=O)(=O)-C", C, "count(C-SO2-C)", "count", calculation_function=calc.calc_sulfone_count),
            DescriptorInfo("SulfoxideCount", "Number of sulfoxide groups C-S(=O)-C", C, "count(C-SO-C)", "count", calculation_function=calc.calc_sulfoxide_count),
            DescriptorInfo("PhosphateCount", "Number of phosphorus atoms carrying three or more oxygens (phosphate/phosphonate)", C, "count(P with >= 3 O)", "count", calculation_function=calc.calc_phosphate_count),
            DescriptorInfo("ImineCount", "Number of non-aromatic C=N double bonds (imines and Schiff bases)", C, "count(C=N)", "count", calculation_function=calc.calc_imine_count),
            DescriptorInfo("AzoCount", "Number of N=N double bonds (azo groups, common in dyes)", C, "count(N=N)", "count", calculation_function=calc.calc_azo_count),
            DescriptorInfo("UreaCount", "Number of urea/carbamate-like fragments N-C(=O)-N", C, "count(N-CO-N)", "count", calculation_function=calc.calc_urea_count),
            DescriptorInfo("GuanidineCount", "Number of guanidine groups; strongly basic and permanently charged at physiological pH", C, "count(C(N)(N)=N)", "count", calculation_function=calc.calc_guanidine_count),
            DescriptorInfo("PhenolCount", "Number of hydroxyls attached to an aromatic ring; more acidic than aliphatic alcohols", C, "count(Ar-OH)", "count", calculation_function=calc.calc_phenol_count),
            DescriptorInfo("AliphaticHydroxylCount", "Number of alcohol hydroxyls that are neither phenolic nor part of a carboxylic acid", C, "count(alkyl-OH)", "count", calculation_function=calc.calc_aliphatic_hydroxyl_count),
            DescriptorInfo("AromaticHalideCount", "Number of halogens attached directly to an aromatic ring", C, "count(Ar-X)", "count", calculation_function=calc.calc_halide_on_aromatic_count),
            DescriptorInfo("AcidicGroupCount", "Ionisable acidic groups: carboxylic, sulfonic and phosphonic acids", C, "count(COOH + SO3H + PO3H)", "count", calculation_function=calc.calc_acidic_group_count),
            DescriptorInfo("BasicGroupCount", "Basic nitrogens: aliphatic amines, amidines and guanidines, excluding amides", C, "count(basic N)", "count", calculation_function=calc.calc_basic_group_count),
            # --- normalised ratios
            DescriptorInfo("RotatableBondFraction", "Rotatable bonds per heavy atom, a size independent measure of conformational freedom", C, "n_rot / n_heavy", "ratio", calculation_function=calc.calc_rotatable_bond_fraction),
            DescriptorInfo("HeteroatomFraction", "Heteroatoms per heavy atom; a compact polarity/complexity indicator", C, "n_hetero / n_heavy", "ratio", calculation_function=calc.calc_heteroatom_ratio_heavy),
        ]

    def _get_topological_descriptors(self) -> List[DescriptorInfo]:
        """Get topological (2D) descriptors."""
        calc = self.topological_calc
        return [
            DescriptorInfo("WienerIndex", "Wiener index: the sum of topological distances over all atom pairs. The oldest topological index, it grows with size and decreases with branching, and correlates with boiling point and molar volume", DescriptorCategory.TOPOLOGICAL, "sum(distance_matrix)", "index", calculation_function=calc.calc_wiener_index),
            DescriptorInfo("ZagrebIndex1", "First Zagreb index, the sum of squared vertex degrees; a branching measure that grows sharply with substitution at a single centre", DescriptorCategory.TOPOLOGICAL, "sum(degree^2)", "index", calculation_function=calc.calc_zagreb_index_1),
            DescriptorInfo("ZagrebIndex2", "Second Zagreb index, the sum of degree products over bonds; complements the first Zagreb index by weighting bonds rather than atoms", DescriptorCategory.TOPOLOGICAL, "sum(degree_i * degree_j)", "index", calculation_function=calc.calc_zagreb_index_2),
            DescriptorInfo("BalabanIndex", "Balaban J index, an average distance-sum connectivity that stays almost constant with molecular size, so it discriminates isomers and ring systems rather than size", DescriptorCategory.TOPOLOGICAL, "J_index", "index", calculation_function=calc.calc_balaban_index),
            DescriptorInfo("HararyIndex", "Harary index, the sum of reciprocal topological distances; unlike the Wiener index it emphasises close atom pairs and therefore local compactness", DescriptorCategory.TOPOLOGICAL, "sum(1/distance)", "index", calculation_function=calc.calc_harary_index),
            DescriptorInfo("RandicIndex", "Randic branching index, the sum of 1/sqrt(d_i*d_j) over bonds; the original molecular connectivity descriptor, lower for more branched skeletons", DescriptorCategory.TOPOLOGICAL, "chi_index", "index", calculation_function=calc.calc_randic_index),
            DescriptorInfo("ConnectivityIndexChi0", "Zeroth order connectivity index computed over bonds (legacy variant kept for backward compatibility; see Chi0n for the standard Kier-Hall definition)", DescriptorCategory.TOPOLOGICAL, "chi_0", "index", calculation_function=calc.calc_connectivity_index_chi0),
            DescriptorInfo("ConnectivityIndexChi1", "First order connectivity index (legacy variant kept for backward compatibility; see Chi1n for the standard Kier-Hall definition)", DescriptorCategory.TOPOLOGICAL, "chi_1", "index", calculation_function=calc.calc_connectivity_index_chi1),
            DescriptorInfo("KappaShapeIndex1", "First order Kappa shape index without the alpha correction, comparing the bond count with the linear and complete graph limits (see Kappa1Alpha for the size corrected form)", DescriptorCategory.TOPOLOGICAL, "kappa_1", "index", calculation_function=calc.calc_kappa_shape_index_1),
            DescriptorInfo("KappaShapeIndex2", "Second order Kappa shape index without the alpha correction; grows with molecular elongation", DescriptorCategory.TOPOLOGICAL, "kappa_2", "index", calculation_function=calc.calc_kappa_shape_index_2),
            DescriptorInfo("KappaShapeIndex3", "Third order Kappa shape index without the alpha correction; reflects where the branch points sit along the skeleton", DescriptorCategory.TOPOLOGICAL, "kappa_3", "index", calculation_function=calc.calc_kappa_shape_index_3),
            DescriptorInfo("HosoyaIndex", "Hosoya Z index, a count of the independent bond matchings of the graph; a classical complexity measure related to thermodynamic stability", DescriptorCategory.TOPOLOGICAL, "Z_index", "index", calculation_function=calc.calc_hosoya_index),
            DescriptorInfo("PlattIndex", "Platt index, the sum of the degrees of the two atoms of every bond; a simple edge-based branching measure correlating with molar volume", DescriptorCategory.TOPOLOGICAL, "sum(edge_degrees)", "index", calculation_function=calc.calc_platt_index),
            DescriptorInfo("PolarityNumber", "Polarity number (Wiener polarity): the count of atom pairs exactly three bonds apart, which is the number of rotatable-bond-like arrangements", DescriptorCategory.TOPOLOGICAL, "count(3-dist_pairs)", "count", calculation_function=calc.calc_polarity_number),
            DescriptorInfo("BertzIndex", "Bertz molecular complexity index, an information measure of the bonding pattern; grows with both size and heterogeneity of the skeleton", DescriptorCategory.TOPOLOGICAL, "complexity_based_index", "index", calculation_function=calc.calc_bertz_index),
            DescriptorInfo("BonchevTrinajsticIndex", "Bonchev-Trinajstic mean information content on the vertex degree distribution, in bits per atom", DescriptorCategory.TOPOLOGICAL, "information_based_index", "index", calculation_function=calc.calc_bonchev_trinajstic_index),
            DescriptorInfo("InformationContent0", "Zeroth order information content: the entropy of the atom equivalence classes defined by degree alone", DescriptorCategory.TOPOLOGICAL, "entropy_based_index", "bits", calculation_function=calc.calc_information_content_0),
            DescriptorInfo("InformationContent1", "First order information content: entropy of the atom classes after one sphere of neighbours is taken into account", DescriptorCategory.TOPOLOGICAL, "entropy_based_index", "bits", calculation_function=calc.calc_information_content_1),
            DescriptorInfo("EccentricConnectivityIndex", "Eccentric connectivity index, the sum over atoms of degree times eccentricity; combines branching with molecular elongation and is widely used in QSAR of drug transport", DescriptorCategory.TOPOLOGICAL, "sum(degree*eccentricity)", "index", calculation_function=calc.calc_eccentric_connectivity_index),
            DescriptorInfo("PathCount3", "Number of distinct paths spanning three atoms (two bonds); a local shape descriptor", DescriptorCategory.TOPOLOGICAL, "count(3-paths)", "count", calculation_function=calc.calc_path_count_3),
            DescriptorInfo("PathCount4", "Number of distinct paths spanning four atoms (three bonds)", DescriptorCategory.TOPOLOGICAL, "count(4-paths)", "count", calculation_function=calc.calc_path_count_4),
            DescriptorInfo("AveragePathLength", "Mean topological distance between connected atom pairs; small for compact fused systems and large for extended chains", DescriptorCategory.TOPOLOGICAL, "avg_path_length", "index", calculation_function=calc.calc_average_path_length),
            DescriptorInfo("LongestShortestPath", "Graph diameter: the longest of all shortest paths, i.e. the topological length of the molecule", DescriptorCategory.TOPOLOGICAL, "graph_diameter", "count", calculation_function=calc.calc_longest_shortest_path),
            DescriptorInfo("CyclomaticNumber", "Cyclomatic number (circuit rank) E - V + C: the number of independent rings, i.e. how many bonds must be cut to make the molecule acyclic", DescriptorCategory.TOPOLOGICAL, "E - V + C", "count", calculation_function=calc.calc_cyclomatic_number),
            DescriptorInfo("ElectrotopologicalStateIndex", "Sum of intrinsic electrotopological states (atomic number over degree); a combined electronic and topological measure of atom accessibility", DescriptorCategory.TOPOLOGICAL, "electronic_topology_index", "index", calculation_function=calc.calc_electrotopological_state_index),
            DescriptorInfo("MolecularTopologicalIndex", "Composite topological index combining the Wiener and first Zagreb indices into a single size-and-branching measure", DescriptorCategory.TOPOLOGICAL, "topological_complexity", "index", calculation_function=calc.calc_molecular_topological_index),
            DescriptorInfo("HyperWienerIndex", "Hyper-Wiener index, which adds squared distances to the Wiener sum; discriminates isomers that the Wiener index cannot separate", DescriptorCategory.TOPOLOGICAL, "sum(distance + distance^2)", "index", calculation_function=calc.calc_hyper_wiener_index),
        ] + self._get_topological_extended()

    def _get_topological_extended(self) -> List[DescriptorInfo]:
        """Connectivity, shape, spectral and autocorrelation descriptors.

        All of them are evaluated on the hydrogen suppressed graph, which is how
        the original definitions are stated. """
        calc = self.topological_calc
        T = DescriptorCategory.TOPOLOGICAL
        return [
            # --- Kier-Hall connectivity (chi) indices
            DescriptorInfo("Chi0n", "Kier-Hall simple connectivity index of order 0: sums 1/sqrt(delta) over atoms, where delta is the heavy atom degree. Grows with size and with the number of unbranched atoms", T, "sum(1/sqrt(d_i))", "index", calculation_function=calc.calc_chi0n),
            DescriptorInfo("Chi1n", "First order simple connectivity index (Randic index) over bonds; the classic branching descriptor, lower for more branched skeletons", T, "sum(1/sqrt(d_i*d_j)) over bonds", "index", calculation_function=calc.calc_chi1n),
            DescriptorInfo("Chi2n", "Second order simple connectivity index over two-bond paths; encodes the local environment one bond further out", T, "sum over 2-paths", "index", calculation_function=calc.calc_chi2n),
            DescriptorInfo("Chi3n", "Third order simple connectivity index over three-bond paths", T, "sum over 3-paths", "index", calculation_function=calc.calc_chi3n),
            DescriptorInfo("Chi4n", "Fourth order simple connectivity index over four-bond paths; sensitive to medium range shape", T, "sum over 4-paths", "index", calculation_function=calc.calc_chi4n),
            DescriptorInfo("Chi0v", "Valence connectivity index of order 0; uses delta_v = (Zv - h)/(Z - Zv - 1), so heteroatoms and their lone pairs are distinguished from carbon", T, "sum(1/sqrt(dv_i))", "index", calculation_function=calc.calc_chi0v),
            DescriptorInfo("Chi1v", "First order valence connectivity index; the electronic counterpart of the Randic index and a strong predictor of logP and boiling point", T, "sum(1/sqrt(dv_i*dv_j))", "index", calculation_function=calc.calc_chi1v),
            DescriptorInfo("Chi2v", "Second order valence connectivity index over two-bond paths", T, "sum over 2-paths (valence)", "index", calculation_function=calc.calc_chi2v),
            DescriptorInfo("Chi3v", "Third order valence connectivity index over three-bond paths", T, "sum over 3-paths (valence)", "index", calculation_function=calc.calc_chi3v),
            DescriptorInfo("Chi4v", "Fourth order valence connectivity index over four-bond paths", T, "sum over 4-paths (valence)", "index", calculation_function=calc.calc_chi4v),
            DescriptorInfo("Chi3Cluster", "Third order cluster connectivity: sums over three-branch stars, so it counts tertiary branch points specifically", T, "sum over 3-stars", "index", calculation_function=calc.calc_chi3_cluster),
            DescriptorInfo("Chi4Cluster", "Fourth order cluster connectivity over four-branch stars, i.e. quaternary centres", T, "sum over 4-stars", "index", calculation_function=calc.calc_chi4_cluster),
            # --- Kier shape indices
            DescriptorInfo("Kappa1Alpha", "First Kier shape index with the alpha correction for atom size and hybridisation; compares the molecule with the linear and the fully connected limits", T, "(A+a)(A+a-1)^2 / (P1+a)^2", "index", calculation_function=calc.calc_kappa1_alpha),
            DescriptorInfo("Kappa2Alpha", "Second Kier shape index; increases with molecular elongation and decreases with branching and cyclisation", T, "(A+a-1)(A+a-2)^2 / (P2+a)^2", "index", calculation_function=calc.calc_kappa2_alpha),
            DescriptorInfo("Kappa3Alpha", "Third Kier shape index, sensitive to the position of branch points along the skeleton", T, "(A+a-1)(A+a-3)^2 / (P3+a)^2", "index", calculation_function=calc.calc_kappa3_alpha),
            DescriptorInfo("KierFlexibility", "Kier molecular flexibility index Phi = kappa1*kappa2/A; roughly 0 for rigid fused systems and large for long flexible chains", T, "kappa1 * kappa2 / A", "index", calculation_function=calc.calc_kier_flexibility),
            # --- modern degree based indices
            DescriptorInfo("ABCIndex", "Atom-bond connectivity index; correlates with the strain energy of branched alkanes", T, "sum(sqrt((d_i+d_j-2)/(d_i*d_j)))", "index", calculation_function=calc.calc_abc_index),
            DescriptorInfo("AugmentedZagrebIndex", "Augmented Zagreb index, a heavily branch-weighted degree index used in heat-of-formation models", T, "sum((d_i*d_j/(d_i+d_j-2))^3)", "index", calculation_function=calc.calc_augmented_zagreb_index),
            DescriptorInfo("ForgottenIndex", "Forgotten topological index, the sum of cubed vertex degrees; emphasises highly connected atoms", T, "sum(d_i^3)", "index", calculation_function=calc.calc_forgotten_index),
            DescriptorInfo("ModifiedZagrebIndex", "Modified first Zagreb index, the sum of inverse squared degrees; weights terminal atoms most", T, "sum(1/d_i^2)", "index", calculation_function=calc.calc_modified_zagreb_index),
            DescriptorInfo("SumConnectivityIndex", "Sum-connectivity index, a variant of the Randic index using the degree sum instead of the product", T, "sum(1/sqrt(d_i+d_j))", "index", calculation_function=calc.calc_sum_connectivity_index),
            DescriptorInfo("GeometricArithmeticIndex", "Geometric-arithmetic index; equals the bond count for regular graphs and drops as degrees become uneven", T, "sum(2*sqrt(d_i*d_j)/(d_i+d_j))", "index", calculation_function=calc.calc_geometric_arithmetic_index),
            DescriptorInfo("HarmonicIndex", "Harmonic index, the harmonic mean counterpart of the Randic index", T, "sum(2/(d_i+d_j))", "index", calculation_function=calc.calc_harmonic_index),
            DescriptorInfo("NarumiSimpleIndex", "Logarithm of the product of all vertex degrees (the product itself overflows for large molecules)", T, "sum(ln d_i)", "index", calculation_function=calc.calc_narumi_simple_index),
            DescriptorInfo("NarumiHarmonicIndex", "Harmonic mean of the vertex degrees", T, "A / sum(1/d_i)", "index", calculation_function=calc.calc_narumi_harmonic_index),
            DescriptorInfo("NarumiGeometricIndex", "Geometric mean of the vertex degrees", T, "(prod d_i)^(1/A)", "index", calculation_function=calc.calc_narumi_geometric_index),
            # --- distance based indices
            DescriptorInfo("SchultzIndex", "Schultz molecular topological index, combining the adjacency and distance matrices weighted by degree", T, "sum(d_i * (A+D)_ij)", "index", calculation_function=calc.calc_schultz_index),
            DescriptorInfo("GutmanIndex", "Gutman (Schultz of the second kind) index: degree weighted sum of topological distances", T, "sum(d_i*d_j*D_ij)", "index", calculation_function=calc.calc_gutman_index),
            DescriptorInfo("SzegedIndex", "Szeged index; for every bond it multiplies the number of atoms closer to each end, generalising the Wiener index to cyclic graphs", T, "sum(n_u * n_v) over bonds", "index", calculation_function=calc.calc_szeged_index),
            DescriptorInfo("PIIndex", "Padmakar-Ivan index, the additive counterpart of the Szeged index; discriminates ring systems well", T, "sum(n_u + n_v) over bonds", "index", calculation_function=calc.calc_pi_index),
            DescriptorInfo("MeanWienerIndex", "Wiener index divided by the number of atom pairs, i.e. the mean topological distance; a size independent compactness measure", T, "W / C(n,2)", "bonds", calculation_function=calc.calc_mean_wiener_index),
            DescriptorInfo("KirchhoffIndex", "Resistance distance (quasi-Wiener) index computed from the Laplacian spectrum; unlike the Wiener index it accounts for every path between two atoms, not just the shortest one", T, "n * sum(1/mu_k)", "index", calculation_function=calc.calc_kirchhoff_index),
            DescriptorInfo("TopologicalRadius", "Smallest atomic eccentricity, the graph radius", T, "min(ecc_i)", "bonds", calculation_function=calc.calc_topological_radius),
            DescriptorInfo("AverageEccentricity", "Mean atomic eccentricity; grows linearly with chain length and slowly for compact ring systems", T, "mean(ecc_i)", "bonds", calculation_function=calc.calc_average_eccentricity),
            DescriptorInfo("PetitjeanIndex", "Topological shape index (diameter - radius)/radius; 0 for perfectly symmetric graphs such as benzene and approaching 1 for long chains", T, "(D - R) / R", "index", calculation_function=calc.calc_petitjean_index),
            DescriptorInfo("EccentricDistanceSum", "Sum over atoms of eccentricity times the total distance to all other atoms; combines shape and size", T, "sum(ecc_i * sum_j D_ij)", "index", calculation_function=calc.calc_eccentric_distance_sum),
            # --- spectral descriptors
            DescriptorInfo("GraphEnergy", "Graph energy, the sum of absolute adjacency eigenvalues; tracks the total pi bonding capacity of the skeleton", T, "sum(|lambda_i|)", "index", calculation_function=calc.calc_graph_energy),
            DescriptorInfo("SpectralRadius", "Largest adjacency eigenvalue; bounded by the maximum degree and larger for densely fused systems", T, "max(|lambda_i|)", "index", calculation_function=calc.calc_spectral_radius),
            DescriptorInfo("EstradaIndex", "Estrada index, the sum of exponentials of the adjacency eigenvalues; a folding/compactness measure dominated by short closed walks", T, "sum(exp(lambda_i))", "index", calculation_function=calc.calc_estrada_index),
            DescriptorInfo("LaplacianSpectralRadius", "Largest Laplacian eigenvalue, an upper bound related to the maximum degree and the connectivity of the graph", T, "max(mu_i)", "index", calculation_function=calc.calc_laplacian_spectral_radius),
            DescriptorInfo("AlgebraicConnectivity", "Fiedler value, the second smallest Laplacian eigenvalue; measures how hard the molecular graph is to cut in two and is 0 for disconnected structures", T, "mu_2", "index", calculation_function=calc.calc_algebraic_connectivity),
            DescriptorInfo("LogSpanningTreeCount", "Log10 of the number of spanning trees (Kirchhoff matrix-tree theorem); 0 for acyclic molecules and rising with ring count and fusion", T, "log10(prod(mu_k)/n)", "index", calculation_function=calc.calc_log_spanning_tree_count),
            DescriptorInfo("WalkCount2", "Number of walks of length 2 in the molecular graph", T, "sum(A^2)", "count", calculation_function=calc.calc_walk_count_2),
            DescriptorInfo("WalkCount3", "Number of walks of length 3; together with the lower orders it forms a classic complexity series", T, "sum(A^3)", "count", calculation_function=calc.calc_walk_count_3),
            DescriptorInfo("WalkCount4", "Number of walks of length 4", T, "sum(A^4)", "count", calculation_function=calc.calc_walk_count_4),
            DescriptorInfo("SelfReturningWalk3", "Trace of A^3, exactly six times the number of three-membered rings", T, "tr(A^3)", "count", calculation_function=calc.calc_self_returning_walk_3),
            # --- 2D autocorrelation
            DescriptorInfo("ATSElectronegativity1", "Moreau-Broto autocorrelation of electronegativity at topological distance 1: sums the product of the (carbon-scaled) electronegativities of bonded atoms", T, "sum(w_i*w_j), d_ij = 1", "index", calculation_function=calc.calc_ats_en_1),
            DescriptorInfo("ATSElectronegativity2", "Electronegativity autocorrelation at distance 2 (atoms separated by two bonds)", T, "sum(w_i*w_j), d_ij = 2", "index", calculation_function=calc.calc_ats_en_2),
            DescriptorInfo("ATSElectronegativity3", "Electronegativity autocorrelation at distance 3", T, "sum(w_i*w_j), d_ij = 3", "index", calculation_function=calc.calc_ats_en_3),
            DescriptorInfo("ATSPolarizability1", "Autocorrelation of atomic polarizability at distance 1; encodes where the easily polarised atoms sit relative to each other", T, "sum(w_i*w_j), d_ij = 1", "index", calculation_function=calc.calc_ats_polarizability_1),
            DescriptorInfo("ATSPolarizability2", "Atomic polarizability autocorrelation at distance 2", T, "sum(w_i*w_j), d_ij = 2", "index", calculation_function=calc.calc_ats_polarizability_2),
            DescriptorInfo("ATSPolarizability3", "Atomic polarizability autocorrelation at distance 3", T, "sum(w_i*w_j), d_ij = 3", "index", calculation_function=calc.calc_ats_polarizability_3),
            DescriptorInfo("ATSVolume1", "Autocorrelation of van der Waals volume at distance 1; a steric analogue of the electronegativity autocorrelations", T, "sum(w_i*w_j), d_ij = 1", "index", calculation_function=calc.calc_ats_volume_1),
            DescriptorInfo("ATSVolume2", "Van der Waals volume autocorrelation at distance 2", T, "sum(w_i*w_j), d_ij = 2", "index", calculation_function=calc.calc_ats_volume_2),
            DescriptorInfo("ATSVolume3", "Van der Waals volume autocorrelation at distance 3", T, "sum(w_i*w_j), d_ij = 3", "index", calculation_function=calc.calc_ats_volume_3),
            DescriptorInfo("MoranElectronegativity1", "Moran spatial autocorrelation of electronegativity at distance 1, normalised to [-1, 1]; positive means neighbouring atoms have similar electronegativity", T, "Moran I, d = 1", "index", calculation_function=calc.calc_moran_en_1),
            DescriptorInfo("MoranElectronegativity2", "Moran autocorrelation of electronegativity at distance 2", T, "Moran I, d = 2", "index", calculation_function=calc.calc_moran_en_2),
            DescriptorInfo("MoranMass1", "Moran autocorrelation of atomic mass between bonded atoms", T, "Moran I, d = 1", "index", calculation_function=calc.calc_moran_mass_1),
            DescriptorInfo("GearyElectronegativity1", "Geary autocorrelation of electronegativity at distance 1; near 0 for smooth property distributions and above 1 when bonded atoms differ strongly", T, "Geary C, d = 1", "index", calculation_function=calc.calc_geary_en_1),
            DescriptorInfo("GearyElectronegativity2", "Geary autocorrelation of electronegativity at distance 2", T, "Geary C, d = 2", "index", calculation_function=calc.calc_geary_en_2),
            DescriptorInfo("GearyMass1", "Geary autocorrelation of atomic mass between bonded atoms", T, "Geary C, d = 1", "index", calculation_function=calc.calc_geary_mass_1),
            # --- Galvez topological charge
            DescriptorInfo("TopologicalCharge1", "Galvez topological charge index at distance 1; measures the net charge transfer implied by the topology between bonded atoms", T, "sum|CT_ij|, d_ij = 1", "index", calculation_function=calc.calc_topological_charge_1),
            DescriptorInfo("TopologicalCharge2", "Galvez topological charge index at distance 2", T, "sum|CT_ij|, d_ij = 2", "index", calculation_function=calc.calc_topological_charge_2),
            DescriptorInfo("TopologicalCharge3", "Galvez topological charge index at distance 3", T, "sum|CT_ij|, d_ij = 3", "index", calculation_function=calc.calc_topological_charge_3),
            DescriptorInfo("MeanTopologicalCharge", "Mean Galvez charge index over distances 1-3, normalised by the number of atoms minus one", T, "sum(GGI_k)/(n-1)", "index", calculation_function=calc.calc_mean_topological_charge),
            # --- information content
            DescriptorInfo("AtomTypeInformation", "Shannon entropy of the element distribution, in bits per atom; 0 for a hydrocarbon skeleton and rising with elemental diversity", T, "-sum(p*log2 p)", "bits", calculation_function=calc.calc_atom_type_information),
            DescriptorInfo("TotalAtomInformation", "Total information content of the element distribution (entropy times atom count)", T, "n * -sum(p*log2 p)", "bits", calculation_function=calc.calc_total_atom_information),
            DescriptorInfo("BondTypeInformation", "Shannon entropy of the bond order distribution (single/double/triple/aromatic)", T, "-sum(p*log2 p)", "bits", calculation_function=calc.calc_bond_type_information),
            DescriptorInfo("DistanceInformation", "Bonchev-Trinajstic information index on the distribution of topological distances; a well established measure of structural complexity", T, "n_pairs * -sum(p*log2 p)", "bits", calculation_function=calc.calc_distance_information),
            DescriptorInfo("VertexDegreeInformation", "Shannon entropy of the vertex degree distribution, a compact branching-diversity measure", T, "-sum(p*log2 p)", "bits", calculation_function=calc.calc_vertex_degree_information),
        ]

    def _get_geometric_descriptors(self) -> List[DescriptorInfo]:
        """Get geometric (3D) descriptors."""
        calc = self.geometric_calc
        return [
            DescriptorInfo("SASA", "Solvent accessible surface area from the Shrake-Rupley algorithm with a 1.4 A probe; the surface a water molecule can touch", DescriptorCategory.GEOMETRIC, "surface_area_calculation", "Å²", calculation_function=calc.calc_sasa),
            DescriptorInfo("MolecularVolume", "Van der Waals volume obtained by numerical integration over the atomic spheres", DescriptorCategory.GEOMETRIC, "volume_calculation", "Å³", calculation_function=calc.calc_molecular_volume),
            DescriptorInfo("RadiusOfGyration", "Root mean square distance of the atoms from the centroid; the compactness measure used for polymers and conformer comparison", DescriptorCategory.GEOMETRIC, "sqrt(sum((r - r_center)^2 / N)", "Å", calculation_function=calc.calc_radius_of_gyration),
            DescriptorInfo("Asphericity", "Asphericity from the gyration tensor eigenvalues: 0 for a spherically symmetric molecule and 1 for a perfectly linear one", DescriptorCategory.GEOMETRIC, "eigenvalue_based", "dimensionless", calculation_function=calc.calc_asphericity),
            DescriptorInfo("Eccentricity", "Geometric eccentricity from the gyration tensor; 0 for a sphere and approaching 1 for elongated shapes", DescriptorCategory.GEOMETRIC, "eigenvalue_based", "dimensionless", calculation_function=calc.calc_eccentricity),
            DescriptorInfo("PrincipalMoment1", "Largest eigenvalue of the (unweighted) gyration tensor, the spatial extent along the principal axis", DescriptorCategory.GEOMETRIC, "eigenvalue1", "Å²", calculation_function=calc.calc_principal_moment_1),
            DescriptorInfo("PrincipalMoment2", "Second eigenvalue of the gyration tensor", DescriptorCategory.GEOMETRIC, "eigenvalue2", "Å²", calculation_function=calc.calc_principal_moment_2),
            DescriptorInfo("PrincipalMoment3", "Smallest eigenvalue of the gyration tensor; near zero for planar molecules", DescriptorCategory.GEOMETRIC, "eigenvalue3", "Å²", calculation_function=calc.calc_principal_moment_3),
            DescriptorInfo("MolecularDiameter", "Largest through-space distance between any two atoms, the diameter of the enclosing sphere", DescriptorCategory.GEOMETRIC, "max_distance", "Å", calculation_function=calc.calc_molecular_diameter),
        ] + self._get_geometric_extended()

    def _get_geometric_extended(self) -> List[DescriptorInfo]:
        """Inertia, shape ratios, size along the principal axes, density.

        These need a geometry; when the structure has none, coordinates are
        generated once with the built-in 3D generator. """
        calc = self.geometric_calc
        G = DescriptorCategory.GEOMETRIC
        return [
            DescriptorInfo("InertiaMomentA", "Smallest principal moment of inertia (mass weighted); small for rod shaped molecules", G, "eigenvalue 1 of the inertia tensor", "amu·Å²", calculation_function=calc.calc_inertia_a),
            DescriptorInfo("InertiaMomentB", "Intermediate principal moment of inertia", G, "eigenvalue 2 of the inertia tensor", "amu·Å²", calculation_function=calc.calc_inertia_b),
            DescriptorInfo("InertiaMomentC", "Largest principal moment of inertia", G, "eigenvalue 3 of the inertia tensor", "amu·Å²", calculation_function=calc.calc_inertia_c),
            DescriptorInfo("NPR1", "Normalised principal moment ratio I1/I3. Together with NPR2 it places a molecule on the rod-disc-sphere triangle: rods sit near (0.5, 0.5), discs near (0, 1) and spheres near (1, 1)", G, "I1 / I3", "ratio", calculation_function=calc.calc_npr1),
            DescriptorInfo("NPR2", "Normalised principal moment ratio I2/I3, the second coordinate of the principal moments of inertia shape triangle", G, "I2 / I3", "ratio", calculation_function=calc.calc_npr2),
            DescriptorInfo("InertialShapeFactor", "I2/(I1*I3); large for elongated molecules and small for compact globular ones", G, "I2 / (I1 * I3)", "1/(amu·Å²)", calculation_function=calc.calc_inertial_shape_factor),
            DescriptorInfo("SpherocityIndex", "Three times the smallest gyration eigenvalue over their sum; 1 for a perfect sphere and 0 for a planar or linear molecule", G, "3*L1 / (L1+L2+L3)", "dimensionless", calculation_function=calc.calc_spherocity_index),
            DescriptorInfo("MolecularLength", "Extent of the molecule along its longest principal axis", G, "max-min projection on axis 1", "Å", calculation_function=calc.calc_molecular_length),
            DescriptorInfo("MolecularWidth", "Extent along the second principal axis", G, "max-min projection on axis 2", "Å", calculation_function=calc.calc_molecular_width),
            DescriptorInfo("MolecularThickness", "Extent along the shortest principal axis; near 0 for flat aromatic systems", G, "max-min projection on axis 3", "Å", calculation_function=calc.calc_molecular_thickness),
            DescriptorInfo("LengthToWidthRatio", "Aspect ratio of the molecule, length divided by width", G, "L / W", "ratio", calculation_function=calc.calc_length_to_width_ratio),
            DescriptorInfo("BoundingBoxVolume", "Volume of the box aligned with the principal axes that encloses the molecule", G, "L * W * T", "Å³", calculation_function=calc.calc_bounding_box_volume),
            DescriptorInfo("Span", "Largest distance from the centroid to any atom, the radius of the enclosing sphere", G, "max|r_i - centroid|", "Å", calculation_function=calc.calc_span),
            DescriptorInfo("PlaneOfBestFit", "Mean distance of the atoms from their best fit plane; 0 for planar molecules and rising with three-dimensional character", G, "mean|(r_i - r0)·n|", "Å", calculation_function=calc.calc_plane_of_best_fit),
            DescriptorInfo("MassWeightedRadiusOfGyration", "Radius of gyration weighted by atomic mass, the physical definition used in polymer and dynamics work", G, "sqrt(sum(m_i*r_i²)/sum(m_i))", "Å", calculation_function=calc.calc_mass_weighted_radius_of_gyration),
            DescriptorInfo("MeanAtomicDistance3D", "Mean through-space distance over all atom pairs", G, "mean(|r_i - r_j|)", "Å", calculation_function=calc.calc_mean_atomic_distance),
            DescriptorInfo("GeometricWienerIndex", "Sum of all through-space interatomic distances, the 3D analogue of the Wiener index", G, "sum(|r_i - r_j|)", "Å", calculation_function=calc.calc_geometric_wiener_index),
            DescriptorInfo("GravitationalIndex", "Katritzky gravitational index over all atom pairs; encodes mass distribution and correlates with boiling point and density", G, "sum(m_i*m_j / r_ij²)", "amu²/Å²", calculation_function=calc.calc_gravitational_index),
            DescriptorInfo("GravitationalIndexBonded", "Gravitational index restricted to bonded atom pairs, a local mass-distribution term", G, "sum over bonds(m_i*m_j / r_ij²)", "amu²/Å²", calculation_function=calc.calc_gravitational_index_bonded),
            DescriptorInfo("SurfaceToVolumeRatio", "Solvent accessible surface divided by molecular volume; high for extended or branched shapes", G, "SASA / V", "1/Å", calculation_function=calc.calc_surface_to_volume_ratio),
            DescriptorInfo("Globularity", "Surface of the sphere of equal volume divided by the real surface; 1 for a sphere and below 1 for anything rugged or elongated", G, "(36*pi*V²)^(1/3) / SASA", "dimensionless", calculation_function=calc.calc_globularity),
            DescriptorInfo("MolecularDensity", "Molecular weight per unit molecular volume", G, "MW / V", "g/mol/Å³", calculation_function=calc.calc_molecular_density),
        ]

    def _get_electronic_descriptors(self) -> List[DescriptorInfo]:
        """Get electronic descriptors."""
        calc = self.electronic_calc
        return [
            DescriptorInfo("TotalCharge", "Sum of the Gasteiger partial charges over the selection; close to the formal charge for a whole molecule", DescriptorCategory.ELECTRONIC, "sum(partial_charges)", "e", calculation_function=calc.calc_total_charge),
            DescriptorInfo("MaxPartialCharge", "Most positive Gasteiger partial charge, marking the strongest electrophilic centre", DescriptorCategory.ELECTRONIC, "max(partial_charges)", "e", calculation_function=calc.calc_max_partial_charge),
            DescriptorInfo("MinPartialCharge", "Most negative Gasteiger partial charge, marking the strongest nucleophilic centre", DescriptorCategory.ELECTRONIC, "min(partial_charges)", "e", calculation_function=calc.calc_min_partial_charge),
            DescriptorInfo("MaxPositiveCharge", "Largest positive partial charge, a hydrogen bond donor strength indicator", DescriptorCategory.ELECTRONIC, "max_positive", "e", calculation_function=calc.calc_max_positive_charge),
            DescriptorInfo("MaxNegativeCharge", "Largest negative partial charge, a hydrogen bond acceptor strength indicator", DescriptorCategory.ELECTRONIC, "max_negative", "e", calculation_function=calc.calc_max_negative_charge),
            DescriptorInfo("DipoleMoment", "Molecular dipole moment computed from the Gasteiger charges and the 3D geometry; drives solubility and crystal packing", DescriptorCategory.ELECTRONIC, "sqrt(sum(charge * position)^2)", "Debye", calculation_function=calc.calc_dipole_moment),
            DescriptorInfo("PolarSurfaceArea", "Solvent accessible surface contributed by the polar atoms (N, O and their hydrogens), computed from the real 3D surface rather than from fragment values", DescriptorCategory.ELECTRONIC, "surface_area_of_polar_atoms", "Å²", calculation_function=calc.calc_polar_surface_area),
            DescriptorInfo("APolarSurfaceArea", "Solvent accessible surface contributed by the non-polar atoms, the hydrophobic part of the surface", DescriptorCategory.ELECTRONIC, "surface_area_of_apolar_atoms", "Å²", calculation_function=calc.calc_apolar_surface_area),
            DescriptorInfo("MeanAbsoluteCharge", "Mean absolute partial charge per atom, a compact measure of overall molecular polarity", DescriptorCategory.ELECTRONIC, "mean(|charge|)", "e", calculation_function=calc.calc_mean_absolute_charge),
        ] + self._get_electronic_extended()

    def _get_electronic_extended(self) -> List[DescriptorInfo]:
        """Charge statistics and the charged partial surface area (CPSA) family.

        Gasteiger charges are computed on demand when the molecule has none. """
        calc = self.electronic_calc
        E = DescriptorCategory.ELECTRONIC
        return [
            DescriptorInfo("TotalAbsoluteCharge", "Sum of the absolute partial charges; the electronic 'size' of the molecule and a measure of overall charge separation", E, "sum(|q_i|)", "e", calculation_function=calc.calc_total_absolute_charge),
            DescriptorInfo("ChargeVariance", "Variance of the atomic partial charges; large when the molecule mixes strongly polarised and neutral regions", E, "var(q_i)", "e²", calculation_function=calc.calc_charge_variance),
            DescriptorInfo("SumSquaredCharges", "Sum of squared partial charges, the electrostatic self-energy term of the charge distribution", E, "sum(q_i²)", "e²", calculation_function=calc.calc_sum_squared_charges),
            DescriptorInfo("SubmolecularPolarity", "Largest difference between any two atomic charges (Todeschini submolecular polarity parameter); a direct measure of internal charge separation", E, "max(q) - min(q)", "e", calculation_function=calc.calc_submolecular_polarity),
            DescriptorInfo("RelativePositiveCharge", "Most positive atomic charge divided by the total positive charge; identifies whether positive charge is concentrated on one atom or spread out", E, "q_max / sum(q > 0)", "fraction", calculation_function=calc.calc_relative_positive_charge),
            DescriptorInfo("RelativeNegativeCharge", "Most negative atomic charge divided by the total negative charge", E, "|q_min| / sum|q < 0|", "fraction", calculation_function=calc.calc_relative_negative_charge),
            DescriptorInfo("MeanElectronegativity", "Mean Pauling electronegativity of the atoms present", E, "mean(EN_i)", "Pauling", calculation_function=calc.calc_mean_electronegativity),
            DescriptorInfo("ElectronegativityVariance", "Variance of the atomic electronegativities, a composition-only polarity indicator", E, "var(EN_i)", "Pauling²", calculation_function=calc.calc_electronegativity_variance),
            DescriptorInfo("MaxBondPolarity", "Largest electronegativity difference across any bond; flags the most polarised (and often most reactive) bond", E, "max|EN_i - EN_j|", "Pauling", calculation_function=calc.calc_max_bond_polarity),
            DescriptorInfo("PolarAtomCharge", "Total absolute charge carried by N, O, S and P atoms, the hydrogen bonding centres", E, "sum|q_i| over N,O,S,P", "e", calculation_function=calc.calc_topological_polar_charge),
            DescriptorInfo("PPSA1", "Partial positive surface area: solvent accessible surface of all positively charged atoms (Stanton-Jurs CPSA family)", E, "sum(SA_i), q_i > 0", "Å²", calculation_function=calc.calc_ppsa1),
            DescriptorInfo("PNSA1", "Partial negative surface area: solvent accessible surface of all negatively charged atoms", E, "sum(SA_i), q_i < 0", "Å²", calculation_function=calc.calc_pnsa1),
            DescriptorInfo("PPSA2", "Total charge weighted positive surface area: the positive surface multiplied by the total positive charge", E, "sum(q+) * sum(SA+)", "e·Å²", calculation_function=calc.calc_ppsa2),
            DescriptorInfo("PNSA2", "Total charge weighted negative surface area", E, "sum(q-) * sum(SA-)", "e·Å²", calculation_function=calc.calc_pnsa2),
            DescriptorInfo("PPSA3", "Atomic charge weighted positive surface area; each positive atom contributes its own charge times its own surface", E, "sum(q_i * SA_i), q_i > 0", "e·Å²", calculation_function=calc.calc_ppsa3),
            DescriptorInfo("PNSA3", "Atomic charge weighted negative surface area, one of the strongest CPSA descriptors for hydrogen bond acceptor strength", E, "sum(q_i * SA_i), q_i < 0", "e·Å²", calculation_function=calc.calc_pnsa3),
            DescriptorInfo("DPSA1", "Difference between the positive and negative partial surface areas", E, "PPSA1 - PNSA1", "Å²", calculation_function=calc.calc_dpsa1),
            DescriptorInfo("FPSA1", "Fractional positive surface area, PPSA1 over the total surface; size independent so it compares across molecules", E, "PPSA1 / SASA", "fraction", calculation_function=calc.calc_fpsa1),
            DescriptorInfo("FNSA1", "Fractional negative surface area, PNSA1 over the total surface", E, "PNSA1 / SASA", "fraction", calculation_function=calc.calc_fnsa1),
            DescriptorInfo("WPSA1", "Surface weighted positive surface area, PPSA1 times the total surface divided by 1000", E, "PPSA1 * SASA / 1000", "Å⁴", calculation_function=calc.calc_wpsa1),
            DescriptorInfo("WNSA1", "Surface weighted negative surface area", E, "PNSA1 * SASA / 1000", "Å⁴", calculation_function=calc.calc_wnsa1),
            DescriptorInfo("RPCS", "Relative positive charged surface area: the surface of the most positive atom scaled by its share of the positive charge; a hydrogen bond donor strength descriptor", E, "SA(q_max) * q_max / sum(q+)", "Å²", calculation_function=calc.calc_rpcs),
            DescriptorInfo("RNCS", "Relative negative charged surface area, the acceptor counterpart of RPCS", E, "SA(q_min) * |q_min| / sum|q-|", "Å²", calculation_function=calc.calc_rncs),
        ]

    def _get_quantum_descriptors(self) -> List[DescriptorInfo]:
        """Quantum descriptors from a simple Huckel MO calculation on the pi system."""
        calc = self.quantum_calc
        Q = DescriptorCategory.QUANTUM
        return [
            DescriptorInfo("HOMOEnergy", "Energy of the highest occupied molecular orbital from a Huckel calculation on the pi system (alpha = -6.15 eV, beta = -3.32 eV). Molecules without a pi system fall back to the ionisation energy of their highest lone pair or sigma bond, so the value stays molecule dependent", Q, "alpha + x_HOMO * beta", "eV", calculation_function=calc.calc_homo_energy),
            DescriptorInfo("LUMOEnergy", "Energy of the lowest unoccupied molecular orbital; the more negative, the easier the molecule accepts an electron", Q, "alpha + x_LUMO * beta", "eV", calculation_function=calc.calc_lumo_energy),
            DescriptorInfo("HOMOLUMOGap", "Frontier orbital gap. Large for saturated and isolated systems, small for extended conjugation; it governs UV absorption and kinetic stability", Q, "E_LUMO - E_HOMO", "eV", calculation_function=calc.calc_homo_lumo_gap),
            DescriptorInfo("ChemicalPotential", "Electronic chemical potential mu = (E_HOMO + E_LUMO)/2, the negative of the Mulliken electronegativity; measures the escaping tendency of electrons", Q, "(E_HOMO + E_LUMO)/2", "eV", calculation_function=calc.calc_chemical_potential),
            DescriptorInfo("ChemicalHardness", "Chemical hardness eta = (E_LUMO - E_HOMO)/2; hard molecules resist charge transfer, soft ones polarise easily", Q, "(E_LUMO - E_HOMO)/2", "eV", calculation_function=calc.calc_chemical_hardness),
            DescriptorInfo("Softness", "Global softness, the reciprocal of the hardness", Q, "1 / eta", "1/eV", calculation_function=calc.calc_softness),
            DescriptorInfo("IonizationPotential", "Koopmans ionisation potential, IP = -E(HOMO). For benzene this gives 9.47 eV against an experimental 9.24 eV", Q, "-E_HOMO", "eV", calculation_function=calc.calc_ionization_potential),
            DescriptorInfo("ElectronAffinity", "Koopmans electron affinity, EA = -E(LUMO)", Q, "-E_LUMO", "eV", calculation_function=calc.calc_electron_affinity),
            DescriptorInfo("MullikenElectronegativity", "Absolute electronegativity chi = (IP + EA)/2; the driving force for charge transfer between two molecules", Q, "(IP + EA)/2", "eV", calculation_function=calc.calc_mulliken_electronegativity),
            DescriptorInfo("ElectrophilicityIndex", "Parr electrophilicity omega = mu²/(2*eta), the energy stabilisation on saturating the molecule with electrons; a standard reactivity ranking for Michael acceptors and electrophilic toxicants", Q, "mu² / (2*eta)", "eV", calculation_function=calc.calc_electrophilicity_index),
            DescriptorInfo("NucleophilicityIndex", "Reciprocal of the electrophilicity index, ranking electron donors", Q, "1 / omega", "1/eV", calculation_function=calc.calc_nucleophilicity_index),
            DescriptorInfo("ElectroacceptingPower", "Gazquez electron accepting power omega+ = (IP + 3*EA)²/(16*(IP - EA)), the propensity to accept charge", Q, "(IP+3EA)² / (16(IP-EA))", "eV", calculation_function=calc.calc_electroaccepting_power),
            DescriptorInfo("ElectrodonatingPower", "Gazquez electron donating power omega- = (3*IP + EA)²/(16*(IP - EA))", Q, "(3IP+EA)² / (16(IP-EA))", "eV", calculation_function=calc.calc_electrodonating_power),
            DescriptorInfo("PiElectronCount", "Number of electrons in the pi system, counting lone pairs donated by heteroatoms in aromatic rings", Q, "count(pi electrons)", "count", calculation_function=calc.calc_pi_electron_count),
            DescriptorInfo("PiSystemSize", "Number of atoms contributing a p orbital to the pi system", Q, "count(pi centres)", "count", calculation_function=calc.calc_pi_system_size),
            DescriptorInfo("PiSystemCount", "Number of separate conjugated fragments in the molecule", Q, "count(conjugated components)", "count", calculation_function=calc.calc_pi_system_count),
            DescriptorInfo("TotalPiEnergy", "Total Huckel pi electron energy in units of beta (8.000 for benzene); larger values mean more pi bonding", Q, "sum over occupied(2*x_i)", "beta", calculation_function=calc.calc_total_pi_energy),
            DescriptorInfo("DelocalizationEnergy", "Resonance energy in units of beta: the extra pi stabilisation over the same number of isolated double bonds. Benzene 2.000, butadiene 0.472, naphthalene 3.683, matching the textbook Huckel values", Q, "E_pi - 2*n_pairs", "beta", calculation_function=calc.calc_delocalization_energy),
            DescriptorInfo("GapPerPiAtom", "HOMO-LUMO gap divided by the size of the pi system; tracks how far conjugation extends rather than how large the molecule is", Q, "gap / n_pi", "eV", calculation_function=calc.calc_homo_lumo_gap_per_atom),
            DescriptorInfo("HardnessRatio", "Chemical hardness over electronegativity, a scale free reactivity ratio", Q, "eta / chi", "ratio", calculation_function=calc.calc_absolute_hardness_ratio),
        ]

    def _get_fingerprint_descriptors(self) -> List[DescriptorInfo]:
        """Get fingerprint descriptors."""
        calc = self.fingerprint_calc
        return [
            DescriptorInfo("MorganFingerprint", "Morgan (extended connectivity) fingerprint bit count, describing circular atom environments", DescriptorCategory.FINGERPRINTS, "circular_fingerprint", "bits", calculation_function=calc.calc_morgan_fingerprint),
            DescriptorInfo("MACCSFingerprint", "MACCS structural key fingerprint: 166 predefined substructure keys", DescriptorCategory.FINGERPRINTS, "structural_keys", "bits", calculation_function=calc.calc_maccs_fingerprint),
            DescriptorInfo("TopologicalFingerprint", "Path-based (Daylight-style) topological fingerprint bit count", DescriptorCategory.FINGERPRINTS, "topological_features", "bits", calculation_function=calc.calc_topological_fingerprint),
            DescriptorInfo("AtomPairFingerprint", "Atom pair fingerprint: hashed pairs of atom types with their topological separation", DescriptorCategory.FINGERPRINTS, "atom_pairs", "bits", calculation_function=calc.calc_atom_pair_fingerprint),
        ]

    def _get_hybrid_descriptors(self) -> List[DescriptorInfo]:
        """Get hybrid/drug-like descriptors."""
        calc = self.hybrid_calc
        return [
            DescriptorInfo("Lipophilicity (LogP)", "Octanol-water partition coefficient from the molecular lipophilicity potential (MLP) integrated over the solvent accessible surface; the primary lipophilicity descriptor", DescriptorCategory.HYBRID, "fragment_based_logP", "logP", calculation_function=calc.calc_lipophilicity),
            DescriptorInfo("Polarizability", "Molecular polarizability as the sum of atomic polarizability increments; governs dispersion interactions and refractive index", DescriptorCategory.HYBRID, "sum(atomic_polarizabilities)", "Å³", calculation_function=calc.calc_polarizability),
            DescriptorInfo("MolarRefractivity", "Molar refractivity from atomic increments, a combined size and polarizability term used by the Ghose filter", DescriptorCategory.HYBRID, "fragment_based_calculation", "cm³/mol", calculation_function=calc.calc_molar_refractivity),
            DescriptorInfo("LipinskiHBA", "Lipinski hydrogen bond acceptor count (nitrogen plus oxygen atoms), the counting rule used by the rule of five", DescriptorCategory.HYBRID, "lipinski_hba", "count", calculation_function=calc.calc_lipinski_hba),
            DescriptorInfo("LipinskiHBD", "Lipinski hydrogen bond donor count (N-H and O-H groups)", DescriptorCategory.HYBRID, "lipinski_hbd", "count", calculation_function=calc.calc_lipinski_hbd),
            DescriptorInfo("LipinskiViolationCount", "Number of Lipinski rule of five criteria violated: MW > 500, logP > 5, HBA > 10, HBD > 5. Two or more violations suggest poor oral absorption", DescriptorCategory.HYBRID, "lipinski_violations", "count", calculation_function=calc.calc_lipinski_violations),
            DescriptorInfo("VeberTPSA", "Polar surface area estimated from per-element increments, the form used in the original Veber analysis (see TPSA for the Ertl fragment version)", DescriptorCategory.HYBRID, "veber_psa", "Å²", calculation_function=calc.calc_veber_tpsa),
            DescriptorInfo("DrugLikenessScore", "Composite drug-likeness score on a 0-1 scale built from Lipinski violations, molecular weight, logP and rotatable bond count", DescriptorCategory.HYBRID, "drug_score", "score", calculation_function=calc.calc_drug_likeness_score),
            DescriptorInfo("SyntheticAccessibility", "Synthetic accessibility estimate on a 1-10 scale (lower is easier), penalising ring systems, flexibility and heteroatom richness", DescriptorCategory.HYBRID, "SA_score", "score", calculation_function=calc.calc_synthetic_accessibility),
            DescriptorInfo("Fsp3", "Fraction of sp3 hybridised carbons; a measure of three-dimensional character that correlates with clinical success", DescriptorCategory.HYBRID, "sp3/total_C", "fraction", calculation_function=calc.calc_fsp3),
        ] + self._get_hybrid_extended()

    def _get_hybrid_extended(self) -> List[DescriptorInfo]:
        """Polar surface area, solubility, volume and drug-likeness filters."""
        calc = self.hybrid_calc
        H = DescriptorCategory.HYBRID
        return [
            DescriptorInfo("TPSA", "Topological polar surface area from Ertl's fragment contributions for nitrogen and oxygen. Reproduces the published values exactly (aspirin 63.60, caffeine 58.44 Å²). Below ~90 Å² is associated with brain penetration and below ~140 Å² with good oral absorption", H, "sum(Ertl fragment contributions)", "Å²", calculation_function=calc.calc_tpsa),
            DescriptorInfo("TPSAWithSP", "Topological polar surface area counting sulfur and phosphorus as polar as well; the variant preferred for organophosphates and sulfur rich drugs", H, "sum(contributions incl. S,P)", "Å²", calculation_function=calc.calc_tpsa_with_sp),
            DescriptorInfo("TPSAPerHeavyAtom", "TPSA divided by the heavy atom count, comparable across molecular sizes", H, "TPSA / n_heavy", "Å²", calculation_function=calc.calc_tpsa_per_heavy_atom),
            DescriptorInfo("ESOLLogS", "Delaney ESOL estimate of aqueous solubility as log10 of the molar solubility: LogS = 0.16 - 0.63*logP - 0.0062*MW + 0.066*RotB - 0.74*AromaticProportion. Values below -6 indicate practically insoluble compounds", H, "0.16 - 0.63logP - 0.0062MW + 0.066RB - 0.74AP", "log(mol/L)", calculation_function=calc.calc_esol_logs),
            DescriptorInfo("McGowanVolume", "McGowan characteristic volume, an additive volume built from atomic increments minus 6.56 per bond. It is the size term of Abraham's solvation equations", H, "sum(V_atoms) - 6.56*n_bonds", "cm³/mol/100", calculation_function=calc.calc_mcgowan_volume),
            DescriptorInfo("HydrophilicFactor", "Todeschini hydrophilic factor, built from the number of hydrophilic groups (OH, NH, SH) relative to carbon and heavy atom counts. Strongly negative for pure hydrocarbons and large for polyols and sugars", H, "Todeschini HF", "index", calculation_function=calc.calc_hydrophilic_factor),
            DescriptorInfo("MolecularFlexibility", "Rotatable bonds per heavy atom; a size independent conformational freedom measure used in bioavailability models", H, "n_rot / n_heavy", "ratio", calculation_function=calc.calc_molecular_flexibility),
            DescriptorInfo("LogPPerHeavyAtom", "LogP divided by heavy atom count, the lipophilicity efficiency scale used to compare fragments with leads", H, "logP / n_heavy", "ratio", calculation_function=calc.calc_ligand_efficiency_scale),
            DescriptorInfo("PolarAtomFraction", "N, O, S and P atoms as a fraction of the heavy atoms", H, "n_polar / n_heavy", "fraction", calculation_function=calc.calc_polar_atom_fraction),
            DescriptorInfo("GhoseViolations", "Number of Ghose filter criteria violated: MW 160-480, logP -0.4 to 5.6, 20-70 heavy atoms and molar refractivity 40-130. Zero means the compound sits inside the drug-like qualifying range", H, "count(violations of 4 criteria)", "count", calculation_function=calc.calc_ghose_violations),
            DescriptorInfo("VeberViolations", "Veber oral bioavailability criteria violated: rotatable bonds <= 10 and TPSA <= 140 Å². Derived from rat oral bioavailability data", H, "count(violations of 2 criteria)", "count", calculation_function=calc.calc_veber_violations),
            DescriptorInfo("EganViolations", "Egan 'absorption egg' criteria violated: TPSA <= 131.6 Å² and logP <= 5.88, a two-property model of passive intestinal absorption", H, "count(violations of 2 criteria)", "count", calculation_function=calc.calc_egan_violations),
            DescriptorInfo("MueggeViolations", "Muegge pharmacophore filter criteria violated (MW, logP, TPSA, rings, carbon and heteroatom counts, rotatable bonds, HBA, HBD). A stricter drug-likeness screen than Lipinski", H, "count(violations of 9 criteria)", "count", calculation_function=calc.calc_muegge_violations),
            DescriptorInfo("RuleOfThreeViolations", "Congreve rule of three criteria violated: MW <= 300, logP <= 3, HBD <= 3, HBA <= 3, rotatable bonds <= 3. The standard filter for fragment screening libraries", H, "count(violations of 5 criteria)", "count", calculation_function=calc.calc_rule_of_three_violations),
            DescriptorInfo("LeadLikenessViolations", "Lead-likeness criteria violated: MW 250-350, logP <= 3.5, rotatable bonds <= 7. Leads need room to grow during optimisation", H, "count(violations of 3 criteria)", "count", calculation_function=calc.calc_lead_likeness_violations),
            DescriptorInfo("BBBScore", "Blood-brain barrier likelihood on a 0-1 scale, combining TPSA, logP, molecular weight, hydrogen bond donors and flexibility. Above ~0.7 is typical of CNS active drugs", H, "weighted property score", "score", calculation_function=calc.calc_bbb_score),
            DescriptorInfo("OralBioavailabilityScore", "Fraction of the six standard oral filters (Lipinski, Ghose, Veber, Egan, Muegge, lead-likeness) the molecule passes cleanly", H, "n_passed / 6", "fraction", calculation_function=calc.calc_oral_bioavailability_score),
        ]

    def _get_custom_descriptors(self) -> List[DescriptorInfo]:
        """Get custom descriptors."""
        calc = self.constitutional_calc
        return [
            DescriptorInfo("SelectionSize", "Number of atoms in the current selection, useful when descriptors are computed for a fragment", DescriptorCategory.CUSTOM, "len(selection)", "count", calculation_function=calc.calc_atom_count),
            DescriptorInfo("SelectionDensity", "Atoms of the selection per unit molecular volume, a packing density measure", DescriptorCategory.CUSTOM, "selection_size / molecular_volume", "atoms/Å³", calculation_function=self._calc_selection_density),
        ]

    def _calc_selection_density(self, molecule, selection) -> float:
        """Calculate selection density."""
        volume = self.geometric_calc.calc_molecular_volume(molecule, selection)
        if volume > 0:
            return len(selection.atom_indices) / volume
        return 0.0

    def calculate_descriptor(self, molecule: Molecule, descriptor_name: str,
                           selection: Optional[AtomSelection] = None) -> DescriptorResult:
        """Calculate a single descriptor value."""
        if selection is None:
            selection = AtomSelection(
                selection_type=SelectionType.ALL,
                atom_indices=list(range(len(molecule.atoms)))
            )

        for category, descriptors in self.descriptors.items():
            for desc in descriptors:
                if desc.name == descriptor_name:
                    if self.cache:
                        cache_key = self.cache.get_cache_key(molecule, selection, descriptor_name)
                        cached_value = self.cache.get(cache_key)
                        if cached_value is not None:
                            return DescriptorResult(
                                descriptor_name=descriptor_name,
                                value=cached_value,
                            )

                    try:
                        value = desc.calculation_function(molecule, selection)
                        if self.cache:
                            self.cache.set(cache_key, value)
                        return DescriptorResult(
                            descriptor_name=descriptor_name,
                            value=value,
                        )
                    except Exception as e:
                        import traceback
                        print(f"[DESCRIPTOR ERROR] {descriptor_name}: {e}")
                        traceback.print_exc()
                        return DescriptorResult(
                            descriptor_name=descriptor_name,
                            value=None,
                        )

        return DescriptorResult(
            descriptor_name=descriptor_name,
            value=None,
        )

    def calculate_all(self, molecule: Molecule,
                   selection: Optional[AtomSelection] = None,
                   categories: Optional[List[DescriptorCategory]] = None,
                   progress_callback: Optional[Callable[[CalculationProgress], None]] = None) -> Dict[str, DescriptorResult]:
        """Calculate all descriptors for a molecule."""
        if selection is None:
            selection = AtomSelection(
                selection_type=SelectionType.ALL,
                atom_indices=list(range(len(molecule.atoms)))
            )

        results = {}
        categories_to_calc = categories or list(DescriptorCategory)

        total = sum(len(self.descriptors.get(cat, [])) for cat in categories_to_calc)
        completed = 0

        for category in categories_to_calc:
            for desc in self.descriptors.get(category, []):
                result = self.calculate_descriptor(molecule, desc.name, selection)
                results[desc.name] = result

                completed += 1
                if progress_callback:
                    progress = CalculationProgress(
                        current_descriptor=desc.name,
                        current_category=category.value,
                        completed=completed,
                        total=total,
                        percentage=100 * completed / total,
                    )
                    progress_callback(progress)

        return results

    def get_descriptor_info(self, descriptor_name: str) -> Optional[DescriptorInfo]:
        """Get information about a specific descriptor."""
        for category, descriptors in self.descriptors.items():
            for desc in descriptors:
                if desc.name == descriptor_name:
                    return desc
        return None

    def list_descriptors(self, category: Optional[DescriptorCategory] = None) -> List[str]:
        """List available descriptors."""
        if category:
            return [desc.name for desc in self.descriptors.get(category, [])]
        return [desc.name for descriptors in self.descriptors.values() for desc in descriptors]

    def clear_cache(self):
        """Clear the calculation cache."""
        if self.cache:
            self.cache.clear()

    def set_progress_callback(self, callback: Optional[Callable[[CalculationProgress], None]]):
        """Set a progress callback for calculations."""
        self.progress_callback = callback

    def calculate_descriptors(self, molecule: Molecule, selection: AtomSelection,
                             categories: List[DescriptorCategory]) -> Dict[str, DescriptorResult]:
        """Calculate descriptors for given categories (GUI compatibility method)."""
        return self.calculate_all(molecule, selection, categories, self.progress_callback)
