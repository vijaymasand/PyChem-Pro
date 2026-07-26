# PyChem-Pro Molecular Descriptors

Reference for every descriptor produced by the **Molecular Descriptor Calculator** (`Tools -> Molecular Descriptor Calculator`).

The descriptor engine registers **318 descriptors across 8 categories**. The calculator additionally emits the systematically generated PyDes vectors described in the appendix, which brings a typical CSV export to roughly 700 columns per molecule.

Everything is computed with numpy and the project's own chemistry code; there are no external cheminformatics dependencies.

## Contents

- [Constitutional](#constitutional) - 105 descriptors
- [Topological](#topological) - 97 descriptors
- [Geometric](#geometric) - 31 descriptors
- [Electronic](#electronic) - 32 descriptors
- [Quantum](#quantum) - 20 descriptors
- [Fingerprints](#fingerprints) - 4 descriptors
- [Hybrid](#hybrid) - 27 descriptors
- [Custom](#custom) - 2 descriptors
- [Appendix: the PyDes vector families](#appendix-the-pydes-vector-families)

---

## Constitutional

Counts and compositional ratios read straight off the molecular graph: elements, bonds, rings, functional groups and atom environments. They need no geometry and no charges, so they are always available, and they are the descriptors that are easiest to interpret when a model turns out to depend on them.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `MolecularWeight` | g/mol | `sum(m_i) + n_H * 1.008` | Average molecular weight from the standard atomic weights, including the hydrogens a structure from SMILES only carries implicitly (aspirin 180.16 g/mol) |
| `AtomCount` | count | `len(atoms)` | Total number of atoms |
| `HeavyAtomCount` | count | `count(non-H)` | Number of non-hydrogen atoms |
| `BondCount` | count | `len(bonds)` | Total number of bonds |
| `RotatableBonds` | count | `count(rotatable)` | Number of rotatable bonds: acyclic single bonds between heavy atoms, excluding terminal ones. The standard conformational flexibility measure used by the Veber rules |
| `RingCount` | count | `count(rings)` | Number of rings in the smallest set of smallest rings (SSSR); the ring count of the molecular skeleton |
| `AromaticRingCount` | count | `count(aromatic_rings)` | Number of aromatic rings |
| `HeteroRingCount` | count | `count(hetero_rings)` | Number of heteroatom-containing rings |
| `Ring5Count` | count | `count(5-rings)` | Number of 5-membered rings |
| `Ring6Count` | count | `count(6-rings)` | Number of 6-membered rings |
| `LargestRingSize` | count | `max(ring_sizes)` | Size of the largest ring |
| `HDonorCount` | count | `count(H-donors)` | Number of hydrogen bond donors: N-H and O-H groups |
| `HAcceptorCount` | count | `count(H-acceptors)` | Number of hydrogen bond acceptors: nitrogen and oxygen atoms with a free lone pair |
| `LipophilicAtomCount` | count | `count(lipophilic)` | Number of atoms classified as lipophilic (carbons and halogens with no polar neighbour) |
| `CarbonCount` | count | `count(C)` | Number of carbon atoms |
| `NitrogenCount` | count | `count(N)` | Number of nitrogen atoms |
| `OxygenCount` | count | `count(O)` | Number of oxygen atoms |
| `SulfurCount` | count | `count(S)` | Number of sulfur atoms |
| `PhosphorusCount` | count | `count(P)` | Number of phosphorus atoms |
| `HalogenCount` | count | `count(F,Cl,Br,I)` | Number of halogen atoms |
| `FluorineCount` | count | `count(F)` | Number of fluorine atoms |
| `ChlorineCount` | count | `count(Cl)` | Number of chlorine atoms |
| `BromineCount` | count | `count(Br)` | Number of bromine atoms |
| `IodineCount` | count | `count(I)` | Number of iodine atoms |
| `HydrogenCount` | count | `count(H)` | Number of hydrogen atoms |
| `SingleBondCount` | count | `count(single_bonds)` | Number of single bonds |
| `DoubleBondCount` | count | `count(double_bonds)` | Number of double bonds |
| `TripleBondCount` | count | `count(triple_bonds)` | Number of triple bonds |
| `AromaticBondCount` | count | `count(aromatic_bonds)` | Number of aromatic bonds |
| `HeteroAtomCount` | count | `count(hetero)` | Number of hetero atoms |
| `AliphaticCarbonCount` | count | `count(aliphatic_C)` | Number of aliphatic carbons |
| `AromaticCarbonCount` | count | `count(aromatic_C)` | Number of aromatic carbons |
| `NonRingCarbonCount` | count | `count(C_non-ring)` | Number of non-ring carbon atoms |
| `SP3CarbonCount` | count | `count(sp3_C)` | Number of sp3 hybridized carbon atoms |
| `SP2CarbonCount` | count | `count(sp2_C)` | Number of sp2 hybridized carbon atoms |
| `SPCarbonCount` | count | `count(sp_C)` | Number of sp hybridized carbon atoms |
| `SP3NitrogenCount` | count | `count(sp3_N)` | Number of sp3 hybridized nitrogen atoms |
| `SP2NitrogenCount` | count | `count(sp2_N)` | Number of sp2 hybridized nitrogen atoms |
| `SP3OxygenCount` | count | `count(sp3_O)` | Number of sp3 hybridized oxygen atoms |
| `SP2OxygenCount` | count | `count(sp2_O)` | Number of sp2 hybridized oxygen atoms |
| `FormalCharge` | e | `sum(formal_charges)` | Sum of the formal charges of the selected atoms, the net charge of the species |
| `UnsaturationCount` | count | `DBE` | Degree of unsaturation (double bond equivalents): rings plus pi bonds, computed from the molecular formula |
| `AromaticProportion` | fraction | `aromatic/total` | Fraction of atoms belonging to an aromatic ring |
| `AliphaticProportion` | fraction | `aliphatic/total` | Fraction of atoms that are not aromatic |
| `SP3Proportion` | fraction | `sp3_C/total_C` | Fraction of carbons that are sp3 hybridised |
| `HeteroProportion` | fraction | `hetero/total` | Fraction of atoms that are neither carbon nor hydrogen |
| `HCRatio` | ratio | `H/C` | Hydrogen to carbon ratio |
| `NCRatio` | ratio | `N/C` | Nitrogen to carbon ratio |
| `OCRatio` | ratio | `O/C` | Oxygen to carbon ratio |
| `AmideBondCount` | count | `count(amide_bonds)` | Number of amide bonds (C(=O)-N), the linkage of peptides and a very common medicinal chemistry motif |
| `EsterBondCount` | count | `count(ester_bonds)` | Number of ester linkages (C(=O)-O-C), frequent metabolic soft spots |
| `CarbonylBondCount` | count | `count(carbonyl_bonds)` | Number of C=O double bonds of any kind (aldehyde, ketone, acid, ester, amide) |
| `HydroxylGroupCount` | count | `count(OH_groups)` | Number of hydroxyl (-OH) groups |
| `CarboxylGroupCount` | count | `count(COOH_groups)` | Number of carboxyl (-COOH) groups |
| `AmineGroupCount` | count | `count(amine_groups)` | Number of nitrogen atoms that are not part of an amide, i.e. amine-like nitrogens |
| `MethylGroupCount` | count | `count(methyl_groups)` | Number of methyl (-CH3) groups |
| `LongestChain` | count | `max(chain_length)` | Number of atoms in the longest carbon chain, the classic measure of molecular length |
| `BranchCount` | count | `count(branches)` | Number of atoms carrying three or more heavy neighbours, i.e. branch points |
| `FragmentCount` | count | `count(fragments)` | Number of disconnected fragments |
| `ExactMass` | g/mol | `sum(m_i) + n_H * 1.008` | Mass computed from the element masses with implicit hydrogens added; differs from MolecularWeight when the structure carries implicit H |
| `HeavyAtomMolWeight` | g/mol | `sum(m_i, i != H)` | Molecular weight of the non-hydrogen atoms only; useful for normalising properties by molecular size |
| `AverageAtomicMass` | g/mol | `MW_heavy / n_heavy` | Mean atomic mass per heavy atom; rises with heavy halogens and metals, a proxy for how 'heavy' the elements are |
| `ValenceElectronCount` | count | `sum(Zv_i)` | Total number of valence electrons; correlates strongly with polarizability and dispersion interactions |
| `VdWVolumeSum` | Å³ | `sum(V_vdW,i)` | Sum of Bondi van der Waals atomic volumes, without any packing correction; a coordinate-free size measure |
| `PrimaryCarbonCount` | count | `count(C with 1 C neighbour)` | Carbons bonded to exactly one other carbon (methyl and chain ends) |
| `SecondaryCarbonCount` | count | `count(C with 2 C neighbours)` | Carbons bonded to two other carbons (chain interior) |
| `TertiaryCarbonCount` | count | `count(C with 3 C neighbours)` | Carbons bonded to three other carbons (branch points) |
| `QuaternaryCarbonCount` | count | `count(C with 4 C neighbours)` | Carbons bonded to four other carbons; sterically congested centres that slow metabolism |
| `TerminalAtomCount` | count | `count(deg = 1)` | Heavy atoms with a single heavy neighbour, i.e. the leaves of the molecular graph |
| `MaxAtomDegree` | count | `max(deg_i)` | Largest number of heavy neighbours on any atom; measures the most branched centre |
| `MeanAtomDegree` | count | `2 * E / V` | Average heavy atom degree, equal to 2*bonds/atoms; grows with ring fusion and branching |
| `RingAtomCount` | count | `\|union(rings)\|` | Number of atoms that belong to at least one ring |
| `RingBondCount` | count | `count(ring bonds)` | Number of bonds that belong to at least one ring |
| `AliphaticRingCount` | count | `count(non-aromatic rings)` | Rings containing at least one non-aromatic atom |
| `SaturatedRingCount` | count | `count(all-single rings)` | Rings in which every bond is single; a marker of three-dimensional character |
| `AromaticHeterocycleCount` | count | `count(aromatic hetero rings)` | Aromatic rings that contain at least one heteroatom (pyridine, imidazole, ...) |
| `AromaticCarbocycleCount` | count | `count(aromatic all-C rings)` | Aromatic rings made only of carbon (benzene rings) |
| `SaturatedHeterocycleCount` | count | `count(saturated hetero rings)` | Fully saturated rings containing a heteroatom (piperidine, morpholine, ...) |
| `SpiroAtomCount` | count | `count(\|R_i ∩ R_j\| = 1)` | Atoms shared by two rings that have no bond in common; spiro centres enforce perpendicular ring planes |
| `BridgeheadAtomCount` | count | `count(shared ring junctions)` | Atoms at the junction of two rings sharing a bond or a bridge; characteristic of rigid cage systems |
| `FusedRingCount` | count | `count(\|R_i ∩ R_j\| >= 2)` | Number of ring pairs that share at least one bond, i.e. the number of fusion sites |
| `MacrocycleCount` | count | `count(ring size > 12)` | Rings larger than 12 atoms, which behave very differently from ordinary rings in ADME terms |
| `RingComplexity` | ratio | `\|union(rings)\| / n_rings` | Ring atoms per ring; approaches the ring size for isolated rings and drops as rings become fused |
| `NitroGroupCount` | count | `count(N with 2 terminal O)` | Number of nitro groups (R-NO2), recognised in both the charge separated and the pentavalent notation |
| `NitrileCount` | count | `count(C#N)` | Number of nitrile (C#N) groups |
| `AldehydeCount` | count | `count(HC=O)` | Number of aldehyde groups: a carbonyl carbon carrying a hydrogen |
| `KetoneCount` | count | `count(C-CO-C)` | Number of ketone groups: a carbonyl carbon flanked by two carbons |
| `EtherCount` | count | `count(C-O-C)` | Number of ether oxygens (C-O-C), excluding esters and acids |
| `ThioetherCount` | count | `count(C-S-C)` | Number of thioether sulfurs (C-S-C) |
| `ThiolCount` | count | `count(SH)` | Number of thiol (-SH) groups; frequent covalent-binding and metabolic liabilities |
| `SulfonamideCount` | count | `count(SO2N)` | Number of sulfonamide groups S(=O)(=O)N |
| `SulfoneCount` | count | `count(C-SO2-C)` | Number of sulfone groups C-S(=O)(=O)-C |
| `SulfoxideCount` | count | `count(C-SO-C)` | Number of sulfoxide groups C-S(=O)-C |
| `PhosphateCount` | count | `count(P with >= 3 O)` | Number of phosphorus atoms carrying three or more oxygens (phosphate/phosphonate) |
| `ImineCount` | count | `count(C=N)` | Number of non-aromatic C=N double bonds (imines and Schiff bases) |
| `AzoCount` | count | `count(N=N)` | Number of N=N double bonds (azo groups, common in dyes) |
| `UreaCount` | count | `count(N-CO-N)` | Number of urea/carbamate-like fragments N-C(=O)-N |
| `GuanidineCount` | count | `count(C(N)(N)=N)` | Number of guanidine groups; strongly basic and permanently charged at physiological pH |
| `PhenolCount` | count | `count(Ar-OH)` | Number of hydroxyls attached to an aromatic ring; more acidic than aliphatic alcohols |
| `AliphaticHydroxylCount` | count | `count(alkyl-OH)` | Number of alcohol hydroxyls that are neither phenolic nor part of a carboxylic acid |
| `AromaticHalideCount` | count | `count(Ar-X)` | Number of halogens attached directly to an aromatic ring |
| `AcidicGroupCount` | count | `count(COOH + SO3H + PO3H)` | Ionisable acidic groups: carboxylic, sulfonic and phosphonic acids |
| `BasicGroupCount` | count | `count(basic N)` | Basic nitrogens: aliphatic amines, amidines and guanidines, excluding amides |
| `RotatableBondFraction` | ratio | `n_rot / n_heavy` | Rotatable bonds per heavy atom, a size independent measure of conformational freedom |
| `HeteroatomFraction` | ratio | `n_hetero / n_heavy` | Heteroatoms per heavy atom; a compact polarity/complexity indicator |

---

## Topological

Graph invariants computed on the **hydrogen suppressed** molecular graph, as the original definitions require. They encode size, branching, shape and heteroatom placement without any 3D information, which makes them reproducible and conformation independent - the reason they dominate classical QSAR.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `WienerIndex` | index | `sum(distance_matrix)` | Wiener index: the sum of topological distances over all atom pairs. The oldest topological index, it grows with size and decreases with branching, and correlates with boiling point and molar volume |
| `ZagrebIndex1` | index | `sum(degree^2)` | First Zagreb index, the sum of squared vertex degrees; a branching measure that grows sharply with substitution at a single centre |
| `ZagrebIndex2` | index | `sum(degree_i * degree_j)` | Second Zagreb index, the sum of degree products over bonds; complements the first Zagreb index by weighting bonds rather than atoms |
| `BalabanIndex` | index | `J_index` | Balaban J index, an average distance-sum connectivity that stays almost constant with molecular size, so it discriminates isomers and ring systems rather than size |
| `HararyIndex` | index | `sum(1/distance)` | Harary index, the sum of reciprocal topological distances; unlike the Wiener index it emphasises close atom pairs and therefore local compactness |
| `RandicIndex` | index | `chi_index` | Randic branching index, the sum of 1/sqrt(d_i*d_j) over bonds; the original molecular connectivity descriptor, lower for more branched skeletons |
| `ConnectivityIndexChi0` | index | `chi_0` | Zeroth order connectivity index computed over bonds (legacy variant kept for backward compatibility; see Chi0n for the standard Kier-Hall definition) |
| `ConnectivityIndexChi1` | index | `chi_1` | First order connectivity index (legacy variant kept for backward compatibility; see Chi1n for the standard Kier-Hall definition) |
| `KappaShapeIndex1` | index | `kappa_1` | First order Kappa shape index without the alpha correction, comparing the bond count with the linear and complete graph limits (see Kappa1Alpha for the size corrected form) |
| `KappaShapeIndex2` | index | `kappa_2` | Second order Kappa shape index without the alpha correction; grows with molecular elongation |
| `KappaShapeIndex3` | index | `kappa_3` | Third order Kappa shape index without the alpha correction; reflects where the branch points sit along the skeleton |
| `HosoyaIndex` | index | `Z_index` | Hosoya Z index, a count of the independent bond matchings of the graph; a classical complexity measure related to thermodynamic stability |
| `PlattIndex` | index | `sum(edge_degrees)` | Platt index, the sum of the degrees of the two atoms of every bond; a simple edge-based branching measure correlating with molar volume |
| `PolarityNumber` | count | `count(3-dist_pairs)` | Polarity number (Wiener polarity): the count of atom pairs exactly three bonds apart, which is the number of rotatable-bond-like arrangements |
| `BertzIndex` | index | `complexity_based_index` | Bertz molecular complexity index, an information measure of the bonding pattern; grows with both size and heterogeneity of the skeleton |
| `BonchevTrinajsticIndex` | index | `information_based_index` | Bonchev-Trinajstic mean information content on the vertex degree distribution, in bits per atom |
| `InformationContent0` | bits | `entropy_based_index` | Zeroth order information content: the entropy of the atom equivalence classes defined by degree alone |
| `InformationContent1` | bits | `entropy_based_index` | First order information content: entropy of the atom classes after one sphere of neighbours is taken into account |
| `EccentricConnectivityIndex` | index | `sum(degree*eccentricity)` | Eccentric connectivity index, the sum over atoms of degree times eccentricity; combines branching with molecular elongation and is widely used in QSAR of drug transport |
| `PathCount3` | count | `count(3-paths)` | Number of distinct paths spanning three atoms (two bonds); a local shape descriptor |
| `PathCount4` | count | `count(4-paths)` | Number of distinct paths spanning four atoms (three bonds) |
| `AveragePathLength` | index | `avg_path_length` | Mean topological distance between connected atom pairs; small for compact fused systems and large for extended chains |
| `LongestShortestPath` | count | `graph_diameter` | Graph diameter: the longest of all shortest paths, i.e. the topological length of the molecule |
| `CyclomaticNumber` | count | `E - V + C` | Cyclomatic number (circuit rank) E - V + C: the number of independent rings, i.e. how many bonds must be cut to make the molecule acyclic |
| `ElectrotopologicalStateIndex` | index | `electronic_topology_index` | Sum of intrinsic electrotopological states (atomic number over degree); a combined electronic and topological measure of atom accessibility |
| `MolecularTopologicalIndex` | index | `topological_complexity` | Composite topological index combining the Wiener and first Zagreb indices into a single size-and-branching measure |
| `HyperWienerIndex` | index | `sum(distance + distance^2)` | Hyper-Wiener index, which adds squared distances to the Wiener sum; discriminates isomers that the Wiener index cannot separate |
| `Chi0n` | index | `sum(1/sqrt(d_i))` | Kier-Hall simple connectivity index of order 0: sums 1/sqrt(delta) over atoms, where delta is the heavy atom degree. Grows with size and with the number of unbranched atoms |
| `Chi1n` | index | `sum(1/sqrt(d_i*d_j)) over bonds` | First order simple connectivity index (Randic index) over bonds; the classic branching descriptor, lower for more branched skeletons |
| `Chi2n` | index | `sum over 2-paths` | Second order simple connectivity index over two-bond paths; encodes the local environment one bond further out |
| `Chi3n` | index | `sum over 3-paths` | Third order simple connectivity index over three-bond paths |
| `Chi4n` | index | `sum over 4-paths` | Fourth order simple connectivity index over four-bond paths; sensitive to medium range shape |
| `Chi0v` | index | `sum(1/sqrt(dv_i))` | Valence connectivity index of order 0; uses delta_v = (Zv - h)/(Z - Zv - 1), so heteroatoms and their lone pairs are distinguished from carbon |
| `Chi1v` | index | `sum(1/sqrt(dv_i*dv_j))` | First order valence connectivity index; the electronic counterpart of the Randic index and a strong predictor of logP and boiling point |
| `Chi2v` | index | `sum over 2-paths (valence)` | Second order valence connectivity index over two-bond paths |
| `Chi3v` | index | `sum over 3-paths (valence)` | Third order valence connectivity index over three-bond paths |
| `Chi4v` | index | `sum over 4-paths (valence)` | Fourth order valence connectivity index over four-bond paths |
| `Chi3Cluster` | index | `sum over 3-stars` | Third order cluster connectivity: sums over three-branch stars, so it counts tertiary branch points specifically |
| `Chi4Cluster` | index | `sum over 4-stars` | Fourth order cluster connectivity over four-branch stars, i.e. quaternary centres |
| `Kappa1Alpha` | index | `(A+a)(A+a-1)^2 / (P1+a)^2` | First Kier shape index with the alpha correction for atom size and hybridisation; compares the molecule with the linear and the fully connected limits |
| `Kappa2Alpha` | index | `(A+a-1)(A+a-2)^2 / (P2+a)^2` | Second Kier shape index; increases with molecular elongation and decreases with branching and cyclisation |
| `Kappa3Alpha` | index | `(A+a-1)(A+a-3)^2 / (P3+a)^2` | Third Kier shape index, sensitive to the position of branch points along the skeleton |
| `KierFlexibility` | index | `kappa1 * kappa2 / A` | Kier molecular flexibility index Phi = kappa1*kappa2/A; roughly 0 for rigid fused systems and large for long flexible chains |
| `ABCIndex` | index | `sum(sqrt((d_i+d_j-2)/(d_i*d_j)))` | Atom-bond connectivity index; correlates with the strain energy of branched alkanes |
| `AugmentedZagrebIndex` | index | `sum((d_i*d_j/(d_i+d_j-2))^3)` | Augmented Zagreb index, a heavily branch-weighted degree index used in heat-of-formation models |
| `ForgottenIndex` | index | `sum(d_i^3)` | Forgotten topological index, the sum of cubed vertex degrees; emphasises highly connected atoms |
| `ModifiedZagrebIndex` | index | `sum(1/d_i^2)` | Modified first Zagreb index, the sum of inverse squared degrees; weights terminal atoms most |
| `SumConnectivityIndex` | index | `sum(1/sqrt(d_i+d_j))` | Sum-connectivity index, a variant of the Randic index using the degree sum instead of the product |
| `GeometricArithmeticIndex` | index | `sum(2*sqrt(d_i*d_j)/(d_i+d_j))` | Geometric-arithmetic index; equals the bond count for regular graphs and drops as degrees become uneven |
| `HarmonicIndex` | index | `sum(2/(d_i+d_j))` | Harmonic index, the harmonic mean counterpart of the Randic index |
| `NarumiSimpleIndex` | index | `sum(ln d_i)` | Logarithm of the product of all vertex degrees (the product itself overflows for large molecules) |
| `NarumiHarmonicIndex` | index | `A / sum(1/d_i)` | Harmonic mean of the vertex degrees |
| `NarumiGeometricIndex` | index | `(prod d_i)^(1/A)` | Geometric mean of the vertex degrees |
| `SchultzIndex` | index | `sum(d_i * (A+D)_ij)` | Schultz molecular topological index, combining the adjacency and distance matrices weighted by degree |
| `GutmanIndex` | index | `sum(d_i*d_j*D_ij)` | Gutman (Schultz of the second kind) index: degree weighted sum of topological distances |
| `SzegedIndex` | index | `sum(n_u * n_v) over bonds` | Szeged index; for every bond it multiplies the number of atoms closer to each end, generalising the Wiener index to cyclic graphs |
| `PIIndex` | index | `sum(n_u + n_v) over bonds` | Padmakar-Ivan index, the additive counterpart of the Szeged index; discriminates ring systems well |
| `MeanWienerIndex` | bonds | `W / C(n,2)` | Wiener index divided by the number of atom pairs, i.e. the mean topological distance; a size independent compactness measure |
| `KirchhoffIndex` | index | `n * sum(1/mu_k)` | Resistance distance (quasi-Wiener) index computed from the Laplacian spectrum; unlike the Wiener index it accounts for every path between two atoms, not just the shortest one |
| `TopologicalRadius` | bonds | `min(ecc_i)` | Smallest atomic eccentricity, the graph radius |
| `AverageEccentricity` | bonds | `mean(ecc_i)` | Mean atomic eccentricity; grows linearly with chain length and slowly for compact ring systems |
| `PetitjeanIndex` | index | `(D - R) / R` | Topological shape index (diameter - radius)/radius; 0 for perfectly symmetric graphs such as benzene and approaching 1 for long chains |
| `EccentricDistanceSum` | index | `sum(ecc_i * sum_j D_ij)` | Sum over atoms of eccentricity times the total distance to all other atoms; combines shape and size |
| `GraphEnergy` | index | `sum(\|lambda_i\|)` | Graph energy, the sum of absolute adjacency eigenvalues; tracks the total pi bonding capacity of the skeleton |
| `SpectralRadius` | index | `max(\|lambda_i\|)` | Largest adjacency eigenvalue; bounded by the maximum degree and larger for densely fused systems |
| `EstradaIndex` | index | `sum(exp(lambda_i))` | Estrada index, the sum of exponentials of the adjacency eigenvalues; a folding/compactness measure dominated by short closed walks |
| `LaplacianSpectralRadius` | index | `max(mu_i)` | Largest Laplacian eigenvalue, an upper bound related to the maximum degree and the connectivity of the graph |
| `AlgebraicConnectivity` | index | `mu_2` | Fiedler value, the second smallest Laplacian eigenvalue; measures how hard the molecular graph is to cut in two and is 0 for disconnected structures |
| `LogSpanningTreeCount` | index | `log10(prod(mu_k)/n)` | Log10 of the number of spanning trees (Kirchhoff matrix-tree theorem); 0 for acyclic molecules and rising with ring count and fusion |
| `WalkCount2` | count | `sum(A^2)` | Number of walks of length 2 in the molecular graph |
| `WalkCount3` | count | `sum(A^3)` | Number of walks of length 3; together with the lower orders it forms a classic complexity series |
| `WalkCount4` | count | `sum(A^4)` | Number of walks of length 4 |
| `SelfReturningWalk3` | count | `tr(A^3)` | Trace of A^3, exactly six times the number of three-membered rings |
| `ATSElectronegativity1` | index | `sum(w_i*w_j), d_ij = 1` | Moreau-Broto autocorrelation of electronegativity at topological distance 1: sums the product of the (carbon-scaled) electronegativities of bonded atoms |
| `ATSElectronegativity2` | index | `sum(w_i*w_j), d_ij = 2` | Electronegativity autocorrelation at distance 2 (atoms separated by two bonds) |
| `ATSElectronegativity3` | index | `sum(w_i*w_j), d_ij = 3` | Electronegativity autocorrelation at distance 3 |
| `ATSPolarizability1` | index | `sum(w_i*w_j), d_ij = 1` | Autocorrelation of atomic polarizability at distance 1; encodes where the easily polarised atoms sit relative to each other |
| `ATSPolarizability2` | index | `sum(w_i*w_j), d_ij = 2` | Atomic polarizability autocorrelation at distance 2 |
| `ATSPolarizability3` | index | `sum(w_i*w_j), d_ij = 3` | Atomic polarizability autocorrelation at distance 3 |
| `ATSVolume1` | index | `sum(w_i*w_j), d_ij = 1` | Autocorrelation of van der Waals volume at distance 1; a steric analogue of the electronegativity autocorrelations |
| `ATSVolume2` | index | `sum(w_i*w_j), d_ij = 2` | Van der Waals volume autocorrelation at distance 2 |
| `ATSVolume3` | index | `sum(w_i*w_j), d_ij = 3` | Van der Waals volume autocorrelation at distance 3 |
| `MoranElectronegativity1` | index | `Moran I, d = 1` | Moran spatial autocorrelation of electronegativity at distance 1, normalised to [-1, 1]; positive means neighbouring atoms have similar electronegativity |
| `MoranElectronegativity2` | index | `Moran I, d = 2` | Moran autocorrelation of electronegativity at distance 2 |
| `MoranMass1` | index | `Moran I, d = 1` | Moran autocorrelation of atomic mass between bonded atoms |
| `GearyElectronegativity1` | index | `Geary C, d = 1` | Geary autocorrelation of electronegativity at distance 1; near 0 for smooth property distributions and above 1 when bonded atoms differ strongly |
| `GearyElectronegativity2` | index | `Geary C, d = 2` | Geary autocorrelation of electronegativity at distance 2 |
| `GearyMass1` | index | `Geary C, d = 1` | Geary autocorrelation of atomic mass between bonded atoms |
| `TopologicalCharge1` | index | `sum\|CT_ij\|, d_ij = 1` | Galvez topological charge index at distance 1; measures the net charge transfer implied by the topology between bonded atoms |
| `TopologicalCharge2` | index | `sum\|CT_ij\|, d_ij = 2` | Galvez topological charge index at distance 2 |
| `TopologicalCharge3` | index | `sum\|CT_ij\|, d_ij = 3` | Galvez topological charge index at distance 3 |
| `MeanTopologicalCharge` | index | `sum(GGI_k)/(n-1)` | Mean Galvez charge index over distances 1-3, normalised by the number of atoms minus one |
| `AtomTypeInformation` | bits | `-sum(p*log2 p)` | Shannon entropy of the element distribution, in bits per atom; 0 for a hydrocarbon skeleton and rising with elemental diversity |
| `TotalAtomInformation` | bits | `n * -sum(p*log2 p)` | Total information content of the element distribution (entropy times atom count) |
| `BondTypeInformation` | bits | `-sum(p*log2 p)` | Shannon entropy of the bond order distribution (single/double/triple/aromatic) |
| `DistanceInformation` | bits | `n_pairs * -sum(p*log2 p)` | Bonchev-Trinajstic information index on the distribution of topological distances; a well established measure of structural complexity |
| `VertexDegreeInformation` | bits | `-sum(p*log2 p)` | Shannon entropy of the vertex degree distribution, a compact branching-diversity measure |

---

## Geometric

Descriptors that need a 3D geometry: surface, volume, moments of inertia, shape ratios and extent along the principal axes. When a structure has no coordinates (for example one parsed from SMILES) PyChem generates them once with its own 3D generator, so these values are never silently zero.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `SASA` | Å² | `surface_area_calculation` | Solvent accessible surface area from the Shrake-Rupley algorithm with a 1.4 A probe; the surface a water molecule can touch |
| `MolecularVolume` | Å³ | `volume_calculation` | Van der Waals volume obtained by numerical integration over the atomic spheres |
| `RadiusOfGyration` | Å | `sqrt(sum((r - r_center)^2 / N)` | Root mean square distance of the atoms from the centroid; the compactness measure used for polymers and conformer comparison |
| `Asphericity` | dimensionless | `eigenvalue_based` | Asphericity from the gyration tensor eigenvalues: 0 for a spherically symmetric molecule and 1 for a perfectly linear one |
| `Eccentricity` | dimensionless | `eigenvalue_based` | Geometric eccentricity from the gyration tensor; 0 for a sphere and approaching 1 for elongated shapes |
| `PrincipalMoment1` | Å² | `eigenvalue1` | Largest eigenvalue of the (unweighted) gyration tensor, the spatial extent along the principal axis |
| `PrincipalMoment2` | Å² | `eigenvalue2` | Second eigenvalue of the gyration tensor |
| `PrincipalMoment3` | Å² | `eigenvalue3` | Smallest eigenvalue of the gyration tensor; near zero for planar molecules |
| `MolecularDiameter` | Å | `max_distance` | Largest through-space distance between any two atoms, the diameter of the enclosing sphere |
| `InertiaMomentA` | amu·Å² | `eigenvalue 1 of the inertia tensor` | Smallest principal moment of inertia (mass weighted); small for rod shaped molecules |
| `InertiaMomentB` | amu·Å² | `eigenvalue 2 of the inertia tensor` | Intermediate principal moment of inertia |
| `InertiaMomentC` | amu·Å² | `eigenvalue 3 of the inertia tensor` | Largest principal moment of inertia |
| `NPR1` | ratio | `I1 / I3` | Normalised principal moment ratio I1/I3. Together with NPR2 it places a molecule on the rod-disc-sphere triangle: rods sit near (0.5, 0.5), discs near (0, 1) and spheres near (1, 1) |
| `NPR2` | ratio | `I2 / I3` | Normalised principal moment ratio I2/I3, the second coordinate of the principal moments of inertia shape triangle |
| `InertialShapeFactor` | 1/(amu·Å²) | `I2 / (I1 * I3)` | I2/(I1*I3); large for elongated molecules and small for compact globular ones |
| `SpherocityIndex` | dimensionless | `3*L1 / (L1+L2+L3)` | Three times the smallest gyration eigenvalue over their sum; 1 for a perfect sphere and 0 for a planar or linear molecule |
| `MolecularLength` | Å | `max-min projection on axis 1` | Extent of the molecule along its longest principal axis |
| `MolecularWidth` | Å | `max-min projection on axis 2` | Extent along the second principal axis |
| `MolecularThickness` | Å | `max-min projection on axis 3` | Extent along the shortest principal axis; near 0 for flat aromatic systems |
| `LengthToWidthRatio` | ratio | `L / W` | Aspect ratio of the molecule, length divided by width |
| `BoundingBoxVolume` | Å³ | `L * W * T` | Volume of the box aligned with the principal axes that encloses the molecule |
| `Span` | Å | `max\|r_i - centroid\|` | Largest distance from the centroid to any atom, the radius of the enclosing sphere |
| `PlaneOfBestFit` | Å | `mean\|(r_i - r0)·n\|` | Mean distance of the atoms from their best fit plane; 0 for planar molecules and rising with three-dimensional character |
| `MassWeightedRadiusOfGyration` | Å | `sqrt(sum(m_i*r_i²)/sum(m_i))` | Radius of gyration weighted by atomic mass, the physical definition used in polymer and dynamics work |
| `MeanAtomicDistance3D` | Å | `mean(\|r_i - r_j\|)` | Mean through-space distance over all atom pairs |
| `GeometricWienerIndex` | Å | `sum(\|r_i - r_j\|)` | Sum of all through-space interatomic distances, the 3D analogue of the Wiener index |
| `GravitationalIndex` | amu²/Å² | `sum(m_i*m_j / r_ij²)` | Katritzky gravitational index over all atom pairs; encodes mass distribution and correlates with boiling point and density |
| `GravitationalIndexBonded` | amu²/Å² | `sum over bonds(m_i*m_j / r_ij²)` | Gravitational index restricted to bonded atom pairs, a local mass-distribution term |
| `SurfaceToVolumeRatio` | 1/Å | `SASA / V` | Solvent accessible surface divided by molecular volume; high for extended or branched shapes |
| `Globularity` | dimensionless | `(36*pi*V²)^(1/3) / SASA` | Surface of the sphere of equal volume divided by the real surface; 1 for a sphere and below 1 for anything rugged or elongated |
| `MolecularDensity` | g/mol/Å³ | `MW / V` | Molecular weight per unit molecular volume |

---

## Electronic

Partial charge statistics and the charged partial surface area (CPSA) family of Stanton and Jurs, which combine Gasteiger charges with the per-atom solvent accessible surface. Gasteiger charges are computed on demand when a molecule does not carry them.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `TotalCharge` | e | `sum(partial_charges)` | Sum of the Gasteiger partial charges over the selection; close to the formal charge for a whole molecule |
| `MaxPartialCharge` | e | `max(partial_charges)` | Most positive Gasteiger partial charge, marking the strongest electrophilic centre |
| `MinPartialCharge` | e | `min(partial_charges)` | Most negative Gasteiger partial charge, marking the strongest nucleophilic centre |
| `MaxPositiveCharge` | e | `max_positive` | Largest positive partial charge, a hydrogen bond donor strength indicator |
| `MaxNegativeCharge` | e | `max_negative` | Largest negative partial charge, a hydrogen bond acceptor strength indicator |
| `DipoleMoment` | Debye | `sqrt(sum(charge * position)^2)` | Molecular dipole moment computed from the Gasteiger charges and the 3D geometry; drives solubility and crystal packing |
| `PolarSurfaceArea` | Å² | `surface_area_of_polar_atoms` | Solvent accessible surface contributed by the polar atoms (N, O and their hydrogens), computed from the real 3D surface rather than from fragment values |
| `APolarSurfaceArea` | Å² | `surface_area_of_apolar_atoms` | Solvent accessible surface contributed by the non-polar atoms, the hydrophobic part of the surface |
| `MeanAbsoluteCharge` | e | `mean(\|charge\|)` | Mean absolute partial charge per atom, a compact measure of overall molecular polarity |
| `TotalAbsoluteCharge` | e | `sum(\|q_i\|)` | Sum of the absolute partial charges; the electronic 'size' of the molecule and a measure of overall charge separation |
| `ChargeVariance` | e² | `var(q_i)` | Variance of the atomic partial charges; large when the molecule mixes strongly polarised and neutral regions |
| `SumSquaredCharges` | e² | `sum(q_i²)` | Sum of squared partial charges, the electrostatic self-energy term of the charge distribution |
| `SubmolecularPolarity` | e | `max(q) - min(q)` | Largest difference between any two atomic charges (Todeschini submolecular polarity parameter); a direct measure of internal charge separation |
| `RelativePositiveCharge` | fraction | `q_max / sum(q > 0)` | Most positive atomic charge divided by the total positive charge; identifies whether positive charge is concentrated on one atom or spread out |
| `RelativeNegativeCharge` | fraction | `\|q_min\| / sum\|q < 0\|` | Most negative atomic charge divided by the total negative charge |
| `MeanElectronegativity` | Pauling | `mean(EN_i)` | Mean Pauling electronegativity of the atoms present |
| `ElectronegativityVariance` | Pauling² | `var(EN_i)` | Variance of the atomic electronegativities, a composition-only polarity indicator |
| `MaxBondPolarity` | Pauling | `max\|EN_i - EN_j\|` | Largest electronegativity difference across any bond; flags the most polarised (and often most reactive) bond |
| `PolarAtomCharge` | e | `sum\|q_i\| over N,O,S,P` | Total absolute charge carried by N, O, S and P atoms, the hydrogen bonding centres |
| `PPSA1` | Å² | `sum(SA_i), q_i > 0` | Partial positive surface area: solvent accessible surface of all positively charged atoms (Stanton-Jurs CPSA family) |
| `PNSA1` | Å² | `sum(SA_i), q_i < 0` | Partial negative surface area: solvent accessible surface of all negatively charged atoms |
| `PPSA2` | e·Å² | `sum(q+) * sum(SA+)` | Total charge weighted positive surface area: the positive surface multiplied by the total positive charge |
| `PNSA2` | e·Å² | `sum(q-) * sum(SA-)` | Total charge weighted negative surface area |
| `PPSA3` | e·Å² | `sum(q_i * SA_i), q_i > 0` | Atomic charge weighted positive surface area; each positive atom contributes its own charge times its own surface |
| `PNSA3` | e·Å² | `sum(q_i * SA_i), q_i < 0` | Atomic charge weighted negative surface area, one of the strongest CPSA descriptors for hydrogen bond acceptor strength |
| `DPSA1` | Å² | `PPSA1 - PNSA1` | Difference between the positive and negative partial surface areas |
| `FPSA1` | fraction | `PPSA1 / SASA` | Fractional positive surface area, PPSA1 over the total surface; size independent so it compares across molecules |
| `FNSA1` | fraction | `PNSA1 / SASA` | Fractional negative surface area, PNSA1 over the total surface |
| `WPSA1` | Å⁴ | `PPSA1 * SASA / 1000` | Surface weighted positive surface area, PPSA1 times the total surface divided by 1000 |
| `WNSA1` | Å⁴ | `PNSA1 * SASA / 1000` | Surface weighted negative surface area |
| `RPCS` | Å² | `SA(q_max) * q_max / sum(q+)` | Relative positive charged surface area: the surface of the most positive atom scaled by its share of the positive charge; a hydrogen bond donor strength descriptor |
| `RNCS` | Å² | `SA(q_min) * \|q_min\| / sum\|q-\|` | Relative negative charged surface area, the acceptor counterpart of RPCS |

---

## Quantum

Frontier orbital properties from a simple Huckel molecular orbital calculation on the pi system, with the ionisation-fitted parametrisation alpha = -6.15 eV and beta = -3.32 eV and Streitwieser heteroatom parameters. The delocalisation energies reproduce the textbook Huckel values exactly (benzene 2.000 beta, butadiene 0.472, naphthalene 3.683) and the ionisation potentials land within a few tenths of an eV of experiment for aromatics. Saturated molecules have no pi system; for those the frontier energies fall back to the ionisation energy of the highest lone pair or sigma bond present, so they stay molecule dependent.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `HOMOEnergy` | eV | `alpha + x_HOMO * beta` | Energy of the highest occupied molecular orbital from a Huckel calculation on the pi system (alpha = -6.15 eV, beta = -3.32 eV). Molecules without a pi system fall back to the ionisation energy of their highest lone pair or sigma bond, so the value stays molecule dependent |
| `LUMOEnergy` | eV | `alpha + x_LUMO * beta` | Energy of the lowest unoccupied molecular orbital; the more negative, the easier the molecule accepts an electron |
| `HOMOLUMOGap` | eV | `E_LUMO - E_HOMO` | Frontier orbital gap. Large for saturated and isolated systems, small for extended conjugation; it governs UV absorption and kinetic stability |
| `ChemicalPotential` | eV | `(E_HOMO + E_LUMO)/2` | Electronic chemical potential mu = (E_HOMO + E_LUMO)/2, the negative of the Mulliken electronegativity; measures the escaping tendency of electrons |
| `ChemicalHardness` | eV | `(E_LUMO - E_HOMO)/2` | Chemical hardness eta = (E_LUMO - E_HOMO)/2; hard molecules resist charge transfer, soft ones polarise easily |
| `Softness` | 1/eV | `1 / eta` | Global softness, the reciprocal of the hardness |
| `IonizationPotential` | eV | `-E_HOMO` | Koopmans ionisation potential, IP = -E(HOMO). For benzene this gives 9.47 eV against an experimental 9.24 eV |
| `ElectronAffinity` | eV | `-E_LUMO` | Koopmans electron affinity, EA = -E(LUMO) |
| `MullikenElectronegativity` | eV | `(IP + EA)/2` | Absolute electronegativity chi = (IP + EA)/2; the driving force for charge transfer between two molecules |
| `ElectrophilicityIndex` | eV | `mu² / (2*eta)` | Parr electrophilicity omega = mu²/(2*eta), the energy stabilisation on saturating the molecule with electrons; a standard reactivity ranking for Michael acceptors and electrophilic toxicants |
| `NucleophilicityIndex` | 1/eV | `1 / omega` | Reciprocal of the electrophilicity index, ranking electron donors |
| `ElectroacceptingPower` | eV | `(IP+3EA)² / (16(IP-EA))` | Gazquez electron accepting power omega+ = (IP + 3*EA)²/(16*(IP - EA)), the propensity to accept charge |
| `ElectrodonatingPower` | eV | `(3IP+EA)² / (16(IP-EA))` | Gazquez electron donating power omega- = (3*IP + EA)²/(16*(IP - EA)) |
| `PiElectronCount` | count | `count(pi electrons)` | Number of electrons in the pi system, counting lone pairs donated by heteroatoms in aromatic rings |
| `PiSystemSize` | count | `count(pi centres)` | Number of atoms contributing a p orbital to the pi system |
| `PiSystemCount` | count | `count(conjugated components)` | Number of separate conjugated fragments in the molecule |
| `TotalPiEnergy` | beta | `sum over occupied(2*x_i)` | Total Huckel pi electron energy in units of beta (8.000 for benzene); larger values mean more pi bonding |
| `DelocalizationEnergy` | beta | `E_pi - 2*n_pairs` | Resonance energy in units of beta: the extra pi stabilisation over the same number of isolated double bonds. Benzene 2.000, butadiene 0.472, naphthalene 3.683, matching the textbook Huckel values |
| `GapPerPiAtom` | eV | `gap / n_pi` | HOMO-LUMO gap divided by the size of the pi system; tracks how far conjugation extends rather than how large the molecule is |
| `HardnessRatio` | ratio | `eta / chi` | Chemical hardness over electronegativity, a scale free reactivity ratio |

---

## Fingerprints

Bit-vector structural keys. Reported here as summary values; use the fingerprint engine directly when the full bit vectors are needed.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `MorganFingerprint` | bits | `circular_fingerprint` | Morgan (extended connectivity) fingerprint bit count, describing circular atom environments |
| `MACCSFingerprint` | bits | `structural_keys` | MACCS structural key fingerprint: 166 predefined substructure keys |
| `TopologicalFingerprint` | bits | `topological_features` | Path-based (Daylight-style) topological fingerprint bit count |
| `AtomPairFingerprint` | bits | `atom_pairs` | Atom pair fingerprint: hashed pairs of atom types with their topological separation |

---

## Hybrid

Physicochemical property models and the medicinal chemistry filters built on them: polar surface area, solubility, volume, and the Lipinski/Ghose/Veber/Egan/Muegge rule sets. Violation counts are always 'number of criteria broken', so **zero is the good outcome**.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `Lipophilicity (LogP)` | logP | `fragment_based_logP` | Octanol-water partition coefficient from the molecular lipophilicity potential (MLP) integrated over the solvent accessible surface; the primary lipophilicity descriptor |
| `Polarizability` | Å³ | `sum(atomic_polarizabilities)` | Molecular polarizability as the sum of atomic polarizability increments; governs dispersion interactions and refractive index |
| `MolarRefractivity` | cm³/mol | `fragment_based_calculation` | Molar refractivity from atomic increments, a combined size and polarizability term used by the Ghose filter |
| `LipinskiHBA` | count | `lipinski_hba` | Lipinski hydrogen bond acceptor count (nitrogen plus oxygen atoms), the counting rule used by the rule of five |
| `LipinskiHBD` | count | `lipinski_hbd` | Lipinski hydrogen bond donor count (N-H and O-H groups) |
| `LipinskiViolationCount` | count | `lipinski_violations` | Number of Lipinski rule of five criteria violated: MW > 500, logP > 5, HBA > 10, HBD > 5. Two or more violations suggest poor oral absorption |
| `VeberTPSA` | Å² | `veber_psa` | Polar surface area estimated from per-element increments, the form used in the original Veber analysis (see TPSA for the Ertl fragment version) |
| `DrugLikenessScore` | score | `drug_score` | Composite drug-likeness score on a 0-1 scale built from Lipinski violations, molecular weight, logP and rotatable bond count |
| `SyntheticAccessibility` | score | `SA_score` | Synthetic accessibility estimate on a 1-10 scale (lower is easier), penalising ring systems, flexibility and heteroatom richness |
| `Fsp3` | fraction | `sp3/total_C` | Fraction of sp3 hybridised carbons; a measure of three-dimensional character that correlates with clinical success |
| `TPSA` | Å² | `sum(Ertl fragment contributions)` | Topological polar surface area from Ertl's fragment contributions for nitrogen and oxygen. Reproduces the published values exactly (aspirin 63.60, caffeine 58.44 Å²). Below ~90 Å² is associated with brain penetration and below ~140 Å² with good oral absorption |
| `TPSAWithSP` | Å² | `sum(contributions incl. S,P)` | Topological polar surface area counting sulfur and phosphorus as polar as well; the variant preferred for organophosphates and sulfur rich drugs |
| `TPSAPerHeavyAtom` | Å² | `TPSA / n_heavy` | TPSA divided by the heavy atom count, comparable across molecular sizes |
| `ESOLLogS` | log(mol/L) | `0.16 - 0.63logP - 0.0062MW + 0.066RB - 0.74AP` | Delaney ESOL estimate of aqueous solubility as log10 of the molar solubility: LogS = 0.16 - 0.63*logP - 0.0062*MW + 0.066*RotB - 0.74*AromaticProportion. Values below -6 indicate practically insoluble compounds |
| `McGowanVolume` | cm³/mol/100 | `sum(V_atoms) - 6.56*n_bonds` | McGowan characteristic volume, an additive volume built from atomic increments minus 6.56 per bond. It is the size term of Abraham's solvation equations |
| `HydrophilicFactor` | index | `Todeschini HF` | Todeschini hydrophilic factor, built from the number of hydrophilic groups (OH, NH, SH) relative to carbon and heavy atom counts. Strongly negative for pure hydrocarbons and large for polyols and sugars |
| `MolecularFlexibility` | ratio | `n_rot / n_heavy` | Rotatable bonds per heavy atom; a size independent conformational freedom measure used in bioavailability models |
| `LogPPerHeavyAtom` | ratio | `logP / n_heavy` | LogP divided by heavy atom count, the lipophilicity efficiency scale used to compare fragments with leads |
| `PolarAtomFraction` | fraction | `n_polar / n_heavy` | N, O, S and P atoms as a fraction of the heavy atoms |
| `GhoseViolations` | count | `count(violations of 4 criteria)` | Number of Ghose filter criteria violated: MW 160-480, logP -0.4 to 5.6, 20-70 heavy atoms and molar refractivity 40-130. Zero means the compound sits inside the drug-like qualifying range |
| `VeberViolations` | count | `count(violations of 2 criteria)` | Veber oral bioavailability criteria violated: rotatable bonds <= 10 and TPSA <= 140 Å². Derived from rat oral bioavailability data |
| `EganViolations` | count | `count(violations of 2 criteria)` | Egan 'absorption egg' criteria violated: TPSA <= 131.6 Å² and logP <= 5.88, a two-property model of passive intestinal absorption |
| `MueggeViolations` | count | `count(violations of 9 criteria)` | Muegge pharmacophore filter criteria violated (MW, logP, TPSA, rings, carbon and heteroatom counts, rotatable bonds, HBA, HBD). A stricter drug-likeness screen than Lipinski |
| `RuleOfThreeViolations` | count | `count(violations of 5 criteria)` | Congreve rule of three criteria violated: MW <= 300, logP <= 3, HBD <= 3, HBA <= 3, rotatable bonds <= 3. The standard filter for fragment screening libraries |
| `LeadLikenessViolations` | count | `count(violations of 3 criteria)` | Lead-likeness criteria violated: MW 250-350, logP <= 3.5, rotatable bonds <= 7. Leads need room to grow during optimisation |
| `BBBScore` | score | `weighted property score` | Blood-brain barrier likelihood on a 0-1 scale, combining TPSA, logP, molecular weight, hydrogen bond donors and flexibility. Above ~0.7 is typical of CNS active drugs |
| `OralBioavailabilityScore` | fraction | `n_passed / 6` | Fraction of the six standard oral filters (Lipinski, Ghose, Veber, Egan, Muegge, lead-likeness) the molecule passes cleanly |

---

## Custom

Selection-aware descriptors, meaningful when descriptors are computed for a sub-selection of atoms rather than the whole molecule.

| Descriptor | Unit | Formula | Description |
| --- | --- | --- | --- |
| `SelectionSize` | count | `len(selection)` | Number of atoms in the current selection, useful when descriptors are computed for a fragment |
| `SelectionDensity` | atoms/Å³ | `selection_size / molecular_volume` | Atoms of the selection per unit molecular volume, a packing density measure |

## Appendix: the PyDes vector families

The calculator window runs the PyDes generators in addition to the descriptor engine, which
adds systematically generated families. They follow fixed naming patterns:

| Pattern | Meaning |
| --- | --- |
| `nAtoms`, `nBonds`, `nHeavyAtoms` | raw counts |
| `nC`, `nN`, `nO`, `nS`, `nP`, `nF`, `nCl`, `nBr`, `nI`, `nB`, `nSi`, `nH` | atom count per element |
| `Ratio_XC` | count of element X divided by the carbon count |
| `Prop_X` | count of element X divided by the total atom count |
| `nDegree1` ... `nDegree6` | number of atoms with that heavy atom degree |
| `nSingleBonds`, `nDoubleBonds`, `nTripleBonds`, `nAromaticBonds` and their `Prop_` forms | bond order distribution |
| `nRings`, `nRing3` ... `nRing12` | ring count by ring size |
| `ATS{k}m`, `ATS{k}Z`, `ATS{k}d` | Moreau-Broto autocorrelation at topological distance k (1-15), weighted by mass, atomic number and degree |
| `MATS{k}*`, `GATS{k}*` | the corresponding Moran and Geary autocorrelation variants |
| `ATS*_norm`, `ATS*_sq` | the same autocorrelations divided by the heavy atom count and squared |
| `DistPower1` ... `DistPower10` | sum of topological distances raised to that power |
| `WDist_m`, `WDist_Z`, `WDist_d` | distance sums weighted by mass, atomic number and degree |

These are deliberately redundant with each other - they exist to give machine learning models a
large, uniformly generated feature space. The descriptor engine entries documented above are the
ones to reach for when a value has to be interpreted, cited, or compared with another program.

---

## Notes on use

* **Violation counts** (Lipinski, Ghose, Veber, Egan, Muegge, rule of three, lead-likeness)
  count broken criteria, so 0 means the molecule passes.
* **Hydrogen suppressed graph.** All topological indices ignore hydrogens, including explicit
  ones read from a file, which is what the published definitions assume.
* **Selections.** Every descriptor accepts an atom selection; counts, graph indices and
  surfaces are then computed for that sub-structure only.
* **Determinism.** Descriptors that require 3D coordinates depend on the generated geometry
  when the input has none. Supply a MOL/SDF/PDB file with coordinates for reproducible
  geometric and CPSA values.
