"""
Descriptor documentation export.

Writes a Markdown reference of every registered descriptor straight from the
engine, so the document can never drift away from the code. Run it after adding
or renaming a descriptor:

    python -m src.features.descriptor_calculator.docs_export
"""
import io
import os

from .descriptor_engine import DescriptorEngine
from .descriptor_types import DescriptorCategory

DEFAULT_OUTPUT = "MOLECULAR_DESCRIPTORS.md"

CATEGORY_NOTES = {
    "Constitutional": (
        "Counts and compositional ratios read straight off the molecular graph: elements, "
        "bonds, rings, functional groups and atom environments. They need no geometry and no "
        "charges, so they are always available, and they are the descriptors that are easiest "
        "to interpret when a model turns out to depend on them."),
    "Topological": (
        "Graph invariants computed on the **hydrogen suppressed** molecular graph, as the "
        "original definitions require. They encode size, branching, shape and heteroatom "
        "placement without any 3D information, which makes them reproducible and conformation "
        "independent - the reason they dominate classical QSAR."),
    "Geometric": (
        "Descriptors that need a 3D geometry: surface, volume, moments of inertia, shape "
        "ratios and extent along the principal axes. When a structure has no coordinates "
        "(for example one parsed from SMILES) PyChem generates them once with its own 3D "
        "generator, so these values are never silently zero."),
    "Electronic": (
        "Partial charge statistics and the charged partial surface area (CPSA) family of "
        "Stanton and Jurs, which combine Gasteiger charges with the per-atom solvent "
        "accessible surface. Gasteiger charges are computed on demand when a molecule does "
        "not carry them."),
    "Quantum": (
        "Frontier orbital properties from a simple Huckel molecular orbital calculation on "
        "the pi system, with the ionisation-fitted parametrisation alpha = -6.15 eV and "
        "beta = -3.32 eV and Streitwieser heteroatom parameters. The delocalisation energies "
        "reproduce the textbook Huckel values exactly (benzene 2.000 beta, butadiene 0.472, "
        "naphthalene 3.683) and the ionisation potentials land within a few tenths of an eV "
        "of experiment for aromatics. Saturated molecules have no pi system; for those the "
        "frontier energies fall back to the ionisation energy of the highest lone pair or "
        "sigma bond present, so they stay molecule dependent."),
    "Fingerprints": (
        "Bit-vector structural keys. Reported here as summary values; use the fingerprint "
        "engine directly when the full bit vectors are needed."),
    "Hybrid": (
        "Physicochemical property models and the medicinal chemistry filters built on them: "
        "polar surface area, solubility, volume, and the Lipinski/Ghose/Veber/Egan/Muegge "
        "rule sets. Violation counts are always 'number of criteria broken', so **zero is "
        "the good outcome**."),
    "Custom": (
        "Selection-aware descriptors, meaningful when descriptors are computed for a "
        "sub-selection of atoms rather than the whole molecule."),
}

PYDES_SECTION = """
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
"""


def build_markdown(engine=None):
    """Return the full descriptor reference as a Markdown string."""
    engine = engine or DescriptorEngine()
    total = sum(len(d) for d in engine.descriptors.values())
    filled = [c for c in engine.descriptors if engine.descriptors[c]]

    lines = ["# PyChem-Pro Molecular Descriptors", ""]
    lines.append("Reference for every descriptor produced by the **Molecular Descriptor "
                 "Calculator** (`Tools -> Molecular Descriptor Calculator`).")
    lines.append("")
    lines.append("The descriptor engine registers **%d descriptors across %d categories**. "
                 "The calculator additionally emits the systematically generated PyDes "
                 "vectors described in the appendix, which brings a typical CSV export to "
                 "roughly 700 columns per molecule." % (total, len(filled)))
    lines.append("")
    lines.append("Everything is computed with numpy and the project's own chemistry code; "
                 "there are no external cheminformatics dependencies.")
    lines.append("")
    lines.append("## Contents")
    lines.append("")
    for category in DescriptorCategory:
        descriptors = engine.descriptors.get(category, [])
        if not descriptors:
            continue
        anchor = category.value.lower().replace(' ', '-')
        lines.append("- [%s](#%s) - %d descriptors"
                     % (category.value, anchor, len(descriptors)))
    lines.append("- [Appendix: the PyDes vector families](#appendix-the-pydes-vector-families)")
    lines.append("")

    for category in DescriptorCategory:
        descriptors = engine.descriptors.get(category, [])
        if not descriptors:
            continue
        lines.append("---")
        lines.append("")
        lines.append("## %s" % category.value)
        lines.append("")
        note = CATEGORY_NOTES.get(category.value)
        if note:
            lines.append(note)
            lines.append("")
        lines.append("| Descriptor | Unit | Formula | Description |")
        lines.append("| --- | --- | --- | --- |")
        for desc in descriptors:
            formula = (desc.formula or "").replace("|", "\\|")
            description = (desc.description or "").replace("|", "\\|")
            lines.append("| `%s` | %s | `%s` | %s |"
                         % (desc.name, desc.unit or "-", formula or "-", description))
        lines.append("")

    lines.append(PYDES_SECTION.strip())
    lines.append("")
    return "\n".join(lines)


def write_markdown(path=DEFAULT_OUTPUT, engine=None):
    """Write the reference document, returning the path it was written to."""
    text = build_markdown(engine)
    with io.open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return os.path.abspath(path)


if __name__ == "__main__":
    import sys
    target = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUTPUT
    written = write_markdown(target)
    print("Descriptor documentation written to %s" % written)
