# MMFF94 Reference Data

This directory ships verbatim copies of Jmol's MMFF94 parameter files.

## Files

- **`mmff94.par.txt`** — Merck Molecular Force Field 1994 parameter tables
  (bond, angle, stretch-bend, torsion, out-of-plane, van der Waals,
  partial bond charge increments, formal-charge adjustments). 5077 lines,
  9 concatenated sections marked by `N. MMFFxxx.PAR:` headers.

- **`mmff94_atom_types.txt`** — MMFF94 atom-type definitions. 260 lines.
  Phase 1 uses only the `mmType` (column 3) and `class` (column 4) fields
  to build the type→class lookup table. The SMARTS column is unused until
  Phase 3.

- **`_cache/`** — gitignored. Holds pickled parsed parameter dicts for
  fast subsequent startup (~5 ms warm vs ~150 ms cold parse).

## Source

`/Users/gauravmasand/Developer/jmol-16.3.33/src/org/jmol/minimize/forcefield/data/`
(Jmol 16.3.33, LGPL 2.1).

## Updating

Re-copy from the Jmol source; do not edit in place. The parser
(`src/services/forcefield/parameters.py`) tolerates Jmol's exact format
and will fail loudly if columns shift.
