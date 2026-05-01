"""
ConversionWorker — Background SMILES-to-3D conversion thread.
"""

from src.shared.qt_compat import QObject, Signal


class ConversionWorker(QObject):
    """Worker for background SMILES conversion."""
    finished = Signal(object, str)  # (molecule, error_msg)
    progress = Signal(int)

    def __init__(self, smiles):
        super().__init__()
        self.smiles = smiles

    def run(self):
        try:
            self.progress.emit(10)

            # Parse SMILES — OASA's text_to_mol(calc_coords=1) already
            # generates high-quality native 2D coordinates stored as
            # x2d/y2d on each domain atom.  Using those directly avoids
            # the lossy domain→OASA→coords_generator round-trip that
            # was corrupting ring geometry and bond orders.
            from src.features.smiles_parser.services.parser import parse_smiles
            mol = parse_smiles(self.smiles, use_bkchem_tokenizer=True)
            self.progress.emit(30)

            # Use native OASA 2D coordinates (x2d/y2d) from parse_smiles.
            # Only fall back to CoordinateGenerator2DSMILES when x2d is missing.
            has_x2d = all(
                hasattr(a, 'x2d') and a.x2d is not None
                for a in mol.atoms if a.symbol != 'H'
            )

            if has_x2d:
                # Use the native OASA coordinates directly
                for atom in mol.atoms:
                    if hasattr(atom, 'x2d') and atom.x2d is not None:
                        atom.x = float(atom.x2d)
                        atom.y = float(atom.y2d)
                        atom.z = 0.0
                    else:
                        atom.x, atom.y, atom.z = 0.0, 0.0, 0.0
            else:
                # Fallback: re-generate coordinates
                from src.features.layout_2d.generators.coordgen2d_smiles_pure_oasa import CoordinateGenerator2DSMILES
                generator = CoordinateGenerator2DSMILES(mol)
                coords_2d = generator.generate()
                for atom in mol.atoms:
                    if atom.index in coords_2d:
                        x, y = coords_2d[atom.index]
                        atom.x = float(x)
                        atom.y = float(y)
                        atom.z = 0.0
                    else:
                        atom.x, atom.y, atom.z = 0.0, 0.0, 0.0
            self.progress.emit(70)

            # Compute charges
            from src.features.cheminformatics.electrostatics.gasteiger import compute_gasteiger_charges
            compute_gasteiger_charges(mol)
            self.progress.emit(90)

            # Assign types for export
            mol.assign_sybyl_types()
            self.progress.emit(100)

            self.finished.emit(mol, "")
        except Exception as e:
            self.finished.emit(None, str(e))
