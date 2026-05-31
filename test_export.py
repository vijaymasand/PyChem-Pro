import sys
from PySide6.QtWidgets import QApplication
from src.features.visualization_3d.ui.mol_viewer_3d import MolViewer3D
from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom

def test():
    app = QApplication(sys.argv)
    viewer = MolViewer3D()
    
    # Create a dummy PDB molecule
    mol = Molecule()
    mol.properties['is_protein'] = True
    atom1 = Atom("C")
    atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
    atom1.res_seq = 1
    atom1.chain_id = "A"
    atom1.res_name = "ALA"
    atom1.pdb_name = "CA"
    atom1.has_coords = True
    mol.atoms.append(atom1)
    
    atom2 = Atom("C")
    atom2.x, atom2.y, atom2.z = 1.0, 1.0, 1.0
    atom2.res_seq = 2
    atom2.chain_id = "A"
    atom2.res_name = "GLY"
    atom2.pdb_name = "CA"
    atom2.has_coords = True
    mol.atoms.append(atom2)
    
    viewer.set_molecule(mol)
    viewer.set_render_mode('cartoon')
    
    # Attempt high quality export
    success = viewer.export_image('test_output.png', 300, True)
    print("Export success:", success)

if __name__ == "__main__":
    test()
