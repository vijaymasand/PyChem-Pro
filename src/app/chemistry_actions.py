"""
Chemistry Actions — Optimize, charges, descriptors, SMILES generation UI actions.

All functions receive the MainWindow instance as ``window``.
"""

from src.shared.qt_compat import QApplication, QMessageBox, QColorDialog, Qt
from src.core.domain.models.bond import BondType

_DEBUG = False


# ── Geometry Optimization ─────────────────────────────────────────

def optimize_geometry(window, checked=False, method=None):
    if not window.molecule:
        return
    if method is None or not isinstance(method, str):
        method = "MMFF94"

    if "AM1" in method:
        window.status_bar.showMessage("Optimizing geometry (AM1)...")
    else:
        window.status_bar.showMessage("Optimizing geometry (MMFF94)...")
    QApplication.processEvents()
    try:
        if "AM1" in method:
            from src.features.cheminformatics.services.am1 import am1_optimize_geometry
            success = am1_optimize_geometry(window.molecule, max_steps=50)
            window.viewer_3d.set_molecule(window.molecule)
            status = "converged" if success else "max steps reached"
            window.status_bar.showMessage(f"Optimized (AM1: {status})")
        else:
            # Use the new MMFF94Service which adds explicit hydrogens,
            # assigns BCI partial charges, and runs a proper force field
            # (bond stretch + angle bending + torsion + VdW).
            from pychem._bridge import get_registry
            service = get_registry().forcefield
            pre_atoms = len(window.molecule.atoms)
            result = service.optimize_geometry(
                window.molecule, max_iters=500, method='steepest_descent'
            )
            added_h = len(window.molecule.atoms) - pre_atoms
            window.input_panel.update_molecule_info(window.molecule)
            window.viewer_3d.set_molecule(window.molecule)
            window.viewer_2d.set_molecule(window.molecule)
            status = "converged" if result.converged else f"stopped ({result.num_steps} steps)"
            window.status_bar.showMessage(
                f"Optimized (MMFF94: {status}, {added_h} H added, E={result.final_energy:.2f})"
            )
    except Exception as e:
        window.status_bar.showMessage(f"Optimization error: {e}")


# ── Charge Computation ───────────────────────────────────────────

def compute_charges(window, checked=False, method=None):
    if not window.molecule:
        return
    if method is None or not isinstance(method, str):
        method = "Gasteiger"

    if "PM3" in method:
        window.status_bar.showMessage("Computing PM3 charges...")
    elif "AM1" in method:
        window.status_bar.showMessage("Computing AM1 charges...")
    elif "MMFF94" in method:
        window.status_bar.showMessage("Computing MMFF94 charges...")
    else:
        window.status_bar.showMessage("Computing Gasteiger charges...")
    QApplication.processEvents()
    try:
        if "PM3" in method:
            from src.features.cheminformatics.services.pm3 import pm3_assign_charges
            success = pm3_assign_charges(window.molecule)
            if success:
                window.input_panel.update_molecule_info(window.molecule)
                window.viewer_3d.update()
                window.status_bar.showMessage("PM3 partial charges assigned successfully")
            else:
                window.status_bar.showMessage("PM3 charge calculation failed")
        elif "AM1" in method:
            from src.features.cheminformatics.services.am1 import am1_assign_charges
            success = am1_assign_charges(window.molecule)
            if success:
                window.input_panel.update_molecule_info(window.molecule)
                window.viewer_3d.update()
                window.status_bar.showMessage("AM1 partial charges assigned successfully")
            else:
                window.status_bar.showMessage("AM1 charge calculation failed")
        elif "MMFF94" in method:
            # Use the new MMFF94Service BCI charge assignment
            from pychem._bridge import get_registry
            try:
                get_registry().forcefield.assign_charges(window.molecule)
                window.input_panel.update_molecule_info(window.molecule)
                window.viewer_3d.update()
                window.status_bar.showMessage("MMFF94 partial charges assigned successfully")
            except Exception as e:
                window.status_bar.showMessage(f"MMFF94 charge calculation failed: {e}")
        else:
            from src.features.cheminformatics.electrostatics.gasteiger import compute_gasteiger_charges
            compute_gasteiger_charges(window.molecule)
            window.input_panel.update_molecule_info(window.molecule)
            window.viewer_3d.update()
            window.status_bar.showMessage("Gasteiger charges computed")
    except Exception as e:
        window.status_bar.showMessage(f"Charge computation error: {e}")


# ── Descriptor Calculator ────────────────────────────────────────

def open_descriptor_calculator(window):
    """Open molecular descriptor calculator dialog."""
    if _DEBUG:
        print(f"[DEBUG] Opening descriptor calculator...")
        print(f"[DEBUG] Molecule available: {window.molecule is not None}")

    if not window.molecule:
        QMessageBox.warning(window, "No Molecule",
                            "Please load a molecule first before calculating descriptors.")
        return

    try:
        if _DEBUG:
            print(f"[DEBUG] Molecule atoms: {len(window.molecule.atoms)}, bonds: {len(window.molecule.bonds)}")

            if hasattr(window.molecule.atoms[0], 'x'):
                print(f"[DEBUG] Molecule has coordinates")

        from src.features.descriptor_calculator.gui import show_descriptor_calculator
        dialog = show_descriptor_calculator(window.molecule, window)

        if _DEBUG:
            print(f"[DEBUG] Dialog created: {type(dialog)}")

        if hasattr(dialog, '__name__') and dialog.__name__ == 'show_descriptor_calculator':
            if _DEBUG:
                print(f"[DEBUG] ERROR: Dummy function was called")
            QMessageBox.critical(window, "GUI Error",
                              "The descriptor calculator GUI is not available.\n"
                              "Qt framework is not working properly.\n"
                              "Please use the alternative Tkinter GUI:\n"
                              "python tkinter_descriptor_gui.py")
            return

        window.descriptor_calculator_dialog = dialog

        print(f"[DEBUG] Dialog stored as instance variable")
        print(f"[DEBUG] Descriptor calculator should now be visible")

    except ImportError as e:
        print(f"[DEBUG] Import Error: {e}")
        QMessageBox.critical(window, "Import Error",
                          f"Could not import descriptor calculator:\n{str(e)}")
    except Exception as e:
        print(f"[DEBUG] General Error: {e}")
        import traceback
        traceback.print_exc()
        QMessageBox.critical(window, "Error",
                          f"Error opening descriptor calculator:\n{str(e)}")


# ── Aromaticity ──────────────────────────────────────────────────

def perceive_aromaticity_action(window):
    """Manually trigger aromaticity perception on the loaded molecule."""
    if not window.molecule:
        return
    window.status_bar.showMessage("Perceiving aromaticity...")
    QApplication.processEvents()
    try:
        from src.features.smiles_parser.rules.aromaticity import perceive_aromaticity
        perceive_aromaticity(window.molecule)
        window.molecule.assign_sybyl_types()
        window.viewer_3d.update()
        window.viewer_2d.update()
        window.status_bar.showMessage("Aromaticity perceived successfully and SYBYL types updated")
    except Exception as e:
        window.status_bar.showMessage(f"Aromaticity perception failed: {e}")
        QMessageBox.warning(window, "Aromaticity Error", f"Perception failed: {e}")


# ── SMILES Generation ───────────────────────────────────────────

def generate_smiles(window):
    """Generate SMILES notation from loaded 3D structure using BKChem algorithms."""
    if not window.molecule:
        QMessageBox.warning(window, "No Molecule", "Please load a molecule first.")
        return

    try:
        try:
            from src.features.smiles_generator import generate_smiles as _gen_smiles
            smiles = _gen_smiles(window.molecule)
            if smiles and len(smiles) > 0:
                print(f"[DEBUG] BKChem SMILES generated: {smiles}")
            else:
                raise ValueError("Empty SMILES from BKChem generator")
        except Exception as bkchem_error:
            print(f"[DEBUG] BKChem SMILES failed: {bkchem_error}, falling back to basic")
            smiles = _molecule_to_smiles_basic(window, window.molecule)

        try:
            from src.features.smiles_parser.services.parser import parse_smiles
            test_mol = parse_smiles(smiles)
            if not test_mol or len(test_mol.atoms) == 0:
                raise ValueError("Generated SMILES could not be parsed")
        except Exception:
            smiles = _generate_simple_formula(window)
            window.status_bar.showMessage("Generated molecular formula instead of SMILES due to complexity")

        if hasattr(window, 'console') and window.console:
            try:
                window.console._append_output(f"\n# Generated SMILES/Formula")
                window.console._append_output(f"smiles = '{smiles}'")
                window.console._append_output(f"# Generated: {smiles}")
            except Exception as console_err:
                print(f"[DEBUG] Could not write to console: {console_err}")

        msg = QMessageBox(window)
        msg.setWindowTitle("SMILES Generated")
        msg.setText(f"SMILES notation has been generated and printed to the console.")
        msg.setDetailedText(f"SMILES: {smiles}")
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()

        window.status_bar.showMessage(f"SMILES generated: {smiles}")

    except Exception as e:
        QMessageBox.critical(window, "SMILES Generation Error",
                           f"Failed to generate SMILES: {str(e)}")
        window.status_bar.showMessage("SMILES generation failed")


def copy_smiles_to_clipboard(window):
    """Generate SMILES for the current molecule and copy to clipboard."""
    if not window.molecule:
        window.status_bar.showMessage("No molecule loaded")
        return
    try:
        from src.vendors.oasa_bridge import domain_to_oasa_mol
        import src.vendors.oasa.smiles as oasa_smiles

        o_mol, _ = domain_to_oasa_mol(window.molecule)
        sm_engine = oasa_smiles.smiles()
        smiles_str = sm_engine.get_smiles(o_mol)

        if smiles_str:
            clipboard = QApplication.clipboard()
            clipboard.setText(smiles_str)
            window.status_bar.showMessage(f"SMILES copied: {smiles_str}")
        else:
            window.status_bar.showMessage("Could not generate SMILES for this molecule")
    except Exception as e:
        window.status_bar.showMessage(f"SMILES generation failed: {e}")


# ── SMILES helpers (internal) ────────────────────────────────────

def _molecule_to_smiles_basic(window, molecule):
    """Basic SMILES generation - attempts to generate SMILES, falls back to formula."""
    if not molecule or not molecule.atoms:
        return ""

    print(f"DEBUG: Attempting SMILES generation for {len(molecule.atoms)} atoms")

    try:
        visited = set()
        fragments = []

        for start_idx in range(len(molecule.atoms)):
            if start_idx not in visited:
                fragment_atoms = _dfs_fragment(molecule, start_idx, visited)
                print(f"DEBUG: Found fragment with {len(fragment_atoms)} atoms")
                if fragment_atoms:
                    fragment_smiles = _fragment_to_smiles(window, molecule, fragment_atoms)
                    print(f"DEBUG: Fragment SMILES: {fragment_smiles}")
                    fragments.append(fragment_smiles)

        smiles_result = '.'.join(fragments)
        print(f"DEBUG: Combined SMILES result: {smiles_result}")

        if not smiles_result:
            print("DEBUG: Empty SMILES result, falling back to formula")
            raise ValueError("Empty SMILES result")

        if any(c in smiles_result for c in ['(', ')', '[', ']']):
            if smiles_result.count('(') != smiles_result.count(')') or smiles_result.count('[') != smiles_result.count(']'):
                print("DEBUG: Unmatched brackets, falling back to formula")
                raise ValueError("Unmatched brackets")

        print(f"DEBUG: Successfully generated SMILES: {smiles_result}")
        return smiles_result

    except Exception as e:
        print(f"DEBUG: Basic SMILES generation failed: {e}. Falling back to formula.")
        return _generate_simple_formula(window)


def _dfs_fragment(molecule, start_idx, visited):
    """DFS to find all atoms in a connected fragment."""
    fragment = []
    stack = [start_idx]

    while stack:
        idx = stack.pop()
        if idx in visited:
            continue

        visited.add(idx)
        fragment.append(idx)

        for neighbor_idx, _ in molecule._adjacency.get(idx, []):
            if neighbor_idx not in visited:
                stack.append(neighbor_idx)

    return fragment


def _fragment_to_smiles(window, molecule, fragment_atoms):
    """Convert a connected fragment to SMILES."""
    if not fragment_atoms:
        return ""

    print(f"DEBUG: Converting fragment with {len(fragment_atoms)} atoms")

    if len(fragment_atoms) <= 2:
        return _linear_smiles(molecule, fragment_atoms)

    if len(fragment_atoms) <= 20:
        return _simple_chain_smiles(molecule, fragment_atoms)

    if len(fragment_atoms) <= 200:
        return _path_based_smiles(window, molecule, fragment_atoms)

    print("DEBUG: Fragment too complex, falling back to formula")
    raise ValueError("Fragment too complex for basic SMILES generation")


def _path_based_smiles(window, molecule, fragment_atoms):
    """Generate SMILES using a simple neighbor-based approach."""
    if len(fragment_atoms) == 1:
        atom = molecule.atoms[fragment_atoms[0]]
        return atom.element.symbol

    connectivity = {}
    for idx in fragment_atoms:
        connectivity[idx] = []
        for neighbor_idx, bond_idx in molecule._adjacency.get(idx, []):
            if neighbor_idx in fragment_atoms:
                bond = molecule.bonds[bond_idx]
                if molecule.atoms[neighbor_idx].element.symbol != 'H':
                    connectivity[idx].append((neighbor_idx, bond))

    start_idx = max(fragment_atoms,
                   key=lambda i: len(connectivity.get(i, [])))

    visited = set()
    result = ""

    def build_smiles(current_idx, parent_idx=None):
        nonlocal result
        if current_idx in visited:
            return ""

        visited.add(current_idx)
        atom = molecule.atoms[current_idx]

        symbol = atom.element.symbol
        if symbol in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']:
            result += symbol
        else:
            result += f'[{symbol}]'

        neighbors = []
        for neighbor_idx, bond in connectivity.get(current_idx, []):
            if neighbor_idx != parent_idx and neighbor_idx not in visited:
                neighbors.append((neighbor_idx, bond))

        if not neighbors:
            return ""

        neighbors.sort(key=lambda x: molecule.atoms[x[0]].element.atomic_number)

        main_neighbor, main_bond = neighbors[0]
        bond_symbol = _bond_to_symbol(main_bond)
        result += bond_symbol
        build_smiles(main_neighbor, current_idx)

        for neighbor_idx, bond in neighbors[1:3]:
            bond_symbol_br = _bond_to_symbol(bond)
            temp_visited = visited.copy()

            branch_smiles = _build_branch(molecule, connectivity, neighbor_idx, current_idx, temp_visited)
            print(f"DEBUG: Branch {neighbor_idx}: {branch_smiles}")

            if branch_smiles and branch_smiles != "" and branch_smiles != "()":
                result += '(' + bond_symbol_br + branch_smiles + ')'
                print(f"DEBUG: Added branch: ({bond_symbol_br}{branch_smiles})")
            else:
                print(f"DEBUG: Skipped empty branch")

        return result

    def _build_branch(molecule, connectivity, start_idx_br, parent_idx, visited_set):
        """Build a branch SMILES string."""
        if start_idx_br in visited_set:
            return ""

        visited_set.add(start_idx_br)
        atom = molecule.atoms[start_idx_br]

        branch = atom.element.symbol
        if atom.element.symbol not in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']:
            branch = f'[{atom.element.symbol}]'

        neighbors = []
        for neighbor_idx, bond in connectivity.get(start_idx_br, []):
            if neighbor_idx != parent_idx and neighbor_idx not in visited_set:
                neighbors.append((neighbor_idx, bond))

        if neighbors:
            neighbor_idx, bond = neighbors[0]
            bond_symbol_br = _bond_to_symbol(bond)
            branch += bond_symbol_br
            branch += _build_branch(molecule, connectivity, neighbor_idx, start_idx_br, visited_set)

        return branch

    current_result = ""
    while len(visited) < len(fragment_atoms):
        if not current_result:
            build_smiles(start_idx)
            current_result = result
        else:
            unvisited = [idx for idx in fragment_atoms if idx not in visited]
            if unvisited:
                next_start = None
                for visited_idx in visited:
                    for neighbor_idx, bond in connectivity.get(visited_idx, []):
                        if neighbor_idx in unvisited:
                            next_start = neighbor_idx
                            bond_symbol_c = _bond_to_symbol(bond)
                            current_result += bond_symbol_c
                            break
                    if next_start:
                        break

                if next_start:
                    build_smiles(next_start)
                    current_result = result
                else:
                    current_result += "."
                    build_smiles(unvisited[0])
                    current_result = result
            else:
                break

    result = current_result

    print(f"DEBUG: Final SMILES result: {result}")
    print(f"DEBUG: Visited {len(visited)} out of {len(fragment_atoms)} atoms")

    return result


def _simple_chain_smiles(molecule, fragment_atoms):
    """Generate SMILES for simple chain molecules."""
    if len(fragment_atoms) == 1:
        atom = molecule.atoms[fragment_atoms[0]]
        return atom.element.symbol

    result = ""
    used_atoms = {fragment_atoms[0]}
    current_atom = fragment_atoms[0]

    while len(used_atoms) < len(fragment_atoms):
        atom = molecule.atoms[current_atom]
        result += atom.element.symbol

        next_atom = None
        for neighbor_idx, bond_idx in molecule._adjacency.get(current_atom, []):
            if neighbor_idx not in used_atoms:
                next_atom = neighbor_idx
                used_atoms.add(neighbor_idx)

                bond = molecule.bonds[bond_idx]
                result += _bond_to_symbol(bond)
                break

        if next_atom is None:
            break
        current_atom = next_atom

    if current_atom is not None:
        result += molecule.atoms[current_atom].element.symbol

    print(f"DEBUG: Simple chain result: {result}")
    return result


def _linear_smiles(molecule, fragment_atoms):
    """Generate SMILES for simple linear molecules."""
    if len(fragment_atoms) == 1:
        atom = molecule.atoms[fragment_atoms[0]]
        return atom.element.symbol

    idx1, idx2 = fragment_atoms
    atom1 = molecule.atoms[idx1]
    atom2 = molecule.atoms[idx2]

    bond_symbol = "-"
    for neighbor_idx, bond_idx in molecule._adjacency.get(idx1, []):
        if neighbor_idx == idx2:
            bond = molecule.bonds[bond_idx]
            bond_symbol = _bond_to_symbol(bond)
            break

    return f"{atom1.element.symbol}{bond_symbol}{atom2.element.symbol}"


def _atom_to_smiles_recursive(window, molecule, atom_idx, parent_idx,
                              visited, fragment_atoms):
    """Recursively convert atom and its neighbors to SMILES."""
    if atom_idx in visited:
        return ""

    visited.add(atom_idx)
    atom = molecule.atoms[atom_idx]

    symbol = atom.element.symbol
    if symbol in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']:
        atom_smiles = symbol
    else:
        atom_smiles = f'[{symbol}]'

    neighbors = []
    for neighbor_idx, bond_idx in molecule._adjacency.get(atom_idx, []):
        if neighbor_idx != parent_idx and neighbor_idx in fragment_atoms:
            bond = molecule.bonds[bond_idx]
            neighbors.append((neighbor_idx, bond))

    neighbors.sort(key=lambda x: molecule.atoms[x[0]].element.atomic_number)

    branches = []
    ring_bonds = []
    next_neighbor = None

    for neighbor_idx, bond in neighbors:
        if neighbor_idx in visited:
            if len(ring_bonds) < 2:
                ring_num = _get_ring_number(atom_idx, neighbor_idx, visited)
                ring_bonds.append((bond, ring_num))
        else:
            if not next_neighbor:
                next_neighbor = (neighbor_idx, bond)
            else:
                branches.append((neighbor_idx, bond))

    for bond, ring_num in ring_bonds:
        bond_symbol = _bond_to_symbol(bond)
        atom_smiles += bond_symbol + str(ring_num)

    if next_neighbor:
        bond_symbol = _bond_to_symbol(next_neighbor[1])
        chain_smiles = _atom_to_smiles_recursive(
            window, molecule, next_neighbor[0], atom_idx, visited, fragment_atoms)
        atom_smiles += bond_symbol + chain_smiles

    for i, (neighbor_idx, bond) in enumerate(branches[:2]):
        bond_symbol = _bond_to_symbol(bond)
        branch_smiles = _atom_to_smiles_recursive(
            window, molecule, neighbor_idx, atom_idx, visited, fragment_atoms)
        atom_smiles += '(' + bond_symbol + branch_smiles + ')'

    return atom_smiles


def _bond_to_symbol(bond):
    """Convert bond type to SMILES symbol."""
    if bond.bond_type == BondType.SINGLE:
        return ''
    elif bond.bond_type == BondType.DOUBLE:
        return '='
    elif bond.bond_type == BondType.TRIPLE:
        return '#'
    else:
        return ''


def _get_ring_number(atom1_idx, atom2_idx, visited):
    """Generate a ring number for ring closure."""
    return (len(visited) % 9) + 1


def _generate_simple_formula(window):
    """Generate molecular formula as fallback."""
    if not window.molecule:
        return ""

    atom_counts = {}
    for atom in window.molecule.atoms:
        symbol = atom.element.symbol
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

    formula_parts = []

    if 'C' in atom_counts:
        count = atom_counts['C']
        formula_parts.append(f"C{count if count > 1 else ''}")
        del atom_counts['C']

    if 'H' in atom_counts:
        count = atom_counts['H']
        formula_parts.append(f"H{count if count > 1 else ''}")
        del atom_counts['H']

    for symbol in sorted(atom_counts.keys()):
        count = atom_counts[symbol]
        formula_parts.append(f"{symbol}{count if count > 1 else ''}")

    return "".join(formula_parts)


# ── Lipophilicity Actions ─────────────────────────────────────

def calculate_logp_action(window):
    """Calculate LogP and print to console."""
    if not window.molecule:
        return
    window.status_bar.showMessage("Calculating LogP...")
    QApplication.processEvents()
    try:
        from src.features.cheminformatics.services.lipophilicity_service import LipophilicityService
        service = LipophilicityService(window.molecule)
        logp = service.calculate_logp()
        print(f"\nCalculated LogP for {getattr(window.molecule, 'name', 'Molecule') or 'Molecule'}: {logp:.4f}")
        window.status_bar.showMessage(f"LogP: {logp:.4f}")
    except Exception as e:
        window.status_bar.showMessage(f"LogP calculation error: {e}")


def show_lipophilic_contributions(window):
    """Label atoms with their lipophilic contributions."""
    if not window.molecule:
        return
    window.status_bar.showMessage("Calculating Lipophilic Contributions...")
    QApplication.processEvents()
    try:
        from src.features.cheminformatics.services.lipophilicity_service import LipophilicityService
        service = LipophilicityService(window.molecule)
        contribs = service.get_atomic_contributions()
        
        labels = {}
        selected = window.viewer_3d.selected_atoms
        for i, contrib in enumerate(contribs):
            if not selected or i in selected:
                labels[i] = f"{contrib:.2f}"
            
        window.viewer_3d.labels = labels
        window.viewer_3d.show_labels = True
        
        # Sync UI checkbox if it exists
        if hasattr(window, 'input_panel') and hasattr(window.input_panel, 'show_labels_check'):
            window.input_panel.show_labels_check.setChecked(True)
            
        window.viewer_3d.update()
        window.status_bar.showMessage("Lipophilic contributions labeled on atoms")
    except Exception as e:
        window.status_bar.showMessage(f"Lipophilicity error: {e}")


def color_by_lipophilicity(window):
    """Apply a gradient color mapping based on lipophilic contributions."""
    if not window.molecule:
        return
    
    # Prompt for two colors
    color_min = QColorDialog.getColor(Qt.GlobalColor.red, window, "Choose Color for Hydrophilic (min)")
    if not color_min.isValid(): return
    
    color_max = QColorDialog.getColor(Qt.GlobalColor.green, window, "Choose Color for Lipophilic (max)")
    if not color_max.isValid(): return
    
    window.status_bar.showMessage("Mapping lipophilicity gradient...")
    QApplication.processEvents()
    
    try:
        from src.features.cheminformatics.services.lipophilicity_service import LipophilicityService
        service = LipophilicityService(window.molecule)
        contribs = service.get_atomic_contributions()
        
        if not contribs: return
        
        c_min = min(contribs)
        c_max = max(contribs)
        c_range = c_max - c_min if c_max != c_min else 1.0
        
        custom_colors = {}
        for i, val in enumerate(contribs):
            # Linear interpolation between color_min and color_max
            t = (val - c_min) / c_range
            r = int(color_min.red() + t * (color_max.red() - color_min.red()))
            g = int(color_min.green() + t * (color_max.green() - color_min.green()))
            b = int(color_min.blue() + t * (color_max.blue() - color_min.blue()))
            custom_colors[i] = (r, g, b)
            
        window.viewer_3d.set_atom_colors(custom_colors)
        window.status_bar.showMessage(f"Gradient mapped (range: {c_min:.2f} to {c_max:.2f})")
    except Exception as e:
        window.status_bar.showMessage(f"Gradient mapping error: {e}")
