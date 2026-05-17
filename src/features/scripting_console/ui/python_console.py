"""
Python Console — Interactive Python interpreter widget.

Provides a REPL-like console embedded in the application.
The current molecule is available as `mol` in the console namespace.
"""

import sys
import io
import traceback
from src.shared.qt_compat import QWidget, QVBoxLayout, QTextEdit, QLineEdit, QLabel, QHBoxLayout, Qt, Signal, QFont, QColor, QTextCursor
from src.shared.ui.theme import COLORS, theme_signals
from src.features.visualization_3d.services.docking_pose_service import DockingPoseService


class PythonConsole(QWidget):
    """
    Embedded Python interpreter with output display and command input.

    The molecule object is available as `mol` in the console namespace.
    Numpy is available as `np`.
    """

    command_executed = Signal(str)  # Emitted after each command with output

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(80)

        self._namespace = {}
        self._history = []
        self._history_idx = 0

        self._init_ui()
        self._init_namespace()

    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Header bar. Widgets are stored on self so _apply_theme()
        # can restyle them when the user switches themes at runtime.
        self._header = QWidget()
        self._header.setFixedHeight(26)
        header_layout = QHBoxLayout(self._header)
        header_layout.setContentsMargins(10, 0, 10, 0)
        self._header_label = QLabel("Python Console")
        header_layout.addWidget(self._header_label)
        header_layout.addStretch()
        self._header_help = QLabel("mol = current molecule  ·  np = numpy")
        header_layout.addWidget(self._header_help)
        layout.addWidget(self._header)

        mono_font = QFont('SF Mono', 10)
        mono_font.setStyleHint(QFont.StyleHint.Monospace)
        fallbacks = ['JetBrains Mono', 'Cascadia Code', 'Consolas', 'Courier New']
        for fb in fallbacks:
            mono_font.setFamily(fb)

        # Output area
        self.output = QTextEdit()
        self.output.setReadOnly(True)
        self.output.setFont(mono_font)
        self.output.setPlaceholderText("Python console ready. Type commands below.")
        layout.addWidget(self.output, 1)

        # Input line
        self._input_container = QWidget()
        input_layout = QHBoxLayout(self._input_container)
        input_layout.setContentsMargins(8, 3, 8, 3)
        input_layout.setSpacing(6)

        self._prompt_label = QLabel(">>>")
        self._prompt_label.setFont(mono_font)
        input_layout.addWidget(self._prompt_label)

        self.input_line = QLineEdit()
        self.input_line.setFont(mono_font)
        self.input_line.setPlaceholderText("Enter Python command...")
        self.input_line.returnPressed.connect(self._execute)
        input_layout.addWidget(self.input_line, 1)

        layout.addWidget(self._input_container)

        # Apply theme styling to all the widgets that need inline
        # stylesheets (header strip colours, prompt colour, etc.)
        self._apply_theme()

        # Subscribe to runtime theme changes so the console repaints
        # when the user picks a different theme from the menu.
        try:
            theme_signals().theme_changed.connect(self._apply_theme)
        except Exception:
            pass

    def _apply_theme(self):
        """(Re)apply theme-dependent inline stylesheets on the console."""
        self._header.setStyleSheet(
            f"background-color: {COLORS['bg_tertiary']};"
            f"border-bottom: 1px solid {COLORS['border']};"
        )
        self._header_label.setStyleSheet(
            f"color: {COLORS['accent']}; font-size: 11px; font-weight: 600;"
            "background: transparent;"
        )
        self._header_help.setStyleSheet(
            f"color: {COLORS['text_muted']}; font-size: 10px;"
            "background: transparent;"
        )
        self.output.setStyleSheet(
            "QTextEdit {"
            f"  background-color: {COLORS['bg_primary']};"
            f"  color: {COLORS['text_primary']};"
            "   border: none;"
            "   padding: 4px 8px;"
            f"  selection-background-color: {COLORS['accent']};"
            "   selection-color: #FFFFFF;"
            "}"
        )
        self._input_container.setStyleSheet(
            f"background-color: {COLORS['bg_secondary']};"
            f"border-top: 1px solid {COLORS['border']};"
        )
        self._prompt_label.setStyleSheet(
            f"color: {COLORS['accent']}; font-weight: bold;"
            "background: transparent;"
        )
        self.input_line.setStyleSheet(
            "QLineEdit {"
            f"  background-color: {COLORS['bg_widget']};"
            f"  color: {COLORS['text_primary']};"
            f"  border: 1px solid {COLORS['border']};"
            "   border-radius: 6px;"
            "   padding: 6px 10px;"
            "}"
            "QLineEdit:focus {"
            f"  border-color: {COLORS['accent']};"
            "}"
        )

    def _init_namespace(self):
        """Set up the console namespace with useful imports and helpers."""
        self._viewer_3d = None
        self._viewer_2d = None

        self._namespace = {
            '__builtins__': __builtins__,
            'mol': None,
            'print': self._console_print,
            # Main selection command — takes natural language strings
            'sele': self._sele,
            's': self._sele,           # Short alias
            'clear': self._clear_selection,
            'selected': self._get_selected,
            # Direct functions still available
            'select': self._select_by_element,
            'select_idx': self._select_by_indices,
            # Residue labeling
            'label_residues': self._label_residues,
            'clear_labels': self._clear_residue_labels,
            # Measurement
            'distance': self._measure_distance,
            'dist': self._measure_distance,
            'angle': self._measure_angle,
            # Fingerprints
            'morgan_fp': self._morgan_fp,
            'topological_fp': self._topological_fp,
            'maccs_keys': self._maccs_keys,
            'similarity': self._similarity,
            # SASA & Aggregate Mapping
            'count': self._count,
            'sum_prop': self._sum_prop,
            'set_sasa_density': self._set_sasa_density,
            # Help
            'help_cmds': self._show_help,
            'pose': self._pose,
            'docking_pose': self._pose,
        }
        # Pre-import numpy
        try:
            import numpy as np
            self._namespace['np'] = np
        except ImportError:
            pass

    def set_viewer(self, viewer_3d, viewer_2d=None):
        """Set viewer references for selection commands."""
        self._viewer_3d = viewer_3d
        self._viewer_2d = viewer_2d

    def set_molecule(self, molecule):
        """Update the molecule reference in the console namespace."""
        self._namespace['mol'] = molecule

    # ─── Natural Language Selection Parser ─────────────────────────

    def _sele(self, expr):
        """
        Select atoms using a natural expression string.
        Case-insensitive. Examples:
            sele('C')               - all carbon atoms
            sele('N or O')          - nitrogen or oxygen
            sele('C and ring')      - carbons in rings
            sele('not H')           - everything except H
            sele('ring')            - all ring atoms
            sele('nonring')         - all non-ring atoms
            sele('bonded N')        - atoms bonded to nitrogen
            sele('within 3.0 N')    - atoms within 3A of nitrogen
            sele('C and not ring')  - non-ring carbons
            sele('N or O or S')     - multiple element union
        """
        mol = self._namespace.get('mol')
        if not mol:
            self._append_text("No molecule loaded.", COLORS['warning'])
            return []

        try:
            result = self._parse_selection_expr(expr.strip(), mol)
            self._apply_selection(result)
            self._append_output(f"sele '{expr}': {len(result)} atom(s)")
            return sorted(result)
        except Exception as e:
            self._append_text(f"Selection error: {e}", COLORS['error'])
            return []

    def _parse_selection_expr(self, expr, mol):
        """Parse a selection expression into a set of atom indices."""
        import re
        expr = expr.strip()
        if not expr:
            return set()

        # Handle enclosing parentheses
        while expr.startswith('(') and expr.endswith(')'):
            depth = 0
            is_enclosing = True
            for i, c in enumerate(expr):
                if c == '(': depth += 1
                elif c == ')': depth -= 1
                if depth == 0 and i < len(expr) - 1:
                    is_enclosing = False
                    break
            if is_enclosing:
                expr = expr[1:-1].strip()
            else:
                break

        # Pre-process for natural language matching
        # 1. Remove the word 'atoms'
        expr = re.sub(r'(?i)\batoms?\b', '', expr).strip()
        
        # 2. Convert implicit 'and' before 'within' safely using token traversal
        tokens = expr.split()
        new_tokens = []
        for i, t in enumerate(tokens):
            if t.lower() == 'within':
                if i > 0 and tokens[i-1].lower() not in ('and', 'or', 'not', '('):
                    new_tokens.append('and')
            new_tokens.append(t)
        expr = " ".join(new_tokens)

        # Handle 'or' — split on 'or' (case-insensitive)
        # We split on ' or ' to avoid matching within words like 'fluorine'
        parts_or = self._split_keyword(expr, 'or')
        if len(parts_or) > 1:
            result = set()
            for part in parts_or:
                result |= self._parse_selection_expr(part.strip(), mol)
            return result

        # Handle 'and' — split on 'and'
        parts_and = self._split_keyword(expr, 'and')
        if len(parts_and) > 1:
            result = None
            for part in parts_and:
                part_set = self._parse_selection_expr(part.strip(), mol)
                if result is None:
                    result = part_set
                else:
                    result &= part_set
            return result if result is not None else set()

        # Handle 'not' prefix
        lower = expr.lower().strip()
        if lower.startswith('not '):
            inner = expr[4:].strip()
            inner_set = self._parse_selection_expr(inner, mol)
            all_atoms = set(a.index for a in mol.atoms)
            return all_atoms - inner_set

        # Handle 'within X.X <expr>' or 'within N bonds <expr>'
        if lower.startswith('within ') or lower.startswith('near ') or lower.startswith('around '):
            if ' bonds ' in lower or ' bond ' in lower:
                return self._parse_within_bonds(expr, mol)
            return self._parse_within(expr, mol)

        # Handle 'bonded <expr>'
        if lower.startswith('bonded '):
            inner = expr[7:].strip()
            inner_set = self._parse_selection_expr(inner, mol)
            return self._compute_bonded(inner_set, mol)

        # Atomic terms
        return self._parse_atomic_term(lower, mol)

    def _split_keyword(self, expr, keyword):
        """Split expression on keyword (case-insensitive), respecting word boundaries and parenthesis depth."""
        lower = expr.lower()
        parts = []
        pattern = f' {keyword} '
        start = 0
        depth = 0
        i = 0
        
        while i < len(expr):
            if expr[i] == '(':
                depth += 1
            elif expr[i] == ')':
                depth -= 1
                
            if depth == 0 and lower.startswith(pattern, i):
                parts.append(expr[start:i].strip())
                start = i + len(pattern)
                i = start - 1
            i += 1
            
        parts.append(expr[start:].strip())
        return [p for p in parts if p]

    def _parse_within_bonds(self, expr, mol):
        """Parse 'within 5 bonds from ring N'"""
        tokens = expr.split()
        if len(tokens) < 4:
            raise ValueError(f"Expected 'within <N> bonds [from] <expr>', got: {expr}")
        try:
            dist = int(tokens[1])
        except ValueError:
            raise ValueError(f"Invalid bond distance: {tokens[1]}")
        
        start_idx = 3
        if len(tokens) > 3 and tokens[3].lower() in ('from', 'of', 'to'):
            start_idx = 4
            
        inner_expr = " ".join(tokens[start_idx:])
        if not inner_expr:
            raise ValueError(f"Missing target expression in: {expr}")
            
        src = self._parse_selection_expr(inner_expr, mol)
        
        # BFS up to `dist`
        result = set(src)
        current_layer = set(src)
        for _ in range(dist):
            next_layer = set()
            for idx in current_layer:
                next_layer.update(mol.get_neighbors(idx))
            next_layer -= result
            if not next_layer:
                break
            result.update(next_layer)
            current_layer = next_layer
            
        return result

    def _parse_within(self, expr, mol):
        """Parse 'within 3.0 N' or 'near 3.0 C'"""
        import math as m
        tokens = expr.split(None, 2)
        if len(tokens) < 3:
            raise ValueError(f"Expected 'within <dist> <expr>', got: {expr}")
        try:
            dist = float(tokens[1])
        except ValueError:
            raise ValueError(f"Invalid distance: {tokens[1]}")
        inner_expr = tokens[2]
        
        if inner_expr.lower() in ('com', 'center of mass', 'center'):
            com = mol.properties.get('center_of_mass')
            if not com:
                raise ValueError("Center of mass not computed.")
            result = set()
            for a in mol.atoms:
                if not a.has_coords: continue
                dx, dy, dz = a.x - com[0], a.y - com[1], (a.z or 0) - com[2]
                if m.sqrt(dx*dx + dy*dy + dz*dz) <= dist:
                    result.add(a.index)
            return result
            
        src = self._parse_selection_expr(inner_expr, mol)

        result = set()
        for i in src:
            a1 = mol.atoms[i]
            if not a1.has_coords:
                continue
            for a2 in mol.atoms:
                if a2.index in src:
                    continue
                if not a2.has_coords:
                    continue
                dx = a1.x - a2.x
                dy = a1.y - a2.y
                dz = (a1.z or 0) - (a2.z or 0)
                d = m.sqrt(dx*dx + dy*dy + dz*dz)
                if d <= dist:
                    result.add(a2.index)
        return result | src

    def _compute_bonded(self, src_set, mol):
        """Get atoms bonded to source set."""
        result = set()
        for idx in src_set:
            for nb in mol.get_neighbors(idx):
                result.add(nb)
        return result | src_set

    def _parse_atomic_term(self, term, mol):
        """Parse a single selection term."""
        
        # Inequalities: charge <= 0, mass > 12, sasa > 5.0, mlpc < 0
        import re
        ineq_match = re.match(r'^(charge|mass|sasa|mlpc|logpc)\s*(==|<=|>=|<|>|!=)\s*([+-]?\d+(?:\.\d+)?)$', term)
        if ineq_match:
            prop, op, val_str = ineq_match.groups()
            val = float(val_str)
            result = set()
            
            # Lazy lipophilicity calculation
            lipo_vals = None
            if prop in ('mlpc', 'logpc'):
                try:
                    from src.features.cheminformatics.services.lipophilicity_service import LipophilicityService
                    service = LipophilicityService(mol)
                    lipo_vals = service.get_atomic_contributions()
                except Exception:
                    lipo_vals = [0.0] * len(mol.atoms)

            for a in mol.atoms:
                if prop == 'charge': p_val = a.partial_charge or 0.0
                elif prop == 'mass': p_val = a.mass
                elif prop == 'sasa': p_val = getattr(a, 'sasa', 0.0)
                elif prop in ('mlpc', 'logpc'): p_val = lipo_vals[a.index] if lipo_vals else 0.0
                
                match = False
                if op == '==': match = p_val == val
                elif op == '!=': match = p_val != val
                elif op == '>': match = p_val > val
                elif op == '<': match = p_val < val
                elif op == '>=': match = p_val >= val
                elif op == '<=': match = p_val <= val
                
                if match: result.add(a.index)
            return result
            
        # Keywords
        if term in ('ring', 'rings'):
            rings = mol.find_rings()
            result = set()
            for r in rings:
                result.update(r)
            return result

        if term in ('nonring', 'non_ring', 'non-ring', 'chain'):
            rings = mol.find_rings()
            ring_set = set()
            for r in rings:
                ring_set.update(r)
            return set(a.index for a in mol.atoms) - ring_set

        if term in ('all', '*', 'everything'):
            return set(a.index for a in mol.atoms)

        if term in ('none', 'nothing'):
            return set()

        if term in ('heavy', 'hetero'):
            return set(a.index for a in mol.atoms if a.symbol != 'H')

        if term in ('heteroatom', 'heteroatoms'):
            return set(a.index for a in mol.atoms if a.symbol not in ('C', 'H'))

        if term in ('halogen', 'halogens'):
            return set(a.index for a in mol.atoms if a.symbol in ('F', 'Cl', 'Br', 'I', 'At'))

        if term == 'h' or term == 'hydrogen':
            return set(a.index for a in mol.atoms if a.symbol == 'H')
            
        if term == 'hetatm':
            return set(a.index for a in mol.atoms if getattr(a, 'is_hetatm', False))

        if term in ('helix', 'alpha', 'alpha-helix'):
            return set(a.index for a in mol.atoms if getattr(a, 'ss_type', '') == 'H')

        if term in ('sheet', 'beta', 'beta-sheet'):
            return set(a.index for a in mol.atoms if getattr(a, 'ss_type', '') == 'E')

        if term in ('coil', 'loop', 'turn'):
            return set(a.index for a in mol.atoms if getattr(a, 'ss_type', '') in ('C', 'T', 'S'))

        # Chain, resname, resid macro prefixes (e.g. 'chain A')
        if term.startswith('chain '):
            ch = term[6:].strip().upper()
            return set(a.index for a in mol.atoms if getattr(a, 'chain_id', '') == ch)
            
        if term.startswith('resname '):
            rn = term[8:].strip().upper()
            return set(a.index for a in mol.atoms if getattr(a, 'res_name', '').upper() == rn)
            
        if term.startswith('resn '):
            rn = term[5:].strip().upper()
            return set(a.index for a in mol.atoms if getattr(a, 'res_name', '').upper() == rn)
            
        if term.startswith('resid ') or term.startswith('resi ') or term.startswith('resnum '):
            rn = term.split()[1].strip()
            return set(a.index for a in mol.atoms if str(getattr(a, 'res_seq', '')) == rn)

        if term in ('organic', 'ligand', 'ligands'):
            return set(a.index for a in mol.atoms if getattr(a, 'is_hetatm', False) and getattr(a, 'res_name', '').upper() not in ('HOH', 'WAT', 'SOL', 'DOD'))

        if term in ('solvent', 'water', 'sol', 'wat', 'hoh'):
            return set(a.index for a in mol.atoms if getattr(a, 'res_name', '').upper() in ('HOH', 'WAT', 'SOL', 'DOD'))

        # Element symbol — case-insensitive matching
        # Try exact match first, then capitalize
        symbol = term[0].upper() + term[1:].lower() if len(term) > 1 else term.upper()

        # Common element names
        element_names = {
            'carbon': 'C', 'nitrogen': 'N', 'oxygen': 'O', 'sulfur': 'S',
            'phosphorus': 'P', 'fluorine': 'F', 'chlorine': 'Cl', 'bromine': 'Br',
            'iodine': 'I', 'boron': 'B', 'silicon': 'Si', 'selenium': 'Se',
        }
        if term in element_names:
            symbol = element_names[term]

        indices = [a.index for a in mol.atoms if a.symbol == symbol]
        if indices:
            return set(indices)

        # Try as-is (e.g., 'Cl', 'Br')
        indices = [a.index for a in mol.atoms if a.symbol.lower() == term]
        if indices:
            return set(indices)
        
        # Property-based selection shortcuts
        property_terms = {
            'donor': 'donor',
            'acc': 'acceptor', 
            'lipo': 'lipophilic'
        }
        
        if term.lower() in property_terms:
            return self._select_by_property(term)
            
        # Support residue name as shortcut (e.g. 'ALA')
        term_upper = term.upper()
        _AMINO_ACIDS = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'MSE', 'SEC', 'PYL', 'HOH', 'WAT', 'SOL', 'DOD'
        }
        mol_res_names = {getattr(a, 'res_name', '').upper() for a in mol.atoms if getattr(a, 'res_name', None)}
        if term_upper in mol_res_names or term_upper in _AMINO_ACIDS:
            return set(a.index for a in mol.atoms if getattr(a, 'res_name', '').upper() == term_upper)
        
        raise ValueError(f"Unknown selection term: '{term}'")

    # ─── Simple Functions ──────────────────────────────────────────

    def _select_by_element(self, element_symbol):
        """select('C') — select by element (backward compatible)."""
        mol = self._namespace.get('mol')
        if not mol:
            self._append_text("No molecule loaded.", COLORS['warning'])
            return []
        indices = [a.index for a in mol.atoms if a.symbol == element_symbol]
        self._apply_selection(indices)
        self._append_output(f"Selected {len(indices)} {element_symbol} atom(s)")
        return indices

    def _select_by_indices(self, indices):
        """select_idx([0,1,2]) — select by index list."""
        self._apply_selection(list(indices))
        self._append_output(f"Selected {len(indices)} atom(s)")
        return list(indices)

    def _clear_selection(self):
        """Clear all atom selections."""
        self._apply_selection([])
        self._append_output("Selection cleared.")
        return []

    def _select_by_property(self, property_type):
        """Select atoms by chemical property.
        
        Usage:
            sele('donor')   - H-bond donors
            sele('acc')     - H-bond acceptors  
            sele('lipo')    - Lipophilic atoms
        """
        mol = self._namespace.get('mol')
        if not mol:
            self._append_text("No molecule loaded.", COLORS['warning'])
            return []
        
        try:
            from src.features.cheminformatics.services.atom_properties import select_atoms_by_property
            
            # Map short commands to full property names
            property_map = {
                'donor': 'donor',
                'acc': 'acceptor', 
                'lipo': 'lipophilic'
            }
            
            prop_name = property_map.get(property_type.lower(), property_type)
            selected_indices = select_atoms_by_property(mol, prop_name)
            
            self._apply_selection(selected_indices)
            
            # Get property name for display
            display_names = {
                'donor': 'H-bond donors',
                'acceptor': 'H-bond acceptors',
                'lipophilic': 'lipophilic atoms'
            }
            
            display_name = display_names.get(prop_name, prop_name)
            self._append_output(f"Selected {len(selected_indices)} {display_name}")
            return sorted(selected_indices)
            
        except Exception as e:
            self._append_text(f"Property selection error: {e}", COLORS['error'])
            return []

    def _get_selected(self):
        """Return list of currently selected atom indices."""
        if self._viewer_3d:
            return sorted(self._viewer_3d.selected_atoms)
        return []

    # ─── Measurement ───────────────────────────────────────────────

    def _count(self, select_expr='all'):
        """count('charge < 0') — get number of atoms matching selection."""
        mol = self._namespace.get('mol')
        if not mol: return 0
        try:
            res = self._parse_selection_expr(select_expr, mol)
            c = len(res)
            self._append_output(f"Count '{select_expr}': {c}")
            return c
        except Exception as e:
            self._append_text(f"Error: {e}", COLORS['error'])
            return 0

    def _sum_prop(self, prop, select_expr='all'):
        """sum_prop('mass', 'heavy') — sum a property over a selection."""
        mol = self._namespace.get('mol')
        if not mol: return 0.0
        try:
            res = self._parse_selection_expr(select_expr, mol)
            total = 0.0
            for idx in res:
                a = mol.atoms[idx]
                if prop == 'mass': total += a.mass
                elif prop == 'charge': total += a.partial_charge or 0.0
                elif prop == 'sasa': total += getattr(a, 'sasa', 0.0)
                elif prop in ('mlpc', 'logpc'):
                    try:
                        from src.features.cheminformatics.services.lipophilicity_service import LipophilicityService
                        service = LipophilicityService(mol)
                        contribs = service.get_atomic_contributions()
                        total = sum(contribs[idx] for idx in res)
                        break # already summed
                    except Exception: pass
                elif hasattr(a, prop):
                    val = getattr(a, prop)
                    if isinstance(val, (int, float)):
                        total += val
            self._append_output(f"Sum of '{prop}' for '{select_expr}': {total:.4f}")
            return total
        except Exception as e:
            self._append_text(f"Error: {e}", COLORS['error'])
            return 0.0

    def _set_sasa_density(self, density=1600):
        """set_sasa_density(1600) — recompute SASA with higher point density."""
        mol = self._namespace.get('mol')
        if not mol: return
        try:
            from src.features.cheminformatics.services.spatial_properties import compute_sasa
            total_sasa = compute_sasa(mol, n_sphere_points=density)
            if self._viewer_3d:
                self._viewer_3d.update()
            self._append_output(f"SASA recomputed with {density} points/atom. Total SASA: {total_sasa:.2f}")
        except Exception as e:
            self._append_text(f"Error: {e}", COLORS['error'])

    def _measure_distance(self, idx1, idx2):
        """distance(0, 1) — measure distance between two atoms."""
        import math
        mol = self._namespace.get('mol')
        if not mol:
            return 0.0
        a1, a2 = mol.atoms[idx1], mol.atoms[idx2]
        if a1.has_coords and a2.has_coords:
            dx = a1.x - a2.x
            dy = a1.y - a2.y
            dz = (a1.z or 0) - (a2.z or 0)
            d = math.sqrt(dx*dx + dy*dy + dz*dz)
            self._append_output(
                f"Distance {a1.symbol}{idx1}-{a2.symbol}{idx2}: {d:.3f} A")
            return d
        self._append_text("Atoms have no coordinates.", COLORS['warning'])
        return 0.0

    def _measure_angle(self, idx1, idx2, idx3):
        """angle(0, 1, 2) — measure angle at vertex idx2."""
        import math
        mol = self._namespace.get('mol')
        if not mol:
            return 0.0
        a1, a2, a3 = mol.atoms[idx1], mol.atoms[idx2], mol.atoms[idx3]
        if a1.has_coords and a2.has_coords and a3.has_coords:
            import numpy as np
            v1 = np.array([a1.x - a2.x, a1.y - a2.y, (a1.z or 0) - (a2.z or 0)])
            v2 = np.array([a3.x - a2.x, a3.y - a2.y, (a3.z or 0) - (a2.z or 0)])
            cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
            ang = math.degrees(math.acos(max(-1, min(1, cos_a))))
            self._append_output(
                f"Angle {a1.symbol}{idx1}-{a2.symbol}{idx2}-{a3.symbol}{idx3}: {ang:.1f} deg")
            return ang
        self._append_text("Atoms have no coordinates.", COLORS['warning'])
        return 0.0

    # ─── Fingerprints ──────────────────────────────────────────────

    def _morgan_fp(self, radius=2, n_bits=2048):
        """morgan_fp(radius=2, n_bits=2048) — compute Morgan fingerprint."""
        from src.features.cheminformatics.topology.fingerprints import get_morgan_fingerprint
        mol = self._namespace.get('mol')
        if not mol:
            self._append_text("No molecule loaded.", COLORS['warning'])
            return None
        fp = get_morgan_fingerprint(mol, radius, n_bits)
        n_set = sum(fp)
        self._append_output(f"Morgan FP (r={radius}, {n_bits} bits) generated with {n_set} bits set.")
        return fp

    def _topological_fp(self, min_path=1, max_path=7, n_bits=2048):
        """topological_fp(max_path=7, n_bits=2048) — compute topological fingerprint."""
        from src.features.cheminformatics.topology.fingerprints import get_topological_fingerprint
        mol = self._namespace.get('mol')
        if not mol:
            self._append_text("No molecule loaded.", COLORS['warning'])
            return None
        fp = get_topological_fingerprint(mol, min_path, max_path, n_bits)
        n_set = sum(fp)
        self._append_output(f"Topological FP (path={min_path}-{max_path}, {n_bits} bits) generated with {n_set} bits set.")
        return fp

    def _maccs_keys(self):
        """maccs_keys() — compute MACCS keys."""
        from src.features.cheminformatics.topology.fingerprints import get_maccs_keys
        mol = self._namespace.get('mol')
        if not mol:
            self._append_text("No molecule loaded.", COLORS['warning'])
            return None
        fp = get_maccs_keys(mol)
        n_set = sum(fp)
        self._append_output(f"MACCS Keys (166 bits) generated with {n_set} bits set.")
        return fp

    def _similarity(self, fp1, fp2):
        """similarity(fp1, fp2) — compute Tanimoto similarity."""
        from src.features.cheminformatics.topology.fingerprints import tanimoto_similarity
        if not fp1 or not fp2:
            self._append_text("Invalid fingerprints provided.", COLORS['warning'])
            return 0.0
        try:
            sim = tanimoto_similarity(fp1, fp2)
            self._append_output(f"Tanimoto Similarity: {sim:.4f}")
            return sim
        except Exception as e:
            self._append_text(f"Error computing similarity: {e}", COLORS['warning'])
            return 0.0

    def _pose(self, dist=5.0, show_hbonds=True, show_salt=True, show_hydro=False):
        """pose(dist=5.0) — Generate and display 3D docking pose."""
        mol = self._namespace.get('mol')
        if not mol or not self._viewer_3d:
            self._append_text("No molecule or 3D viewer available.", COLORS['warning'])
            return
            
        service = DockingPoseService(mol)
        ligands = service.find_ligands()
        if not ligands:
            self._append_text("No ligands detected.", COLORS['warning'])
            return
            
        ligand_indices = ligands[0] # Use largest ligand by default
        nearby_res = service.find_nearby_residues(ligand_indices, dist)
        interactions = service.detect_interactions(ligand_indices, nearby_res)
        
        # Update viewer
        v = self._viewer_3d
        v.show_ligands_in_cartoon = True  # Ensure ligands are visible
        v.visible_sidechains = nearby_res
        
        interaction_colors = {
            "Hydrogen Bond": "#00FF00",
            "Salt Bridge": "#FF00FF",
            "Hydrophobic": "#FFFF00",
            "Pi-Stacking": "#00FFFF"
        }
        
        lines = []
        for inter in interactions:
            if inter.type == "Hydrogen Bond" and not show_hbonds: continue
            if inter.type == "Salt Bridge" and not show_salt: continue
            if inter.type == "Hydrophobic" and not show_hydro: continue
            color = interaction_colors.get(inter.type, "#FFFFFF")
            lines.append((inter.atom1_idx, inter.atom2_idx, inter.type, color))
        
        v.interaction_lines = lines
        
        cam = {}
        labels = {}
        processed = set()
        for idx in ligand_indices:
            cam[idx] = 'ball_and_stick'
        for atom in mol.atoms:
            if atom.res_seq in nearby_res:
                pdb_name = getattr(atom, 'pdb_name', '').strip()
                if pdb_name not in ('CA', 'C', 'N', 'O'):
                    cam[atom.index] = 'stick'
                # Disable automatic residue labeling to reduce visual clutter
                # if pdb_name == 'CA' and atom.res_seq not in processed:
                #     labels[atom.index] = f"{getattr(atom, 'res_name', 'UNK')}{atom.res_seq}"
                #     processed.add(atom.res_seq)
        
        v.custom_atom_modes = cam
        v.labels = labels
        v.update()
        
        self._append_output(f"Docking pose generated: {len(ligand_indices)} ligand atoms, "
                            f"{len(nearby_res)} residues, {len(interactions)} interactions.")
        return interactions

    # ─── Internals ─────────────────────────────────────────────────

    def _apply_selection(self, indices):
        """Apply selection to viewers."""
        idx_set = set(indices)
        if self._viewer_3d:
            self._viewer_3d.set_selected(idx_set)
        if self._viewer_2d:
            self._viewer_2d.set_selected(idx_set)

    def _show_help(self):
        """Show available console commands."""
        help_text = (
            "=== Selection (natural syntax) ===\n"
            "  sele('C')              Select carbons\n"
            "  sele('N or O')         Nitrogen or oxygen\n"
            "  sele('C and ring')     Carbons in rings\n"
            "  sele('not H')          All except H\n"
            "  sele('ring')           All ring atoms\n"
            "  sele('nonring')        Non-ring atoms\n"
            "  sele('bonded N')       Bonded to nitrogen\n"
            "  sele('within 3.0 N')   Within 3A of nitrogen\n"
            "  sele('within 2 bonds O') Within 2 bonds of oxygen\n"
            "  sele('C within 5 bonds N') Carbons within 5 bonds of N\n"
            "  sele('C and not ring') Non-ring carbons\n"
            "  sele('charge < 0')     Negative partial charge\n"
            "  sele('mass > 12')      Mass greater than 12 Da\n"
            "  sele('sasa > 5.0')     Solvent exposed atoms\n"
            "  sele('mlpc < 0')       Hydrophilic atoms (numeric)\n"
            "  sele('within 5.0 COM') Within 5A of Center of Mass\n"
            "  sele('heavy')          Non-hydrogen atoms\n"
            "  sele('halogen')        F, Cl, Br, I\n"
            "  sele('heteroatom')     Non-C, Non-H atoms\n"
            "  sele('chain A')        PDB chain A\n"
            "  sele('resname ALA')    PDB residue name\n"
            "  sele('resnum 10')      PDB residue number\n"
            "  sele('organic')        Ligands (non-water hetatms)\n"
            "  sele('solvent')        Water molecules (HOH/WAT)\n"
            "  s('N or O or S')       s() = short for sele()\n"
            "  clear()                Clear selection\n"
            "  label_residues()       Label visible residues\n"
            "  label_residues('chain A') Label residues in selection\n"
            "  clear_labels()         Clear all residue labels\n"
            "  selected()             Get selected indices\n"
            "  count('expr')          Number of atoms in selection\n"
            "  sum_prop('sasa', 'C')  Sum property across selection\n"
            "  set_sasa_density(500)  Recompute SASA with N points\n"
            "=== Measurement ===\n"
            "  distance(0, 1)         Distance in Angstroms\n"
            "  angle(0, 1, 2)         Angle in degrees\n"
            "=== Fingerprints ===\n"
            "  fp = morgan_fp(r=2)    Morgan (ECFP) fingerprint\n"
            "  fp = topological_fp()  Daylight-like topological\n"
            "  fp = maccs_keys()      MACCS structural keys (subset)\n"
            "  similarity(fp1, fp2)   Tanimoto similarity\n"
            "=== Molecule ===\n"
            "  mol.atoms / bonds      Atom/bond lists\n"
            "  mol.molecular_formula()\n"
            "  np                     NumPy module"
        )
        self._append_text(help_text, COLORS['accent2'])

    def _console_print(self, *args, **kwargs):
        """Custom print that outputs to the console widget."""
        output = io.StringIO()
        kwargs['file'] = output
        print(*args, **kwargs)
        self._append_output(output.getvalue().rstrip('\n'))

    def _execute(self):
        """Execute the current input line."""
        cmd = self.input_line.text().strip()
        if not cmd:
            return

        # Show command
        self._append_text(f">>> {cmd}", COLORS['accent2'])

        # Add to history
        self._history.append(cmd)
        self._history_idx = len(self._history)
        self.input_line.clear()

        # Execute
        try:
            # Try eval first (expression)
            try:
                result = eval(cmd, self._namespace)
                if result is not None:
                    self._append_output(repr(result))
            except SyntaxError:
                # Fall back to exec (statement)
                exec(cmd, self._namespace)
        except Exception:
            tb = traceback.format_exc()
            # Show only the last few lines
            lines = tb.strip().split('\n')
            if len(lines) > 4:
                lines = lines[-3:]
            self._append_text('\n'.join(lines), COLORS['error'])

        self.command_executed.emit(cmd)

    def _append_output(self, text):
        """Append normal output."""
        self._append_text(text, COLORS['text_primary'])

    def _append_text(self, text, color):
        """Append colored text to output."""
        cursor = self.output.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)

        fmt = cursor.charFormat()
        fmt.setForeground(QColor(color))
        cursor.setCharFormat(fmt)
        cursor.insertText(text + '\n')

        self.output.setTextCursor(cursor)
        self.output.ensureCursorVisible()

    def _label_residues(self, selection_expr=None):
        """
        Label residues for selected atoms or visible residues.
        
        Usage:
            label_residues()                    - Label all visible residues
            label_residues('within 5.0 ligand') - Label residues matching selection
            label_residues('chain A')           - Label residues in chain A
        """
        mol = self._namespace.get('mol')
        v = self._viewer_3d
        if not mol or not v:
            self._append_text("No molecule or viewer available.", COLORS['warning'])
            return []
        
        try:
            # Determine which residues to label
            if selection_expr:
                # Use selection expression to find residues
                selected_atoms = self._parse_selection_expr(selection_expr, mol)
                residue_seqs = set()
                for atom_idx in selected_atoms:
                    atom = mol.atoms[atom_idx]
                    rs = getattr(atom, 'res_seq', None)
                    if rs is not None:
                        residue_seqs.add(rs)
            else:
                # Use currently visible residues
                residue_seqs = v.visible_sidechains.copy()
                
                # Also include residues with custom atom modes (like ligands)
                if hasattr(v, 'custom_atom_modes'):
                    for atom_idx in v.custom_atom_modes:
                        atom = mol.atoms[atom_idx]
                        rs = getattr(atom, 'res_seq', None)
                        if rs is not None:
                            residue_seqs.add(rs)
            
            # Create labels for CA atoms of selected residues
            labels_count = 0
            if not hasattr(v, 'labeled_residues'):
                v.labeled_residues = {}
                
            for atom in mol.atoms:
                rs = getattr(atom, 'res_seq', None)
                if rs is not None and rs in residue_seqs:
                    if hasattr(atom, 'pdb_name') and atom.pdb_name.strip() == 'CA':
                        if rs not in v.labeled_residues:
                            v.labeled_residues[rs] = QColor(Qt.GlobalColor.black)
                            labels_count += 1
            
            v.update()
            self._append_output(f"Labeled {labels_count} residues.")
            return list(residue_seqs)
            
        except Exception as e:
            self._append_text(f"Error labeling residues: {e}", COLORS['error'])
            return []

    def _clear_residue_labels(self):
        """Clear all residue labels from the viewer."""
        v = self._viewer_3d
        if not v:
            self._append_text("No viewer available.", COLORS['warning'])
            return
        
        # Clear all labels
        v.labels.clear()
        v.labeled_residues.clear()
        v.update()
        self._append_output("All residue labels cleared.")

    def keyPressEvent(self, event):
        """Handle history navigation."""
        if event.key() == Qt.Key.Key_Up and self._history:
            self._history_idx = max(0, self._history_idx - 1)
            self.input_line.setText(self._history[self._history_idx])
        elif event.key() == Qt.Key.Key_Down and self._history:
            self._history_idx = min(len(self._history), self._history_idx + 1)
            if self._history_idx < len(self._history):
                self.input_line.setText(self._history[self._history_idx])
            else:
                self.input_line.clear()
        else:
            super().keyPressEvent(event)
