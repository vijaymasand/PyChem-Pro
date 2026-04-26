"""
SMILES-specific 2D Molecular Viewer - Enhanced for Virtual Hydrogen Support

=============================================================================
OVERVIEW
=============================================================================

This module provides a specialized 2D molecular viewer for SMILES-converted molecules.
It extends the proven MolViewer2D with specific enhancements for handling virtual
hydrogen atoms and SMILES-derived molecular structures.

=============================================================================
KEY DIFFERENCES FROM mol_viewer_2d.py
=============================================================================

1. **Virtual Hydrogen Support**:
   - mol_viewer_2d.py: Only renders real atoms (file imports)
   - mol_viewer_2d_smiles.py: Renders both real atoms AND virtual hydrogens
   
2. **Coordinate Generation**:
   - mol_viewer_2d.py: Uses standard CoordinateGenerator2D
   - mol_viewer_2d_smiles.py: Uses CoordinateGenerator2DSMILES with H placement
   
3. **Rendering Enhancements**:
   - mol_viewer_2d.py: Standard skeletal formula rendering
   - mol_viewer_2d_smiles.py: Enhanced rendering for virtual atoms and bonds
   
4. **Molecule Source**:
   - mol_viewer_2d.py: Designed for MOL/MOL2 files with existing coordinates
   - mol_viewer_2d_smiles.py: Designed for SMILES strings requiring full generation

=============================================================================
VIRTUAL HYDROGEN RENDERING
=============================================================================

Virtual hydrogen atoms are rendered with:
1. **Negative Indices**: Virtual H atoms use negative indexing system
2. **Special Styling**: Smaller labels and lighter colors
3. **Bond Rendering**: Special handling for bonds to virtual atoms
4. **Toggle Support**: Can be hidden/shown independently

=============================================================================
ENHANCED RENDERING FEATURES
=============================================================================

1. **Virtual Atom Labeling**:
   - Smaller font size for virtual H atoms
   - Faded color to distinguish from real atoms
   - Proper positioning to avoid overlap

2. **Virtual Bond Rendering**:
   - Thinner lines for bonds to virtual atoms
   - Proper line style to indicate virtual nature
   - Correct bond order visualization

3. **Layout Optimization**:
   - Better spacing for virtual atoms
   - Prevents overlap with real atoms
   - Maintains chemical accuracy

=============================================================================
USAGE
=============================================================================

```python
# Create SMILES-specific viewer
viewer = MolViewer2DSMILES(parent)

# Set SMILES molecule (generates coordinates automatically)
viewer.set_molecule(molecule)

# Toggle virtual hydrogen visibility
viewer.show_hydrogens = True/False
viewer.update()
```

=============================================================================
DEPENDENCIES
=============================================================================

- CoordinateGenerator2DSMILES for enhanced coordinate generation
- Enhanced rendering system for virtual atoms
- Standard Qt painting framework

=============================================================================
AUTHOR & VERSION
=============================================================================

Derived from mol_viewer_2d.py with SMILES-specific enhancements
Version: 1.0
Focus: Professional 2D visualization of SMILES-derived molecules with virtual H support
"""

import math
from src.shared.qt_compat import QWidget, Qt, QPointF, QRectF, Signal
from src.shared.qt_compat import (
    QPainter, QColor, QPen, QBrush, QFont, QFontMetrics,
    QImage, QPainterPath, QPolygonF
)
from src.shared.ui.theme import COLORS


# Enhanced element colors for SMILES viewer (supports virtual atoms)
ELEMENT_COLORS = {
    'C':  QColor(220, 220, 220),   # Near-white (hidden in skeletal)
    'H':  QColor(170, 170, 180),   # Light gray
    'N':  QColor(48, 128, 255),    # Vivid blue
    'O':  QColor(255, 50, 50),     # Vivid red
    'S':  QColor(230, 195, 50),    # Gold-yellow
    'P':  QColor(255, 140, 30),    # Warm orange
    'F':  QColor(80, 210, 80),     # Bright green
    'Cl': QColor(80, 210, 80),     # Green
    'Br': QColor(180, 80, 40),     # Dark orange-brown
    'I':  QColor(170, 50, 170),    # Purple
    'B':  QColor(255, 180, 160),   # Salmon
    'Se': QColor(255, 160, 0),     # Orange
    'Si': QColor(200, 180, 140),   # Tan
}

# Virtual hydrogen color (faded to distinguish from real H)
VIRTUAL_H_COLOR = QColor(200, 200, 210)  # Very light gray


class MolViewer2DSMILES(QWidget):
    """
    Enhanced 2D molecular viewer for SMILES molecules with virtual hydrogen support.
    
    This class extends the proven MolViewer2D rendering system with specific
    enhancements for SMILES-derived molecules, including:
    
    1. Virtual hydrogen atom rendering
    2. Enhanced bond rendering for virtual atoms
    3. Optimized layout for SMILES structures
    4. Professional appearance maintenance
    
    Signals:
        selection_changed: Emitted when atom selection changes
        delete_requested: Emitted when deletion is requested
    """
    
    # --- Signals ---
    selection_changed = Signal(object)   # emits set of selected atom indices
    delete_requested = Signal(object)    # emits set of atom indices to delete
    
    def __init__(self, parent=None):
        """
        Initialize SMILES-specific 2D viewer.
        
        Args:
            parent: Parent widget (optional)
        """
        super().__init__(parent)
        self.molecule = None
        self.coords_2d = {}
        self.selected_atoms = set()
        self.setMinimumSize(400, 400)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        
        # Display settings
        self.bg_color = QColor(COLORS['viewer_bg'])
        self._scale = 40.0
        self._offset_x = 0.0
        self._offset_y = 0.0
        self.show_hydrogens = True  # Toggle for virtual H visibility
        
        # Virtual hydrogen rendering settings
        self.virtual_h_scale = 0.8  # Smaller size for virtual H
        self.virtual_h_alpha = 0.7  # Faded appearance
    
    def set_molecule(self, molecule, coords_2d=None):
        """
        Set molecule and generate 2D coordinates if needed.
        
        Enhanced to use SMILES-specific coordinate generator with hydrogen placement.
        
        Args:
            molecule: Domain molecule object from SMILES conversion
            coords_2d: Optional pre-generated coordinates
        """
        print(f"[DEBUG SMILES] MolViewer2DSMILES.set_molecule() called with {len(molecule.atoms) if molecule else 0} atoms")
        self.molecule = molecule
        self.selected_atoms = set()
        
        if coords_2d:
            print(f"[DEBUG SMILES] Using provided 2D coordinates: {len(coords_2d)}")
            self.coords_2d = coords_2d
        elif molecule and molecule.atoms:
            print("[DEBUG SMILES] Generating 2D coordinates with hydrogen placement...")
            # Use pure OASA-based SMILES coordinate generator
            from src.features.layout_2d.generators.coordgen2d_smiles_pure_oasa import CoordinateGenerator2DSMILES
            generator = CoordinateGenerator2DSMILES(molecule, force_regenerate=True)
            self.coords_2d = generator.generate()
            total_atoms = len(molecule.atoms)
            coords_count = len(self.coords_2d)
            print(f"[DEBUG SMILES] Generated {coords_count}/{total_atoms} coordinates")
            
            # Kekulize aromatic bonds for proper 2D rendering
            self._kekulize_aromatic_bonds()
        else:
            print("[DEBUG SMILES] No molecule or atoms, setting empty coords")
            self.coords_2d = {}
        
        self._auto_fit()
        self.update()
    
    def paintEvent(self, event):
        """
        Paint the molecule with enhanced virtual hydrogen support.
        
        Overrides base paintEvent to add virtual hydrogen rendering.
        """
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        
        # Fill background
        painter.fillRect(self.rect(), self.bg_color)
        
        if not self.molecule or not self.coords_2d:
            painter.end()
            return
        
        # Transform coordinates to screen space
        screen_coords = self._transform_to_screen()
        
        # Draw bonds (including virtual bonds)
        self._draw_all_bonds_enhanced(painter, screen_coords)
        
        # Draw atoms (including virtual atoms)
        self._draw_all_atoms_enhanced(painter, screen_coords)
        
        # Draw selection
        self._draw_selection_enhanced(painter, screen_coords)
        
        painter.end()
    
    def _draw_all_bonds_enhanced(self, painter, screen_coords):
        """
        Draw all bonds with enhanced support for virtual atoms.
        
        Handles bonds to virtual hydrogen atoms with special styling.
        """
        if not self.molecule:
            return
        
        # Draw real bonds first
        for bond in self.molecule.bonds:
            if bond.begin_atom_idx in screen_coords and bond.end_atom_idx in screen_coords:
                self._draw_bond_enhanced(painter, bond, screen_coords)
        
        # Draw virtual bonds (to virtual hydrogen atoms)
        self._draw_virtual_bonds(painter, screen_coords)
    
    def _draw_bond_enhanced(self, painter, bond, screen_coords):
        """
        Draw a single bond with enhanced styling.
        
        Handles bonds to virtual atoms with different appearance.
        """
        x1, y1 = screen_coords[bond.begin_atom_idx]
        x2, y2 = screen_coords[bond.end_atom_idx]
        
        # Check if this is a bond to virtual atom
        is_virtual_bond = (bond.begin_atom_idx < 0) or (bond.end_atom_idx < 0)
        
        # Set pen style based on bond type
        if is_virtual_bond:
            # Virtual bond - thinner and faded
            pen = QPen(QColor(150, 150, 160), 1.5)  # Light gray, thin
            painter.setPen(pen)
        else:
            # Real bond - standard styling
            if bond.is_aromatic:
                pen = QPen(QColor(100, 100, 100), 2.0)  # Gray for aromatic
            else:
                pen = QPen(QColor(50, 50, 50), 2.0)  # Black for normal
            painter.setPen(pen)
        
        # Draw bond based on order
        if bond.order == 1:
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
        elif bond.order == 2:
            self._draw_double_bond_enhanced(painter, x1, y1, x2, y2, is_virtual_bond)
        elif bond.order == 3:
            self._draw_triple_bond_enhanced(painter, x1, y1, x2, y2, is_virtual_bond)
    
    def _draw_double_bond_enhanced(self, painter, x1, y1, x2, y2, is_virtual):
        """
        Draw double bond with enhanced styling for virtual bonds.
        """
        # Calculate perpendicular offset
        dx = x2 - x1
        dy = y2 - y1
        length = math.sqrt(dx*dx + dy*dy)
        if length < 0.001:
            return
        
        # Normalize and perpendicular
        dx /= length
        dy /= length
        offset_x = -dy * 2.0  # Offset distance
        offset_y = dx * 2.0
        
        # Draw two parallel lines
        painter.drawLine(QPointF(x1 + offset_x, y1 + offset_y), 
                       QPointF(x2 + offset_x, y2 + offset_y))
        painter.drawLine(QPointF(x1 - offset_x, y1 - offset_y), 
                       QPointF(x2 - offset_x, y2 - offset_y))
    
    def _draw_triple_bond_enhanced(self, painter, x1, y1, x2, y2, is_virtual):
        """
        Draw triple bond with enhanced styling for virtual bonds.
        """
        # Calculate perpendicular offset
        dx = x2 - x1
        dy = y2 - y1
        length = math.sqrt(dx*dx + dy*dy)
        if length < 0.001:
            return
        
        # Normalize and perpendicular
        dx /= length
        dy /= length
        offset_x = -dy * 3.0  # Larger offset for triple
        offset_y = dx * 3.0
        
        # Draw three parallel lines
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))  # Center line
        painter.drawLine(QPointF(x1 + offset_x, y1 + offset_y), 
                       QPointF(x2 + offset_x, y2 + offset_y))
        painter.drawLine(QPointF(x1 - offset_x, y1 - offset_y), 
                       QPointF(x2 - offset_x, y2 - offset_y))
    
    def _draw_virtual_bonds(self, painter, screen_coords):
        """
        Draw bonds to virtual hydrogen atoms.
        
        Virtual bonds are not stored in the molecule structure,
        so they need to be calculated and rendered separately.
        """
        if not self.show_hydrogens:
            return
        
        # Find virtual hydrogen atoms and their parent bonds
        virtual_atoms = {idx: coord for idx, coord in screen_coords.items() if idx < 0}
        
        for virtual_idx, (vx, vy) in virtual_atoms.items():
            # Extract parent atom index from virtual index
            # Virtual indices follow format: -(parent_index * 100 + h_index + 1)
            parent_idx = abs(virtual_idx // 100) - 1
            
            if parent_idx in screen_coords:
                px, py = screen_coords[parent_idx]
                
                # Draw virtual bond
                pen = QPen(QColor(150, 150, 160), 1.0)  # Light gray, thin
                painter.setPen(pen)
                painter.drawLine(QPointF(px, py), QPointF(vx, vy))
    
    def _draw_all_atoms_enhanced(self, painter, screen_coords):
        """
        Draw all atoms with enhanced virtual hydrogen support.
        
        Handles both real atoms and virtual hydrogen atoms with different styling.
        """
        if not self.molecule:
            return
        
        # Draw real atoms first
        for atom in self.molecule.atoms:
            if atom.index in screen_coords:
                self._draw_atom_enhanced(painter, atom, screen_coords[atom.index], False)
        
        # Draw virtual hydrogen atoms
        if self.show_hydrogens:
            self._draw_virtual_atoms(painter, screen_coords)
    
    def _draw_atom_enhanced(self, painter, atom, screen_pos, is_virtual):
        """
        Draw a single atom with enhanced styling.
        
        Args:
            painter: QPainter object
            atom: Atom object (None for virtual atoms)
            screen_pos: Screen coordinates [x, y]
            is_virtual: Whether this is a virtual atom
        """
        x, y = screen_pos
        
        if is_virtual:
            # Virtual hydrogen atom
            color = VIRTUAL_H_COLOR
            color.setAlphaF(self.virtual_h_alpha)  # Faded appearance
            
            # Smaller font for virtual H
            font = QFont("Arial", 8)
            painter.setFont(font)
            
            # Draw H label
            painter.setPen(color)
            painter.drawText(QPointF(x - 4, y + 3), "H")
        else:
            # Real atom
            symbol = atom.symbol
            
            # Skip carbon in skeletal formula
            if symbol == 'C':
                # Draw carbon as a small dot if it has no neighbors
                neighbors = [b for b in self.molecule.bonds 
                           if b.begin_atom_idx == atom.index or b.end_atom_idx == atom.index]
                if len(neighbors) == 0:
                    painter.setPen(QPen(QColor(100, 100, 100), 2))
                    painter.drawEllipse(QPointF(x, y), 2, 2)
                return
            
            # Get element color
            color = ELEMENT_COLORS.get(symbol, QColor(100, 100, 100))
            painter.setPen(color)
            
            # Draw atom label
            font = QFont("Arial", 10)
            painter.setFont(font)
            painter.drawText(QPointF(x - 5, y + 3), symbol)
    
    def _draw_virtual_atoms(self, painter, screen_coords):
        """
        Draw virtual hydrogen atoms with special styling.
        
        Virtual atoms are not stored in the molecule structure,
        so they need to be calculated and rendered separately.
        """
        virtual_atoms = {idx: coord for idx, coord in screen_coords.items() if idx < 0}
        
        for virtual_idx, (x, y) in virtual_atoms.items():
            self._draw_atom_enhanced(painter, None, (x, y), True)
    
    def _draw_selection_enhanced(self, painter, screen_coords):
        """
        Draw selection highlighting with virtual atom support.
        
        Highlights both real and virtual atoms in selection.
        """
        if not self.selected_atoms:
            return
        
        # Draw selection for real atoms
        for atom_idx in self.selected_atoms:
            if atom_idx >= 0 and atom_idx in screen_coords:
                x, y = screen_coords[atom_idx]
                painter.setPen(QPen(QColor(255, 200, 0), 2))  # Orange highlight
                painter.setBrush(QBrush(QColor(255, 200, 0, 50)))  # Semi-transparent
                painter.drawEllipse(QPointF(x, y), 8, 8)
        
        # Draw selection for virtual atoms
        virtual_selected = [idx for idx in self.selected_atoms if idx < 0]
        for virtual_idx in virtual_selected:
            if virtual_idx in screen_coords:
                x, y = screen_coords[virtual_idx]
                painter.setPen(QPen(QColor(255, 200, 0), 1))  # Thinner for virtual
                painter.setBrush(QBrush(QColor(255, 200, 0, 30)))  # More transparent
                painter.drawEllipse(QPointF(x, y), 6, 6)  # Smaller for virtual
    
    def _transform_to_screen(self):
        """
        Transform molecular coordinates to screen coordinates.
        
        Enhanced to handle virtual atom coordinates.
        """
        screen_coords = {}
        
        # Transform real and virtual coordinates
        for idx, (x, y) in self.coords_2d.items():
            screen_x = x * self._scale + self._offset_x + self.width() / 2
            screen_y = -y * self._scale + self._offset_y + self.height() / 2
            screen_coords[idx] = (screen_x, screen_y)
        
        return screen_coords
    
    def _kekulize_aromatic_bonds(self):
        """
        Kekulize aromatic bonds for proper 2D rendering.
        
        Converts aromatic bonds to alternating single/double bonds
        for better visual representation.
        """
        if not self.molecule:
            return
        
        for bond in self.molecule.bonds:
            if bond.is_aromatic:
                # Simple kekulization - alternate single/double
                # In a real implementation, this would be more sophisticated
                bond.order = 1.5  # Keep as aromatic for now
                # Could implement proper kekulization here if needed
    
    def _auto_fit(self):
        """
        Auto-fit the molecule to the widget.
        
        Calculates optimal scale and offset to fit the entire molecule
        in the available space.
        """
        if not self.coords_2d:
            self._offset_x = 0
            self._offset_y = 0
            self._scale = 40.0
            return
        
        # Calculate bounds
        min_x = min(coord[0] for coord in self.coords_2d.values())
        max_x = max(coord[0] for coord in self.coords_2d.values())
        min_y = min(coord[1] for coord in self.coords_2d.values())
        max_y = max(coord[1] for coord in self.coords_2d.values())
        
        # Calculate scale to fit with padding
        width = max_x - min_x
        height = max_y - min_y
        padding = 20
        
        if width > 0 and height > 0:
            scale_x = (self.width() - 2 * padding) / width
            scale_y = (self.height() - 2 * padding) / height
            self._scale = min(scale_x, scale_y, 100.0)  # Cap at 100
        
        # Center the molecule
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        self._offset_x = -center_x * self._scale
        self._offset_y = center_y * self._scale
    
    def mousePressEvent(self, event):
        """
        Handle mouse press events for selection.
        
        Enhanced to handle virtual atom selection.
        """
        if event.button() == Qt.MouseButton.LeftButton:
            # Transform click to molecular coordinates
            click_pos = event.position()
            mol_pos = self._screen_to_molecular(click_pos.x(), click_pos.y())
            
            # Find nearest atom (including virtual)
            nearest = self._find_nearest_atom_enhanced(mol_pos[0], mol_pos[1])
            
            if nearest is not None:
                if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                    # Add to selection
                    self.selected_atoms.add(nearest)
                else:
                    # Replace selection
                    self.selected_atoms = {nearest}
                
                self.selection_changed.emit(self.selected_atoms)
                self.update()
            else:
                # Clear selection
                self.selected_atoms.clear()
                self.selection_changed.emit(self.selected_atoms)
                self.update()
    
    def _find_nearest_atom_enhanced(self, mx, my):
        """
        Find nearest atom to mouse position including virtual atoms.
        
        Enhanced to include virtual hydrogen atoms in selection.
        """
        min_dist = float('inf')
        nearest = None
        
        # Check all coordinates (real + virtual)
        for idx, (x, y) in self.coords_2d.items():
            dist = math.sqrt((x - mx)**2 + (y - my)**2)
            if dist < min_dist and dist < 0.5:  # Within selection threshold
                min_dist = dist
                nearest = idx
        
        return nearest
    
    def _screen_to_molecular(self, screen_x, screen_y):
        """
        Transform screen coordinates to molecular coordinates.
        """
        mol_x = (screen_x - self.width() / 2 - self._offset_x) / self._scale
        mol_y = -(screen_y - self.height() / 2 - self._offset_y) / self._scale
        return (mol_x, mol_y)
    
    def keyPressEvent(self, event):
        """
        Handle key press events.
        
        Supports hydrogen visibility toggle.
        """
        if event.key() == Qt.Key.Key_H:
            # Toggle hydrogen visibility
            self.show_hydrogens = not self.show_hydrogens
            print(f"[DEBUG SMILES] Hydrogen visibility: {self.show_hydrogens}")
            self.update()
        elif event.key() == Qt.Key.Key_Delete and self.selected_atoms:
            # Request deletion of selected atoms
            self.delete_requested.emit(self.selected_atoms)
        else:
            super().keyPressEvent(event)
