"""
GUI-based color selection interface similar to PyMol.

Provides visual color selection for atoms, spheres, and bonds.
"""

from typing import Dict, List, Optional, Tuple
from src.shared.ui.theme import COLORS


class ColorSelectionGUI:
    """
    GUI-based color selection interface.
    """
    
    def __init__(self):
        self.color_presets = self._create_color_presets()
        self.selected_colors = {}
        
    def _create_color_presets(self) -> Dict[str, str]:
        """Create color presets similar to PyMol."""
        return {
            'red': '#ff0000',
            'green': '#00ff00',
            'blue': '#0000ff',
            'yellow': '#ffff00',
            'orange': '#ff8000',
            'purple': '#8000ff',
            'pink': '#ff00ff',
            'cyan': '#00ffff',
            'white': '#ffffff',
            'black': '#000000',
            'gray': '#808080',
            'brown': '#8b4513',
            'lime': '#00ff00',
            'navy': '#000080',
            'teal': '#008080',
            'maroon': '#800000',
            'olive': '#808000',
            'silver': '#c0c0c0',
            'gold': '#ffd700',
            'indigo': '#4b0082'
        }
    
    def show_color_menu(self) -> Dict[str, str]:
        """Show color selection menu."""
        print("\n" + "="*60)
        print("SMILES Molecular Toolkit - Color Selection")
        print("="*60)
        
        print("\n" + "-"*40)
        print("COLOR SELECTION OPTIONS")
        print("-"*40)
        print("1. Atom Colors - Set individual atom colors")
        print("2. Sphere Colors - Set sphere and dummy atom colors")
        print("3. Bond Colors - Set bond and stick colors")
        print("4. Show Current Colors - Display current settings")
        print("5. Apply Colors - Apply and exit")
        print("0. Cancel")
        
        while True:
            try:
                choice = input("\nEnter your choice (0-5): ").strip()
                
                if choice == '0':
                    print("Color selection cancelled.")
                    return {}
                elif choice == '1':
                    self._set_atom_colors()
                elif choice == '2':
                    self._set_sphere_colors()
                elif choice == '3':
                    self._set_bond_colors()
                elif choice == '4':
                    self._show_current_colors()
                elif choice == '5':
                    return self.selected_colors
                else:
                    print("Invalid choice. Please enter 0-5.")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input. Please try again.")
                break
        
        return {}
    
    def _set_atom_colors(self):
        """Set colors for individual atoms."""
        print("\n" + "-"*40)
        print("ATOM COLOR SELECTION")
        print("-"*40)
        
        print("\nAvailable atoms: C, H, O, N, S, P, F, Cl, Br, I")
        print("Available colors:")
        
        # Display color options
        color_list = list(self.color_presets.items())
        for i, (name, hex_code) in enumerate(color_list, 1):
            print(f"{i:2d}. {name:8s} - {hex_code}")
        
        print("\nFormat: atom_name color_name")
        print("Examples: C red, O blue, N green")
        print("Type 'done' when finished")
        
        while True:
            try:
                selection = input("Atom color (or 'done'): ").strip()
                if selection.lower() == 'done':
                    break
                
                if ' ' in selection:
                    atom, color = selection.split(' ', 1)
                    atom = atom.strip().upper()
                    color = color.strip().lower()
                    
                    # Validate atom
                    valid_atoms = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']
                    if atom not in valid_atoms:
                        print(f"Invalid atom: {atom}. Valid atoms: {', '.join(valid_atoms)}")
                        continue
                    
                    # Validate color
                    if color not in self.color_presets:
                        print(f"Invalid color: {color}")
                        continue
                    
                    # Set the color
                    color_key = f'atom_{atom.lower()}'
                    color_value = self.color_presets[color]
                    self.selected_colors[color_key] = color_value
                    print(f"✅ Set {atom} to {color} ({color_value})")
                else:
                    print("Format: atom_name color_name")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input.")
                break
    
    def _set_sphere_colors(self):
        """Set colors for spheres and dummy atoms."""
        print("\n" + "-"*40)
        print("SPHERE COLOR SELECTION")
        print("-"*40)
        
        print("\nSphere types:")
        print("- default (general spheres)")
        print("- com (center of mass spheres)")
        print("- centroid (geometric centroid spheres)")
        print("- custom (user-defined spheres)")
        
        print("\nAvailable colors:")
        color_list = list(self.color_presets.items())
        for i, (name, hex_code) in enumerate(color_list, 1):
            print(f"{i:2d}. {name:8s} - {hex_code}")
        
        print("\nFormat: sphere_type color_name")
        print("Examples: default red, com blue, custom green")
        print("Type 'done' when finished")
        
        while True:
            try:
                selection = input("Sphere color (or 'done'): ").strip()
                if selection.lower() == 'done':
                    break
                
                if ' ' in selection:
                    sphere_type, color = selection.split(' ', 1)
                    sphere_type = sphere_type.strip().lower()
                    color = color.strip().lower()
                    
                    # Validate sphere type
                    valid_spheres = ['default', 'com', 'centroid', 'custom']
                    if sphere_type not in valid_spheres:
                        print(f"Invalid sphere type: {sphere_type}")
                        print(f"Valid types: {', '.join(valid_spheres)}")
                        continue
                    
                    # Validate color
                    if color not in self.color_presets:
                        print(f"Invalid color: {color}")
                        continue
                    
                    # Set the color
                    color_key = f'sphere_{sphere_type}'
                    color_value = self.color_presets[color]
                    self.selected_colors[color_key] = color_value
                    print(f"✅ Set {sphere_type} sphere to {color} ({color_value})")
                else:
                    print("Format: sphere_type color_name")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input.")
                break
    
    def _set_bond_colors(self):
        """Set colors for bonds and sticks."""
        print("\n" + "-"*40)
        print("BOND COLOR SELECTION")
        print("-"*40)
        
        print("\nBond types:")
        print("- default (general bonds)")
        print("- single (single bonds)")
        print("- double (double bonds)")
        print("- triple (triple bonds)")
        print("- selected (selected bonds)")
        print("- highlight (highlighted bonds)")
        
        print("\nAvailable colors:")
        color_list = list(self.color_presets.items())
        for i, (name, hex_code) in enumerate(color_list, 1):
            print(f"{i:2d}. {name:8s} - {hex_code}")
        
        print("\nFormat: bond_type color_name")
        print("Examples: default red, single blue, double green")
        print("Type 'done' when finished")
        
        while True:
            try:
                selection = input("Bond color (or 'done'): ").strip()
                if selection.lower() == 'done':
                    break
                
                if ' ' in selection:
                    bond_type, color = selection.split(' ', 1)
                    bond_type = bond_type.strip().lower()
                    color = color.strip().lower()
                    
                    # Validate bond type
                    valid_bonds = ['default', 'single', 'double', 'triple', 'selected', 'highlight']
                    if bond_type not in valid_bonds:
                        print(f"Invalid bond type: {bond_type}")
                        print(f"Valid types: {', '.join(valid_bonds)}")
                        continue
                    
                    # Validate color
                    if color not in self.color_presets:
                        print(f"Invalid color: {color}")
                        continue
                    
                    # Set the color
                    color_key = f'stick_{bond_type}'
                    color_value = self.color_presets[color]
                    self.selected_colors[color_key] = color_value
                    print(f"✅ Set {bond_type} bonds to {color} ({color_value})")
                else:
                    print("Format: bond_type color_name")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input.")
                break
    
    def _show_current_colors(self):
        """Display current color settings."""
        print("\n" + "-"*40)
        print("CURRENT COLOR SETTINGS")
        print("-"*40)
        
        if not self.selected_colors:
            print("No colors selected yet.")
            return
        
        print("\nAtom Colors:")
        for key, value in self.selected_colors.items():
            if key.startswith('atom_'):
                atom = key[5:].upper()
                print(f"  {atom}: {value}")
        
        print("\nSphere Colors:")
        for key, value in self.selected_colors.items():
            if key.startswith('sphere_'):
                sphere_type = key[7:]
                print(f"  {sphere_type}: {value}")
        
        print("\nBond Colors:")
        for key, value in self.selected_colors.items():
            if key.startswith('stick_'):
                bond_type = key[6:]
                print(f"  {bond_type}: {value}")
        
        print(f"\nTotal colors selected: {len(self.selected_colors)}")


def show_color_selection_gui() -> Dict[str, str]:
    """Show the color selection GUI."""
    gui = ColorSelectionGUI()
    return gui.show_color_menu()


def apply_selected_colors(colors: Dict[str, str]):
    """Apply selected colors to the theme."""
    from src.shared.ui.theme import COLORS
    COLORS.update(colors)
    print(f"Applied {len(colors)} colors to theme")


def get_color_presets() -> Dict[str, str]:
    """Get available color presets."""
    gui = ColorSelectionGUI()
    return gui.color_presets.copy()
