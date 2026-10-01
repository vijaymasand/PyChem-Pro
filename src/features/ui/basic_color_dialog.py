"""
Basic color dialog that works without complex dependencies.

Simple, direct solution for color selection.
"""

from typing import Dict, Optional


class BasicColorDialog:
    """
    Basic color dialog for simple color selection.
    """
    
    def __init__(self):
        self.selected_colors = {}
        
        # Basic color options
        self.basic_colors = {
            'Red': '#ff0000',
            'Green': '#00ff00', 
            'Blue': '#0000ff',
            'Yellow': '#ffff00',
            'Orange': '#ff8000',
            'Purple': '#8000ff',
            'Pink': '#ff00ff',
            'Cyan': '#00ffff',
            'White': '#ffffff',
            'Black': '#000000',
            'Gray': '#808080',
            'Brown': '#8b4513'
        }
        
    def get_color_selection(self) -> Dict[str, str]:
        """Get color selection from user."""
        print("\n" + "="*50)
        print("Basic Color Selection")
        print("="*50)
        
        print("\nAvailable colors:")
        for i, (name, hex_code) in enumerate(self.basic_colors.items(), 1):
            print(f"{i:2d}. {name:8s} - {hex_code}")
        
        print("\nSelect color to apply:")
        print("1. Atoms (C, H, O, N, etc.)")
        print("2. Spheres (COM, Centroid, Custom)")
        print("3. Bonds/Sticks")
        print("0. Cancel")
        
        try:
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == '0':
                return {}
            elif choice == '1':
                return self._select_atom_colors()
            elif choice == '2':
                return self._select_sphere_colors()
            elif choice == '3':
                return self._select_bond_colors()
            else:
                print("Invalid choice.")
                return {}
                
        except (ValueError, KeyboardInterrupt):
            return {}
    
    def _select_atom_colors(self) -> Dict[str, str]:
        """Select colors for atoms."""
        print("\nSelect atom color:")
        print("Enter: atom_number color_number")
        print("Example: 1 1 (Carbon = Red)")
        print("Enter 'done' when finished")
        
        colors = {}
        
        while True:
            try:
                selection = input("Atom selection (or 'done'): ").strip()
                if selection.lower() == 'done':
                    break
                
                if ' ' in selection:
                    atom_num, color_num = selection.split(' ', 1)
                    atom_num = int(atom_num)
                    color_num = int(color_num)
                    
                    atoms = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']
                    color_names = list(self.basic_colors.keys())
                    
                    if 1 <= atom_num <= len(atoms) and 1 <= color_num <= len(color_names):
                        atom = atoms[atom_num - 1]
                        color_name = color_names[color_num - 1]
                        color_hex = self.basic_colors[color_name]
                        
                        colors[f'atom_{atom.lower()}'] = color_hex
                        print(f"✅ Set {atom} to {color_name}")
                    else:
                        print("❌ Invalid atom or color number")
                else:
                    print("❌ Format: atom_number color_number")
                    
            except (ValueError, KeyboardInterrupt):
                break
        
        return colors
    
    def _select_sphere_colors(self) -> Dict[str, str]:
        """Select colors for spheres."""
        print("\nSelect sphere color:")
        print("Enter: sphere_type color_number")
        print("Sphere types: default, com, centroid, custom")
        print("Example: com 1 (COM = Red)")
        print("Enter 'done' when finished")
        
        colors = {}
        
        while True:
            try:
                selection = input("Sphere selection (or 'done'): ").strip()
                if selection.lower() == 'done':
                    break
                
                if ' ' in selection:
                    sphere_type, color_num = selection.split(' ', 1)
                    color_num = int(color_num)
                    
                    color_names = list(self.basic_colors.keys())
                    
                    if 1 <= color_num <= len(color_names):
                        color_name = color_names[color_num - 1]
                        color_hex = self.basic_colors[color_name]
                        
                        colors[f'sphere_{sphere_type.lower()}'] = color_hex
                        print(f"✅ Set {sphere_type} sphere to {color_name}")
                    else:
                        print("❌ Invalid color number")
                else:
                    print("❌ Format: sphere_type color_number")
                    
            except (ValueError, KeyboardInterrupt):
                break
        
        return colors
    
    def _select_bond_colors(self) -> Dict[str, str]:
        """Select colors for bonds/sticks."""
        print("\nSelect bond/stick color:")
        print("Enter: bond_type color_number")
        print("Bond types: default, single, double, triple")
        print("Example: default 1 (Default bonds = Red)")
        print("Enter 'done' when finished")
        
        colors = {}
        
        while True:
            try:
                selection = input("Bond selection (or 'done'): ").strip()
                if selection.lower() == 'done':
                    break
                
                if ' ' in selection:
                    bond_type, color_num = selection.split(' ', 1)
                    color_num = int(color_num)
                    
                    color_names = list(self.basic_colors.keys())
                    
                    if 1 <= color_num <= len(color_names):
                        color_name = color_names[color_num - 1]
                        color_hex = self.basic_colors[color_name]
                        
                        colors[f'stick_{bond_type.lower()}'] = color_hex
                        print(f"✅ Set {bond_type} bonds to {color_name}")
                    else:
                        print("❌ Invalid color number")
                else:
                    print("❌ Format: bond_type color_number")
                    
            except (ValueError, KeyboardInterrupt):
                break
        
        return colors


def show_basic_color_dialog() -> Dict[str, str]:
    """Show basic color dialog and return selected colors."""
    dialog = BasicColorDialog()
    return dialog.get_color_selection()


def apply_basic_colors(colors: Dict[str, str]):
    """Apply basic colors to theme."""
    from src.shared.ui.theme import COLORS
    COLORS.update(colors)
    print(f"Applied {len(colors)} colors")
