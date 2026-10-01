"""
Simple GUI color dialog that works without PySide6 dependencies.

Provides a basic color selection interface using console input as fallback.
"""

from typing import Dict, List, Optional, Tuple
from src.shared.ui.theme import COLORS


class SimpleColorDialog:
    """
    Simple color dialog that provides color selection options.
    """
    
    def __init__(self):
        self.selected_colors = {}
        self.color_presets = self._create_color_presets()
        
    def _create_color_presets(self) -> Dict[str, str]:
        """Create predefined color presets."""
        return {
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
            'Brown': '#8b4513',
            'Navy': '#000080',
            'Teal': '#008080',
            'Maroon': '#800000',
            'Lime': '#00ff00',
            'Aqua': '#00ffff',
            'Fuchsia': '#ff00ff',
            'Silver': '#c0c0c0',
            'Olive': '#808000',
        }
    
    def show_color_selection_menu(self) -> Dict[str, str]:
        """Show color selection menu and return selected colors."""
        print("\n" + "="*60)
        print("SMILES Molecular Toolkit - Color Selection")
        print("="*60)
        
        print("\nColor Selection Options:")
        print("-" * 40)
        print("1. Atom Colors - Set colors for different atom types")
        print("2. Sphere Colors - Set colors for different sphere types")
        print("3. Stick/Bond Colors - Set colors for bonds/sticks")
        print("4. Use Preset Colors - Choose from predefined colors")
        print("5. Custom Hex Colors - Enter custom hex color codes")
        print("6. Show Current Colors - Display current color settings")
        print("7. Apply to 3D Viewer - Apply colors and exit")
        print("0. Exit - Cancel")
        
        while True:
            try:
                choice = input("\nEnter your choice (0-7): ").strip()
                
                if choice == '0':
                    print("Color selection cancelled.")
                    return {}
                elif choice == '1':
                    self._set_atom_colors()
                elif choice == '2':
                    self._set_sphere_colors()
                elif choice == '3':
                    self._set_stick_colors()
                elif choice == '4':
                    self._use_preset_colors()
                elif choice == '5':
                    self._set_custom_hex_colors()
                elif choice == '6':
                    self._show_current_colors()
                elif choice == '7':
                    print("Applying colors to 3D viewer...")
                    return self.selected_colors
                else:
                    print("Invalid choice. Please enter 0-7.")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input. Please try again.")
                break
        
        return {}
    
    def _set_atom_colors(self):
        """Set colors for different atom types."""
        print("\n" + "-"*40)
        print("Atom Color Settings")
        print("-"*40)
        
        atom_types = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']
        
        print("Available atom types: " + ", ".join(atom_types))
        print("Enter atom type and color (e.g., 'C red', 'O #00ff00')")
        print("Type 'done' when finished")
        
        while True:
            try:
                color_input = input("Atom color (or 'done'): ").strip()
                if color_input.lower() == 'done':
                    break
                
                if ' ' in color_input:
                    atom_type, color = color_input.split(' ', 1)
                    atom_type = atom_type.strip().upper()
                    color = color.strip()
                    
                    if atom_type in atom_types:
                        # Validate color
                        if self._validate_color(color):
                            color_key = f'atom_{atom_type.lower()}'
                            self.selected_colors[color_key] = color
                            print(f"✅ Set {atom_type} color to {color}")
                        else:
                            print(f"❌ Invalid color: {color}")
                    else:
                        print(f"❌ Unknown atom type: {atom_type}")
                else:
                    print("❌ Format: atom_type color (e.g., 'C red')")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input.")
                break
    
    def _set_sphere_colors(self):
        """Set colors for different sphere types."""
        print("\n" + "-"*40)
        print("Sphere Color Settings")
        print("-"*40)
        
        sphere_types = ['default', 'selected', 'highlight', 'com', 'centroid', 'custom']
        
        print("Available sphere types: " + ", ".join(sphere_types))
        print("Enter sphere type and color (e.g., 'default red', 'com #00ff00')")
        print("Type 'done' when finished")
        
        while True:
            try:
                color_input = input("Sphere color (or 'done'): ").strip()
                if color_input.lower() == 'done':
                    break
                
                if ' ' in color_input:
                    sphere_type, color = color_input.split(' ', 1)
                    sphere_type = sphere_type.strip().lower()
                    color = color.strip()
                    
                    if sphere_type in sphere_types:
                        # Validate color
                        if self._validate_color(color):
                            color_key = f'sphere_{sphere_type}'
                            self.selected_colors[color_key] = color
                            print(f"✅ Set {sphere_type} sphere color to {color}")
                        else:
                            print(f"❌ Invalid color: {color}")
                    else:
                        print(f"❌ Unknown sphere type: {sphere_type}")
                else:
                    print("❌ Format: sphere_type color (e.g., 'default red')")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input.")
                break
    
    def _set_stick_colors(self):
        """Set colors for different stick/bond types."""
        print("\n" + "-"*40)
        print("Stick/Bond Color Settings")
        print("-"*40)
        
        stick_types = ['default', 'selected', 'highlight', 'single', 'double', 'triple']
        
        print("Available stick types: " + ", ".join(stick_types))
        print("Enter stick type and color (e.g., 'default red', 'single #00ff00')")
        print("Type 'done' when finished")
        
        while True:
            try:
                color_input = input("Stick color (or 'done'): ").strip()
                if color_input.lower() == 'done':
                    break
                
                if ' ' in color_input:
                    stick_type, color = color_input.split(' ', 1)
                    stick_type = stick_type.strip().lower()
                    color = color.strip()
                    
                    if stick_type in stick_types:
                        # Validate color
                        if self._validate_color(color):
                            color_key = f'stick_{stick_type}'
                            self.selected_colors[color_key] = color
                            print(f"✅ Set {stick_type} stick color to {color}")
                        else:
                            print(f"❌ Invalid color: {color}")
                    else:
                        print(f"❌ Unknown stick type: {stick_type}")
                else:
                    print("❌ Format: stick_type color (e.g., 'default red')")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input.")
                break
    
    def _use_preset_colors(self):
        """Use predefined color presets."""
        print("\n" + "-"*40)
        print("Preset Color Selection")
        print("-"*40)
        
        print("Available preset colors:")
        preset_list = list(self.color_presets.items())
        for i, (name, hex_code) in enumerate(preset_list):
            print(f"{i+1:2d}. {name:10s} - {hex_code}")
        
        print("\nSelect preset color number to apply to:")
        print("1. All atoms (same color)")
        print("2. All spheres (same color)")
        print("3. All sticks (same color)")
        print("4. Specific item")
        
        try:
            preset_choice = input("Preset color number (1-20): ").strip()
            if preset_choice.isdigit():
                preset_idx = int(preset_choice) - 1
                if 0 <= preset_idx < len(preset_list):
                    preset_name, preset_color = preset_list[preset_idx]
                    
                    apply_to = input("Apply to (atoms/spheres/sticks/specific): ").strip().lower()
                    
                    if apply_to == 'atoms':
                        # Apply to all atom types
                        atom_types = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']
                        for atom in atom_types:
                            self.selected_colors[f'atom_{atom.lower()}'] = preset_color
                        print(f"✅ Applied {preset_name} to all atoms")
                    elif apply_to == 'spheres':
                        # Apply to all sphere types
                        sphere_types = ['default', 'selected', 'highlight']
                        for sphere in sphere_types:
                            self.selected_colors[f'sphere_{sphere}'] = preset_color
                        print(f"✅ Applied {preset_name} to all spheres")
                    elif apply_to == 'sticks':
                        # Apply to all stick types
                        stick_types = ['default', 'selected', 'highlight']
                        for stick in stick_types:
                            self.selected_colors[f'stick_{stick}'] = preset_color
                        print(f"✅ Applied {preset_name} to all sticks")
                    elif apply_to == 'specific':
                        # Apply to specific item
                        item_type = input("Item type (atom/sphere/stick): ").strip().lower()
                        item_name = input("Item name: ").strip()
                        
                        if item_type == 'atom':
                            color_key = f'atom_{item_name.lower()}'
                        elif item_type == 'sphere':
                            color_key = f'sphere_{item_name.lower()}'
                        elif item_type == 'stick':
                            color_key = f'stick_{item_name.lower()}'
                        else:
                            print("❌ Invalid item type")
                            return
                        
                        self.selected_colors[color_key] = preset_color
                        print(f"✅ Applied {preset_name} to {item_type} {item_name}")
                    else:
                        print("❌ Invalid apply option")
                else:
                    print("❌ Invalid preset number")
            else:
                print("❌ Please enter a valid number")
                
        except (ValueError, KeyboardInterrupt):
            print("\nInvalid input.")
    
    def _set_custom_hex_colors(self):
        """Set custom hex color codes."""
        print("\n" + "-"*40)
        print("Custom Hex Color Settings")
        print("-"*40)
        
        print("Enter color name and hex code (e.g., 'my_red #ff0000')")
        print("Hex format: #RRGGBB (e.g., #ff0000, #00ff00, #0000ff)")
        print("Type 'done' when finished")
        
        while True:
            try:
                color_input = input("Custom color (or 'done'): ").strip()
                if color_input.lower() == 'done':
                    break
                
                if ' ' in color_input:
                    color_name, hex_code = color_input.split(' ', 1)
                    color_name = color_name.strip()
                    hex_code = hex_code.strip()
                    
                    # Validate hex code
                    if self._validate_hex_color(hex_code):
                        self.selected_colors[color_name] = hex_code
                        print(f"✅ Set {color_name} to {hex_code}")
                    else:
                        print(f"❌ Invalid hex color: {hex_code}")
                        print("   Format: #RRGGBB (e.g., #ff0000)")
                else:
                    print("❌ Format: color_name hex_code (e.g., 'my_red #ff0000')")
                    
            except (ValueError, KeyboardInterrupt):
                print("\nInvalid input.")
                break
    
    def _show_current_colors(self):
        """Display current color settings."""
        print("\n" + "-"*40)
        print("Current Color Settings")
        print("-"*40)
        
        if not self.selected_colors:
            print("No colors selected yet.")
            return
        
        print("Selected colors:")
        for color_key, color_value in self.selected_colors.items():
            print(f"  {color_key}: {color_value}")
        
        print(f"\nTotal colors selected: {len(self.selected_colors)}")
    
    def _validate_color(self, color: str) -> bool:
        """Validate color format."""
        # Check if it's a hex color
        if color.startswith('#'):
            return self._validate_hex_color(color)
        
        # Check if it's a predefined color name
        return color.lower() in [name.lower() for name in self.color_presets.keys()]
    
    def _validate_hex_color(self, hex_code: str) -> bool:
        """Validate hex color code."""
        if not hex_code.startswith('#'):
            return False
        
        hex_part = hex_code[1:]
        if len(hex_part) != 6:
            return False
        
        # Check if all characters are valid hex digits
        try:
            int(hex_part, 16)
            return True
        except ValueError:
            return False


def show_color_dialog() -> Dict[str, str]:
    """Show the color dialog and return selected colors."""
    dialog = SimpleColorDialog()
    return dialog.show_color_selection_menu()


def apply_colors_to_theme(colors: Dict[str, str]):
    """Apply selected colors to the global theme."""
    from src.shared.ui.theme import COLORS
    
    # Update global colors
    COLORS.update(colors)
    
    print(f"Applied {len(colors)} colors to theme")


def get_current_theme_colors() -> Dict[str, str]:
    """Get current theme colors."""
    from src.shared.ui.theme import COLORS
    return {k: v for k, v in COLORS.items() if k.startswith('atom_') or k.startswith('sphere_') or k.startswith('stick_')}
