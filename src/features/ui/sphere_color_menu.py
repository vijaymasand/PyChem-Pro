"""
Interactive color menu for spheres and sticks customization.

Provides a comprehensive color selection interface with preset colors,
custom color picker, and individual sphere coloring capabilities.
"""

from typing import Dict, List, Optional, Tuple
from src.shared.ui.theme import COLORS


class SphereColorMenu:
    """
    Interactive color menu for sphere and stick customization.
    """
    
    def __init__(self):
        self.current_colors = {}
        self.sphere_colors = {}  # Individual sphere colors
        self.stick_colors = {}    # Individual stick colors
        
        # Predefined color palettes
        self.color_palettes = self._create_color_palettes()
        self.current_palette = 0
        
    def _create_color_palettes(self) -> List[Dict[str, str]]:
        """Create predefined color palettes."""
        palettes = []
        
        # 1. Scientific Palette
        scientific = {
            'sphere_default': '#ff6b6b',    # Coral red
            'sphere_selected': '#4ecdc4',  # Cyan
            'sphere_highlight': '#feca57',  # Yellow
            'stick_default': '#95a5a6',    # Gray
            'stick_selected': '#e74c3c',    # Red
            'stick_highlight': '#3498db',  # Blue
        }
        palettes.append(("Scientific", "Professional scientific colors", scientific))
        
        # 2. Ocean Palette
        ocean = {
            'sphere_default': '#006994',    # Deep blue
            'sphere_selected': '#00a8cc',  # Light blue
            'sphere_highlight': '#76d7c4',  # Aqua
            'stick_default': '#1e3a5f',    # Navy
            'stick_selected': '#2e86ab',    # Sky blue
            'stick_highlight': '#48cae4',  # Light cyan
        }
        palettes.append(("Ocean", "Ocean-inspired colors", ocean))
        
        # 3. Forest Palette
        forest = {
            'sphere_default': '#2d6a4f',    # Dark green
            'sphere_selected': '#52b788',  # Medium green
            'sphere_highlight': '#95e1d3',  # Light green
            'stick_default': '#1b5e20',    # Forest green
            'stick_selected': '#4a7c59',    # Olive green
            'stick_highlight': '#a8e6cf',  # Mint green
        }
        palettes.append(("Forest", "Forest and nature colors", forest))
        
        # 4. Sunset Palette
        sunset = {
            'sphere_default': '#ff6b35',    # Orange
            'sphere_selected': '#f7931e',  # Dark orange
            'sphere_highlight': '#fdc830',  # Yellow-orange
            'stick_default': '#c44536',    # Brown
            'stick_selected': '#e74c3c',    # Red
            'stick_highlight': '#f39c12',  # Gold
        }
        palettes.append(("Sunset", "Warm sunset colors", sunset))
        
        # 5. Galaxy Palette
        galaxy = {
            'sphere_default': '#6c5ce7',    # Purple
            'sphere_selected': '#a29bfe',  # Light purple
            'sphere_highlight': '#c471ed',  # Violet
            'stick_default': '#2d3561',    # Dark blue
            'stick_selected': '#4834d4',    # Indigo
            'stick_highlight': '#667eea',  # Periwinkle
        }
        palettes.append(("Galaxy", "Space and galaxy colors", galaxy))
        
        # 6. Candy Palette
        candy = {
            'sphere_default': '#ff6b9d',    # Pink
            'sphere_selected': '#c44569',  # Purple
            'sphere_highlight': '#f8b500',  # Orange
            'stick_default': '#00bcd4',    # Light blue
            'stick_selected': '#ff6b6b',    # Coral
            'stick_highlight': '#4ecdc4',  # Cyan
        }
        palettes.append(("Candy", "Sweet candy colors", candy))
        
        # 7. Monochrome Palette
        monochrome = {
            'sphere_default': '#2c3e50',    # Dark gray
            'sphere_selected': '#34495e',  # Medium gray
            'sphere_highlight': '#95a5a6',  # Light gray
            'stick_default': '#000000',    # Black
            'stick_selected': '#7f8c8d',    # Gray
            'stick_highlight': '#bdc3c7',  # Light gray
        }
        palettes.append(("Monochrome", "Professional grayscale", monochrome))
        
        # 8. Neon Palette
        neon = {
            'sphere_default': '#ff00ff',    # Magenta
            'sphere_selected': '#00ff00',  # Lime green
            'sphere_highlight': '#00ffff',  # Cyan
            'stick_default': '#ff0080',    # Neon pink
            'stick_selected': '#ffff00',    # Yellow
            'stick_highlight': '#ff8000',  # Orange
        }
        palettes.append(("Neon", "High-contrast neon colors", neon))
        
        return palettes
    
    def get_current_palette(self) -> Tuple[str, str, Dict[str, str]]:
        """Get current palette name, description, and colors."""
        if self.current_palette < len(self.color_palettes):
            name, description, colors = self.color_palettes[self.current_palette]
            return name, description, colors
        return "Unknown", "No palette available", {}
    
    def cycle_palette(self) -> Dict[str, str]:
        """Cycle to the next color palette."""
        self.current_palette = (self.current_palette + 1) % len(self.color_palettes)
        name, description, colors = self.get_current_palette()
        self.current_colors.update(colors)
        return colors.copy()
    
    def set_palette_by_name(self, palette_name: str) -> bool:
        """Set a specific palette by name."""
        for i, (name, description, colors) in enumerate(self.color_palettes):
            if name == palette_name:
                self.current_palette = i
                self.current_colors.update(colors)
                return True
        return False
    
    def get_all_palettes(self) -> List[Tuple[str, str]]:
        """Get list of all palette names and descriptions."""
        return [(name, desc) for name, desc, colors in self.color_palettes]
    
    def set_sphere_color(self, sphere_id: str, color: str) -> None:
        """Set color for a specific sphere."""
        self.sphere_colors[sphere_id] = color
    
    def get_sphere_color(self, sphere_id: str) -> str:
        """Get color for a specific sphere."""
        return self.sphere_colors.get(sphere_id, self.current_colors.get('sphere_default', '#ffffff'))
    
    def set_stick_color(self, stick_id: str, color: str) -> None:
        """Set color for a specific stick/bond."""
        self.stick_colors[stick_id] = color
    
    def get_stick_color(self, stick_id: str) -> str:
        """Get color for a specific stick/bond."""
        return self.stick_colors.get(stick_id, self.current_colors.get('stick_default', '#808080'))
    
    def get_current_colors(self) -> Dict[str, str]:
        """Get all current colors."""
        return self.current_colors.copy()
    
    def get_color_by_type(self, color_type: str) -> str:
        """Get color by type (sphere_default, sphere_selected, etc.)."""
        return self.current_colors.get(color_type, '#ffffff')
    
    def apply_custom_colors(self, custom_colors: Dict[str, str]) -> None:
        """Apply custom color definitions."""
        self.current_colors.update(custom_colors)
    
    def reset_to_defaults(self) -> None:
        """Reset to default palette."""
        self.current_palette = 0
        name, description, colors = self.color_palettes[0]
        self.current_colors = colors.copy()
        self.sphere_colors.clear()
        self.stick_colors.clear()
    
    def get_palette_info(self) -> str:
        """Get formatted information about current palette."""
        name, description, colors = self.get_current_palette()
        info = f"Current Palette: {name}\n"
        info += f"Description: {description}\n"
        info += f"Colors: {len(colors)} color definitions\n"
        info += f"Individual spheres: {len(self.sphere_colors)}\n"
        info += f"Individual sticks: {len(self.stick_colors)}"
        return info


class InteractiveColorMenu:
    """
    Interactive color menu with user input capabilities.
    """
    
    def __init__(self):
        self.color_menu = SphereColorMenu()
        
    def show_menu(self) -> None:
        """Display the interactive color menu."""
        print("\n" + "="*70)
        print("SMILES Molecular Toolkit - Interactive Color Menu")
        print("="*70)
        
        # Show current palette
        name, description, colors = self.color_menu.get_current_palette()
        print(f"Current Palette: {name}")
        print(f"Description: {description}")
        print()
        
        # Show available palettes
        palettes = self.color_menu.get_all_palettes()
        print("Available Color Palettes:")
        print("-" * 40)
        for i, (palette_name, palette_desc) in enumerate(palettes):
            marker = "→" if i == self.color_menu.current_palette else " "
            print(f"{marker} {i+1}. {palette_name}")
            print(f"     {palette_desc}")
        
        print()
        print("Color Options:")
        print("-" * 40)
        print("1. Cycle Palette - Switch to next color palette")
        print("2. Set Palette - Choose specific palette by name")
        print("3. Set Sphere Color - Color individual sphere")
        print("4. Set Stick Color - Color individual bonds/sticks")
        print("5. Custom Colors - Define custom colors")
        print("6. Show Current - Display current color settings")
        print("7. Reset Defaults - Reset to default palette")
        print("0. Exit - Return to main application")
        print()
        
    def get_user_choice(self) -> int:
        """Get user menu choice."""
        try:
            choice = input("Enter your choice (0-7): ").strip()
            return int(choice) if choice.isdigit() else -1
        except (ValueError, KeyboardInterrupt):
            return -1
    
    def cycle_palette_option(self) -> None:
        """Handle palette cycling."""
        colors = self.color_menu.cycle_palette()
        name, description, _ = self.color_menu.get_current_palette()
        
        print(f"\n✅ Switched to {name} palette")
        print(f"   Description: {description}")
        print(f"   Colors: {len(colors)} definitions applied")
        
        # Apply to global theme
        from src.shared.ui.theme import COLORS
        COLORS.update(colors)
        
        print("   Colors applied to 3D viewer")
    
    def set_palette_option(self) -> None:
        """Handle specific palette selection."""
        palettes = self.color_menu.get_all_palettes()
        
        print("\nAvailable Palettes:")
        for i, (name, desc) in enumerate(palettes):
            print(f"   {i+1}. {name} - {desc}")
        
        try:
            choice = input("Enter palette number: ").strip()
            if choice.isdigit():
                palette_idx = int(choice) - 1
                if 0 <= palette_idx < len(palettes):
                    palette_name = palettes[palette_idx][0]
                    if self.color_menu.set_palette_by_name(palette_name):
                        colors = self.color_menu.get_current_colors()
                        
                        print(f"\n✅ Applied {palette_name} palette")
                        print(f"   Colors: {len(colors)} definitions applied")
                        
                        # Apply to global theme
                        from src.shared.ui.theme import COLORS
                        COLORS.update(colors)
                        
                        print("   Colors applied to 3D viewer")
                    else:
                        print(f"\n❌ Failed to apply {palette_name} palette")
                else:
                    print(f"\n❌ Invalid palette number: {choice}")
            else:
                print(f"\n❌ Please enter a valid number")
        except (ValueError, KeyboardInterrupt):
            print("\n❌ Invalid input")
    
    def set_sphere_color_option(self) -> None:
        """Handle individual sphere coloring."""
        try:
            sphere_id = input("Enter sphere ID (e.g., 'sphere_1', 'COM', 'centroid'): ").strip()
            if not sphere_id:
                print("\n❌ Sphere ID cannot be empty")
                return
            
            color = input("Enter color (hex format, e.g., '#ff0000' or 'red'): ").strip()
            if not color:
                print("\n❌ Color cannot be empty")
                return
            
            # Validate color format
            if not color.startswith('#') and color.lower() not in ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 'pink', 'cyan', 'magenta', 'white', 'black', 'gray']:
                print(f"\n❌ Invalid color format: {color}")
                print("   Use hex format (#ff0000) or color name (red)")
                return
            
            self.color_menu.set_sphere_color(sphere_id, color)
            print(f"\n✅ Set {sphere_id} color to {color}")
            
        except (ValueError, KeyboardInterrupt):
            print("\n❌ Invalid input")
    
    def set_stick_color_option(self) -> None:
        """Handle individual stick/bond coloring."""
        try:
            stick_id = input("Enter stick/bond ID (e.g., 'bond_1', 'stick_1'): ").strip()
            if not stick_id:
                print("\n❌ Stick ID cannot be empty")
                return
            
            color = input("Enter color (hex format, e.g., '#00ff00' or 'green'): ").strip()
            if not color:
                print("\n❌ Color cannot be empty")
                return
            
            # Validate color format
            if not color.startswith('#') and color.lower() not in ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 'pink', 'cyan', 'magenta', 'white', 'black', 'gray']:
                print(f"\n❌ Invalid color format: {color}")
                print("   Use hex format (#00ff00) or color name (green)")
                return
            
            self.color_menu.set_stick_color(stick_id, color)
            print(f"\n✅ Set {stick_id} color to {color}")
            
        except (ValueError, KeyboardInterrupt):
            print("\n❌ Invalid input")
    
    def custom_colors_option(self) -> None:
        """Handle custom color definition."""
        print("\nCustom Color Definition")
        print("-" * 30)
        print("Enter custom colors for different elements:")
        print("Format: element_name=color (e.g., carbon=#ff0000)")
        print("Type 'done' when finished")
        
        custom_colors = {}
        
        while True:
            try:
                color_def = input("Color definition (or 'done'): ").strip()
                if color_def.lower() == 'done':
                    break
                
                if '=' in color_def:
                    element, color = color_def.split('=', 1)
                    element = element.strip()
                    color = color.strip()
                    
                    # Validate color
                    if color.startswith('#') or color.lower() in ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 'pink', 'cyan', 'magenta', 'white', 'black', 'gray']:
                        custom_colors[element] = color
                        print(f"   ✅ Set {element} = {color}")
                    else:
                        print(f"   ❌ Invalid color: {color}")
                else:
                    print(f"   ❌ Invalid format. Use: element=color")
                    
            except (ValueError, KeyboardInterrupt):
                print("\n❌ Invalid input")
                break
        
        if custom_colors:
            self.color_menu.apply_custom_colors(custom_colors)
            print(f"\n✅ Applied {len(custom_colors)} custom colors")
            
            # Apply to global theme
            from src.shared.ui.theme import COLORS
            COLORS.update(custom_colors)
            
            print("   Colors applied to 3D viewer")
    
    def show_current_option(self) -> None:
        """Display current color settings."""
        info = self.color_menu.get_palette_info()
        print(f"\n{info}")
        
        if self.color_menu.sphere_colors:
            print("\nIndividual Sphere Colors:")
            for sphere_id, color in self.color_menu.sphere_colors.items():
                print(f"   {sphere_id}: {color}")
        
        if self.color_menu.stick_colors:
            print("\nIndividual Stick Colors:")
            for stick_id, color in self.color_menu.stick_colors.items():
                print(f"   {stick_id}: {color}")
    
    def reset_defaults_option(self) -> None:
        """Reset to default colors."""
        self.color_menu.reset_to_defaults()
        
        # Apply to global theme
        from src.shared.ui.theme import COLORS
        COLORS.update(self.color_menu.get_current_colors())
        
        name, description, _ = self.color_menu.get_current_palette()
        print(f"\n✅ Reset to {name} palette")
        print(f"   Description: {description}")
        print("   Colors applied to 3D viewer")
    
    def run_interactive_menu(self) -> None:
        """Run the interactive color menu."""
        while True:
            self.show_menu()
            choice = self.get_user_choice()
            
            if choice == 0:
                print("\n👋 Exiting color menu...")
                break
            elif choice == 1:
                self.cycle_palette_option()
            elif choice == 2:
                self.set_palette_option()
            elif choice == 3:
                self.set_sphere_color_option()
            elif choice == 4:
                self.set_stick_color_option()
            elif choice == 5:
                self.custom_colors_option()
            elif choice == 6:
                self.show_current_option()
            elif choice == 7:
                self.reset_defaults_option()
            else:
                print("\n❌ Invalid choice. Please enter 0-7.")
            
            input("\nPress Enter to continue...")


# Global color menu instance
_color_menu = None

def get_color_menu() -> InteractiveColorMenu:
    """Get the global color menu instance."""
    global _color_menu
    if _color_menu is None:
        _color_menu = InteractiveColorMenu()
    return _color_menu


def show_interactive_color_menu():
    """Show the interactive color menu."""
    menu = get_color_menu()
    menu.run_interactive_menu()


def get_current_palette_colors() -> Dict[str, str]:
    """Get current palette colors."""
    menu = get_color_menu()
    return menu.color_menu.get_current_colors()


def set_sphere_color(sphere_id: str, color: str) -> None:
    """Set color for a specific sphere."""
    menu = get_color_menu()
    menu.color_menu.set_sphere_color(sphere_id, color)


def set_stick_color(stick_id: str, color: str) -> None:
    """Set color for a specific stick/bond."""
    menu = get_color_menu()
    menu.color_menu.set_stick_color(stick_id, color)


def cycle_color_palette() -> Dict[str, str]:
    """Cycle to the next color palette."""
    menu = get_color_menu()
    return menu.color_menu.cycle_palette()
