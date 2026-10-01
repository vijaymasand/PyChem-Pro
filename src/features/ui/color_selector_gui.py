"""
GUI-based color selector for spheres and sticks.

Provides a proper GUI interface for color selection without requiring console input.
"""

from typing import Dict, List, Optional, Tuple
from src.shared.ui.theme import COLORS


class ColorSelectorGUI:
    """
    GUI-based color selector that works without console input.
    """
    
    def __init__(self):
        self.current_palette_index = 0
        self.sphere_colors = {}  # Individual sphere colors
        self.stick_colors = {}    # Individual stick colors
        
        # Predefined color palettes
        self.color_palettes = self._create_color_palettes()
        
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
        if self.current_palette_index < len(self.color_palettes):
            name, description, colors = self.color_palettes[self.current_palette_index]
            return name, description, colors
        return "Unknown", "No palette available", {}
    
    def cycle_palette(self) -> Dict[str, str]:
        """Cycle to the next color palette."""
        self.current_palette_index = (self.current_palette_index + 1) % len(self.color_palettes)
        name, description, colors = self.get_current_palette()
        return colors.copy()
    
    def set_palette_by_index(self, index: int) -> bool:
        """Set a specific palette by index."""
        if 0 <= index < len(self.color_palettes):
            self.current_palette_index = index
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
        return self.sphere_colors.get(sphere_id, self.get_color_by_type('sphere_default'))
    
    def set_stick_color(self, stick_id: str, color: str) -> None:
        """Set color for a specific stick/bond."""
        self.stick_colors[stick_id] = color
    
    def get_stick_color(self, stick_id: str) -> str:
        """Get color for a specific stick/bond."""
        return self.stick_colors.get(stick_id, self.get_color_by_type('stick_default'))
    
    def get_current_colors(self) -> Dict[str, str]:
        """Get all current colors."""
        name, description, colors = self.get_current_palette()
        return colors.copy()
    
    def get_color_by_type(self, color_type: str) -> str:
        """Get color by type (sphere_default, sphere_selected, etc.)."""
        name, description, colors = self.get_current_palette()
        return colors.get(color_type, '#ffffff')
    
    def apply_custom_colors(self, custom_colors: Dict[str, str]) -> None:
        """Apply custom color definitions."""
        # This would modify the current palette
        pass
    
    def reset_to_defaults(self) -> None:
        """Reset to default palette."""
        self.current_palette_index = 0
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
    
    def get_status_message(self) -> str:
        """Get a user-friendly status message."""
        name, description, _ = self.get_current_palette()
        return f"Applied {name}: {description}"


class QuickColorActions:
    """
    Quick color actions for GUI integration.
    """
    
    def __init__(self):
        self.color_selector = ColorSelectorGUI()
        
    def cycle_palette_action(self) -> Dict[str, str]:
        """Action: Cycle to next color palette."""
        colors = self.color_selector.cycle_palette()
        return colors
    
    def set_scientific_palette(self) -> Dict[str, str]:
        """Action: Set Scientific palette."""
        if self.color_selector.set_palette_by_index(0):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_ocean_palette(self) -> Dict[str, str]:
        """Action: Set Ocean palette."""
        if self.color_selector.set_palette_by_index(1):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_forest_palette(self) -> Dict[str, str]:
        """Action: Set Forest palette."""
        if self.color_selector.set_palette_by_index(2):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_sunset_palette(self) -> Dict[str, str]:
        """Action: Set Sunset palette."""
        if self.color_selector.set_palette_by_index(3):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_galaxy_palette(self) -> Dict[str, str]:
        """Action: Set Galaxy palette."""
        if self.color_selector.set_palette_by_index(4):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_candy_palette(self) -> Dict[str, str]:
        """Action: Set Candy palette."""
        if self.color_selector.set_palette_by_index(5):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_monochrome_palette(self) -> Dict[str, str]:
        """Action: Set Monochrome palette."""
        if self.color_selector.set_palette_by_index(6):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_neon_palette(self) -> Dict[str, str]:
        """Action: Set Neon palette."""
        if self.color_selector.set_palette_by_index(7):
            return self.color_selector.get_current_colors()
        return {}
    
    def set_sphere_colors(self, sphere_colors: Dict[str, str]) -> None:
        """Set multiple sphere colors at once."""
        for sphere_id, color in sphere_colors.items():
            self.color_selector.set_sphere_color(sphere_id, color)
    
    def set_stick_colors(self, stick_colors: Dict[str, str]) -> None:
        """Set multiple stick colors at once."""
        for stick_id, color in stick_colors.items():
            self.color_selector.set_stick_color(stick_id, color)
    
    def get_current_status(self) -> str:
        """Get current status message."""
        return self.color_selector.get_status_message()


# Global color selector instance
_color_selector = None

def get_color_selector() -> ColorSelectorGUI:
    """Get the global color selector instance."""
    global _color_selector
    if _color_selector is None:
        _color_selector = ColorSelectorGUI()
    return _color_selector


def get_quick_actions() -> QuickColorActions:
    """Get quick color actions for GUI integration."""
    return QuickColorActions()


def cycle_color_palette() -> Dict[str, str]:
    """Cycle to the next color palette."""
    selector = get_color_selector()
    return selector.cycle_palette()


def get_current_palette_colors() -> Dict[str, str]:
    """Get current palette colors."""
    selector = get_color_selector()
    return selector.get_current_colors()


def set_sphere_color(sphere_id: str, color: str) -> None:
    """Set color for a specific sphere."""
    selector = get_color_selector()
    selector.set_sphere_color(sphere_id, color)


def get_sphere_color(sphere_id: str) -> str:
    """Get color for a specific sphere."""
    selector = get_color_selector()
    return selector.get_sphere_color(sphere_id)


def set_stick_color(stick_id: str, color: str) -> None:
    """Set color for a specific stick/bond."""
    selector = get_color_selector()
    selector.set_stick_color(stick_id, color)


def get_stick_color(stick_id: str) -> str:
    """Get color for a specific stick/bond."""
    selector = get_color_selector()
    return selector.get_stick_color(stick_id)


def get_color_status() -> str:
    """Get current color status message."""
    selector = get_color_selector()
    return selector.get_status_message()
