"""
Quick color picker that works without console input.

Simple dropdown-style color selection for GUI applications.
"""

from typing import Dict, List, Optional


class QuickColorPicker:
    """
    Quick color picker with predefined color schemes.
    """
    
    def __init__(self):
        self.current_scheme = 0
        self.color_schemes = self._create_color_schemes()
        
    def _create_color_schemes(self) -> List[Dict[str, str]]:
        """Create predefined color schemes."""
        schemes = []
        
        # Scheme 1: Default
        default = {
            'atom_c': '#55ff7f', 'atom_h': '#d0d0d0', 'atom_o': '#ff0d0d', 'atom_n': '#3050f8',
            'sphere_default': '#ff00ff', 'sphere_com': '#ff00ff', 'sphere_centroid': '#00ff00', 'sphere_custom': '#ffff00',
            'stick_default': '#808080'
        }
        schemes.append(("Default", default))
        
        # Scheme 2: Red Theme
        red_theme = {
            'atom_c': '#ff0000', 'atom_h': '#ffffff', 'atom_o': '#cc0000', 'atom_n': '#990000',
            'sphere_default': '#ff6666', 'sphere_com': '#ff3333', 'sphere_centroid': '#ff0000', 'sphere_custom': '#cc0000',
            'stick_default': '#ff9999'
        }
        schemes.append(("Red Theme", red_theme))
        
        # Scheme 3: Blue Theme
        blue_theme = {
            'atom_c': '#0066cc', 'atom_h': '#ffffff', 'atom_o': '#004499', 'atom_n': '#003366',
            'sphere_default': '#6699ff', 'sphere_com': '#3366ff', 'sphere_centroid': '#0066ff', 'sphere_custom': '#0033ff',
            'stick_default': '#99ccff'
        }
        schemes.append(("Blue Theme", blue_theme))
        
        # Scheme 4: Green Theme
        green_theme = {
            'atom_c': '#00aa00', 'atom_h': '#ffffff', 'atom_o': '#008800', 'atom_n': '#006600',
            'sphere_default': '#66ff66', 'sphere_com': '#33ff33', 'sphere_centroid': '#00ff00', 'sphere_custom': '#00cc00',
            'stick_default': '#99ff99'
        }
        schemes.append(("Green Theme", green_theme))
        
        # Scheme 5: Purple Theme
        purple_theme = {
            'atom_c': '#9933cc', 'atom_h': '#ffffff', 'atom_o': '#772299', 'atom_n': '#551177',
            'sphere_default': '#cc99ff', 'sphere_com': '#ff66ff', 'sphere_centroid': '#cc33ff', 'sphere_custom': '#9933ff',
            'stick_default': '#ddbbff'
        }
        schemes.append(("Purple Theme", purple_theme))
        
        # Scheme 6: Orange Theme
        orange_theme = {
            'atom_c': '#ff8800', 'atom_h': '#ffffff', 'atom_o': '#cc6600', 'atom_n': '#994400',
            'sphere_default': '#ffcc66', 'sphere_com': '#ff9933', 'sphere_centroid': '#ff6600', 'sphere_custom': '#cc4400',
            'stick_default': '#ffdd99'
        }
        schemes.append(("Orange Theme", orange_theme))
        
        return schemes
    
    def get_current_scheme(self) -> str:
        """Get current color scheme name."""
        if self.current_scheme < len(self.color_schemes):
            return self.color_schemes[self.current_scheme][0]
        return "Unknown"
    
    def cycle_scheme(self) -> Dict[str, str]:
        """Cycle to next color scheme."""
        self.current_scheme = (self.current_scheme + 1) % len(self.color_schemes)
        name, colors = self.color_schemes[self.current_scheme]
        return colors.copy()
    
    def get_scheme_by_name(self, scheme_name: str) -> Optional[Dict[str, str]]:
        """Get a specific color scheme by name."""
        for name, colors in self.color_schemes:
            if name == scheme_name:
                return colors.copy()
        return None
    
    def get_available_schemes(self) -> List[str]:
        """Get list of available scheme names."""
        return [name for name, colors in self.color_schemes]
    
    def apply_scheme(self, colors: Dict[str, str]) -> None:
        """Apply colors to the theme."""
        from src.shared.ui.theme import COLORS
        COLORS.update(colors)
    
    def get_status_message(self) -> str:
        """Get status message for current scheme."""
        scheme_name = self.get_current_scheme()
        return f"Applied {scheme_name}"


# Global instance
_color_picker = None

def get_color_picker() -> QuickColorPicker:
    """Get the global color picker instance."""
    global _color_picker
    if _color_picker is None:
        _color_picker = QuickColorPicker()
    return _color_picker


def cycle_color_scheme() -> Dict[str, str]:
    """Cycle to next color scheme."""
    picker = get_color_picker()
    colors = picker.cycle_scheme()
    picker.apply_scheme(colors)
    return colors


def get_color_status() -> str:
    """Get current color status."""
    picker = get_color_picker()
    return picker.get_status_message()


def get_available_schemes() -> List[str]:
    """Get available color schemes."""
    picker = get_color_picker()
    return picker.get_available_schemes()


def apply_specific_scheme(scheme_name: str) -> bool:
    """Apply a specific color scheme by name."""
    picker = get_color_picker()
    colors = picker.get_scheme_by_name(scheme_name)
    if colors:
        picker.apply_scheme(colors)
        return True
    return False
