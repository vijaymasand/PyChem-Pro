"""
Predefined color schemes for molecular visualization.

Provides various color schemes for different visualization needs.
"""

from typing import Dict, List, Tuple
from src.shared.ui.theme import COLORS


class ColorScheme:
    """Represents a color scheme with name and colors."""
    
    def __init__(self, name: str, description: str, colors: Dict[str, str]):
        self.name = name
        self.description = description
        self.colors = colors
    
    def apply(self):
        """Apply this color scheme to the global COLORS."""
        COLORS.update(self.colors)


class ColorSchemeManager:
    """Manages predefined and custom color schemes."""
    
    def __init__(self):
        self.schemes = self._create_predefined_schemes()
        self.current_scheme = "Default"
    
    def _create_predefined_schemes(self) -> List[ColorScheme]:
        """Create predefined color schemes."""
        schemes = []
        
        # 1. Default Scheme
        default_colors = {
            'atom_h': '#d0d0d0',      # White
            'atom_c': '#55ff7f',      # Gray
            'atom_n': '#3050f8',      # Blue
            'atom_o': '#ff0d0d',      # Red
            'atom_f': '#90e050',      # Light green
            'atom_p': '#ff8000',      # Orange
            'atom_s': '#ffff30',      # Yellow
            'atom_cl': '#1ff01f',     # Green
            'atom_br': '#a62929',     # Dark red
            'atom_i': '#940094',      # Purple
            'atom_selected': '#ff00ff',    # Magenta
            'atom_highlight': '#ffff00',   # Yellow
            'atom_positive': '#0000ff',    # Blue
            'atom_negative': '#ff0000',    # Red
        }
        schemes.append(ColorScheme("Default", "Standard CPK colors", default_colors))
        
        # 2. Pastel Scheme
        pastel_colors = {
            'atom_h': '#f0f0f0',      # Light gray
            'atom_c': '#ff6b6b',      # Soft red
            'atom_n': '#4ecdc4',      # Soft cyan
            'atom_o': '#45b7d1',      # Soft blue
            'atom_f': '#96ceb4',      # Soft green
            'atom_p': '#feca57',      # Soft yellow
            'atom_s': '#ff9ff3',      # Soft pink
            'atom_cl': '#54a0ff',     # Soft blue
            'atom_br': '#48dbfb',     # Light blue
            'atom_i': '#a29bfe',      # Light purple
            'atom_selected': '#ee5a6f',    # Soft red
            'atom_highlight': '#feca57',   # Soft yellow
            'atom_positive': '#54a0ff',    # Soft blue
            'atom_negative': '#ff6b6b',    # Soft red
        }
        schemes.append(ColorScheme("Pastel", "Soft pastel colors", pastel_colors))
        
        # 3. Vibrant Scheme
        vibrant_colors = {
            'atom_h': '#ffffff',      # White
            'atom_c': '#2ecc71',      # Green
            'atom_n': '#3498db',      # Blue
            'atom_o': '#e74c3c',      # Red
            'atom_f': '#1abc9c',      # Teal
            'atom_p': '#f39c12',      # Orange
            'atom_s': '#f1c40f',      # Yellow
            'atom_cl': '#16a085',     # Dark teal
            'atom_br': '#27ae60',     # Dark green
            'atom_i': '#8e44ad',      # Purple
            'atom_selected': '#e74c3c',    # Red
            'atom_highlight': '#f1c40f',   # Yellow
            'atom_positive': '#3498db',    # Blue
            'atom_negative': '#e74c3c',    # Red
        }
        schemes.append(ColorScheme("Vibrant", "High contrast vibrant colors", vibrant_colors))
        
        # 4. Earth Tones Scheme
        earth_colors = {
            'atom_h': '#d4d4d4',      # Light gray
            'atom_c': '#8b7355',      # Brown
            'atom_n': '#5f8a8b',      # Teal
            'atom_o': '#a0522d',      # Sienna
            'atom_f': '#778899',      # Light slate gray
            'atom_p': '#cd853f',      # Peru
            'atom_s': '#daa520',      # Goldenrod
            'atom_cl': '#6b8e23',     # Olive drab
            'atom_br': '#8b4513',     # Saddle brown
            'atom_i': '#708090',      # Slate gray
            'atom_selected': '#cd5c5c',    # Indian red
            'atom_highlight': '#ffd700',   # Gold
            'atom_positive': '#5f8a8b',    # Teal
            'atom_negative': '#a0522d',    # Sienna
        }
        schemes.append(ColorScheme("Earth Tones", "Natural earth colors", earth_colors))
        
        # 5. High Contrast Scheme
        high_contrast_colors = {
            'atom_h': '#ffffff',      # White
            'atom_c': '#000000',      # Black
            'atom_n': '#0000ff',      # Pure blue
            'atom_o': '#ff0000',      # Pure red
            'atom_f': '#00ff00',      # Pure green
            'atom_p': '#ff00ff',      # Magenta
            'atom_s': '#ffff00',      # Yellow
            'atom_cl': '#00ffff',      # Cyan
            'atom_br': '#ff8800',      # Orange
            'atom_i': '#8800ff',      # Purple
            'atom_selected': '#ff00ff',    # Magenta
            'atom_highlight': '#ffff00',   # Yellow
            'atom_positive': '#0000ff',    # Blue
            'atom_negative': '#ff0000',    # Red
        }
        schemes.append(ColorScheme("High Contrast", "Maximum contrast colors", high_contrast_colors))
        
        # 6. Monochrome Scheme
        monochrome_colors = {
            'atom_h': '#ffffff',      # White
            'atom_c': '#808080',      # Medium gray
            'atom_n': '#606060',      # Dark gray
            'atom_o': '#404040',      # Darker gray
            'atom_f': '#a0a0a0',      # Light gray
            'atom_p': '#909090',      # Medium light gray
            'atom_s': '#707070',      # Medium dark gray
            'atom_cl': '#505050',      # Dark gray
            'atom_br': '#303030',      # Very dark gray
            'atom_i': '#202020',      # Very very dark gray
            'atom_selected': '#000000',    # Black
            'atom_highlight': '#ffffff',   # White
            'atom_positive': '#303030',    # Dark gray
            'atom_negative': '#606060',    # Medium gray
        }
        schemes.append(ColorScheme("Monochrome", "Grayscale colors", monochrome_colors))
        
        # 7. Protein Scheme (for biomolecules)
        protein_colors = {
            'atom_h': '#ffffff',      # White
            'atom_c': '#808080',      # Gray
            'atom_n': '#0000ff',      # Blue (basic residues)
            'atom_o': '#ff0000',      # Red (acidic residues)
            'atom_f': '#00ff00',      # Green (hydrophobic)
            'atom_p': '#ff8800',      # Orange (phosphate)
            'atom_s': '#ffff00',      # Yellow (sulfur)
            'atom_cl': '#00ffff',      # Cyan
            'atom_br': '#ff00ff',      # Magenta
            'atom_i': '#8800ff',      # Purple
            'atom_selected': '#ffff00',    # Yellow
            'atom_highlight': '#00ff00',   # Green
            'atom_positive': '#0000ff',    # Blue (positive)
            'atom_negative': '#ff0000',    # Red (negative)
        }
        schemes.append(ColorScheme("Protein", "Optimized for protein structures", protein_colors))
        
        # 8. Drug Discovery Scheme
        drug_colors = {
            'atom_h': '#ffffff',      # White
            'atom_c': '#ff6b6b',      # Red (carbon backbone)
            'atom_n': '#4ecdc4',      # Cyan (nitrogen)
            'atom_o': '#45b7d1',      # Blue (oxygen)
            'atom_f': '#96ceb4',      # Green (fluorine)
            'atom_p': '#feca57',      # Yellow (phosphorus)
            'atom_s': '#ff9ff3',      # Pink (sulfur)
            'atom_cl': '#54a0ff',     # Blue (chlorine)
            'atom_br': '#48dbfb',     # Light blue (bromine)
            'atom_i': '#a29bfe',      # Purple (iodine)
            'atom_selected': '#ee5a6f',    # Red
            'atom_highlight': '#feca57',   # Yellow
            'atom_positive': '#4ecdc4',    # Cyan
            'atom_negative': '#45b7d1',    # Blue
        }
        schemes.append(ColorScheme("Drug Discovery", "Optimized for drug molecules", drug_colors))
        
        return schemes
    
    def get_scheme_names(self) -> List[str]:
        """Get list of all scheme names."""
        return [scheme.name for scheme in self.schemes]
    
    def get_scheme(self, name: str) -> ColorScheme:
        """Get a specific color scheme by name."""
        for scheme in self.schemes:
            if scheme.name == name:
                return scheme
        return self.schemes[0]  # Return default if not found
    
    def apply_scheme(self, name: str):
        """Apply a color scheme by name."""
        scheme = self.get_scheme(name)
        scheme.apply()
        self.current_scheme = name
        return scheme.colors
    
    def get_current_scheme(self) -> str:
        """Get the current scheme name."""
        return self.current_scheme
    
    def get_scheme_description(self, name: str) -> str:
        """Get description of a scheme."""
        scheme = self.get_scheme(name)
        return scheme.description
    
    def cycle_scheme(self) -> Dict[str, str]:
        """Cycle to the next color scheme."""
        current_index = self.get_scheme_names().index(self.current_scheme)
        next_index = (current_index + 1) % len(self.schemes)
        next_name = self.get_scheme_names()[next_index]
        return self.apply_scheme(next_name)
    
    def add_custom_scheme(self, name: str, description: str, colors: Dict[str, str]):
        """Add a custom color scheme."""
        custom_scheme = ColorScheme(name, description, colors)
        self.schemes.append(custom_scheme)
        return custom_scheme
    
    def remove_scheme(self, name: str) -> bool:
        """Remove a custom color scheme."""
        if name == "Default":
            return False  # Cannot remove default scheme
        
        for i, scheme in enumerate(self.schemes):
            if scheme.name == name:
                del self.schemes[i]
                return True
        return False


# Global color scheme manager instance
color_scheme_manager = ColorSchemeManager()


def get_color_scheme_manager() -> ColorSchemeManager:
    """Get the global color scheme manager."""
    return color_scheme_manager


def apply_color_scheme(name: str) -> Dict[str, str]:
    """Apply a color scheme by name."""
    return color_scheme_manager.apply_scheme(name)


def cycle_color_scheme() -> Dict[str, str]:
    """Cycle to the next color scheme."""
    return color_scheme_manager.cycle_scheme()


def get_available_schemes() -> List[str]:
    """Get list of available color schemes."""
    return color_scheme_manager.get_scheme_names()
