"""
Standalone color customization without any GUI dependencies.

Provides color scheme cycling and management without requiring PySide6.
"""

from typing import Dict, List, Optional
from src.shared.ui.theme import COLORS


class StandaloneColorManager:
    """
    Standalone color manager that works without any GUI dependencies.
    """
    
    def __init__(self):
        self.current_index = 0
        self.color_schemes = self._create_color_schemes()
        
    def _create_color_schemes(self) -> List[Dict[str, str]]:
        """Create standalone color schemes."""
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
        schemes.append(("Default", "Standard CPK colors", default_colors))
        
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
        schemes.append(("Pastel", "Soft, easy-on-the-eyes colors", pastel_colors))
        
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
        schemes.append(("Vibrant", "High contrast vivid colors", vibrant_colors))
        
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
        schemes.append(("Earth Tones", "Natural earth colors", earth_colors))
        
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
        schemes.append(("High Contrast", "Maximum contrast colors", high_contrast_colors))
        
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
        schemes.append(("Monochrome", "Grayscale colors", monochrome_colors))
        
        # 7. Protein Scheme
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
        schemes.append(("Protein", "Optimized for protein structures", protein_colors))
        
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
        schemes.append(("Drug Discovery", "Optimized for drug molecules", drug_colors))
        
        return schemes
    
    def cycle_scheme(self) -> Dict[str, str]:
        """Cycle to the next color scheme."""
        if not self.color_schemes:
            return {}
        
        # Get current scheme
        name, description, colors = self.color_schemes[self.current_index]
        
        # Apply current scheme
        COLORS.update(colors)
        
        # Move to next scheme
        self.current_index = (self.current_index + 1) % len(self.color_schemes)
        
        return colors.copy()
    
    def get_current_scheme_info(self) -> tuple:
        """Get current scheme name and description."""
        if not self.color_schemes:
            return "Unknown", "No schemes available"
        
        name, description, colors = self.color_schemes[self.current_index]
        return name, description
    
    def get_all_scheme_names(self) -> List[str]:
        """Get list of all scheme names."""
        return [scheme[0] for scheme in self.color_schemes]
    
    def get_scheme_count(self) -> int:
        """Get number of available schemes."""
        return len(self.color_schemes)
    
    def apply_scheme_by_name(self, name: str) -> bool:
        """Apply a specific scheme by name."""
        for i, (scheme_name, description, colors) in enumerate(self.color_schemes):
            if scheme_name == name:
                COLORS.update(colors)
                self.current_index = i
                return True
        return False
    
    def get_current_colors(self) -> Dict[str, str]:
        """Get current color settings."""
        return {k: v for k, v in COLORS.items() if k.startswith('atom_')}
    
    def print_scheme_info(self):
        """Print information about all schemes."""
        print("\n" + "="*60)
        print("SMILES Molecular Toolkit - Color Schemes")
        print("="*60)
        
        for i, (name, description, colors) in enumerate(self.color_schemes):
            if i == self.current_index:
                marker = "->"
            else:
                marker = "  "
            print(f"{marker} {i+1}. {name}")
            print(f"     {description}")
            print(f"     Colors: {len(colors)} atom colors defined")
        
        print(f"\nTotal schemes: {len(self.color_schemes)}")
        print("Click 'Colors' button to cycle through schemes")
        print("="*60)


# Global standalone color manager
_standalone_manager = None

def get_standalone_manager():
    """Get the global standalone color manager."""
    global _standalone_manager
    if _standalone_manager is None:
        _standalone_manager = StandaloneColorManager()
    return _standalone_manager


def cycle_standalone_colors() -> Dict[str, str]:
    """Cycle through standalone color schemes."""
    manager = get_standalone_manager()
    return manager.cycle_scheme()


def get_current_scheme_info() -> tuple:
    """Get current scheme information."""
    manager = get_standalone_manager()
    return manager.get_current_scheme_info()


def print_all_schemes():
    """Print all available color schemes."""
    manager = get_standalone_manager()
    manager.print_scheme_info()
