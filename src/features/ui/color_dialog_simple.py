"""
Simple color customization dialog without PySide6 dependencies.

Provides basic color scheme cycling without GUI components.
"""

from typing import Dict, List, Optional
from src.shared.ui.theme import COLORS
from src.features.ui.color_schemes import get_color_scheme_manager


class SimpleColorCustomizer:
    """
    Simple color customizer that works without GUI components.
    """
    
    def __init__(self):
        self.scheme_manager = get_color_scheme_manager()
        self.current_index = 0
        
    def get_available_schemes(self) -> List[str]:
        """Get list of available color schemes."""
        return self.scheme_manager.get_scheme_names()
    
    def get_current_scheme(self) -> str:
        """Get current color scheme name."""
        return self.scheme_manager.get_current_scheme()
    
    def cycle_scheme(self) -> Dict[str, str]:
        """Cycle to the next color scheme."""
        return self.scheme_manager.cycle_scheme()
    
    def apply_scheme(self, scheme_name: str) -> Dict[str, str]:
        """Apply a specific color scheme."""
        return self.scheme_manager.apply_scheme(scheme_name)
    
    def get_scheme_description(self, scheme_name: str) -> str:
        """Get description of a color scheme."""
        return self.scheme_manager.get_scheme_description(scheme_name)
    
    def get_current_colors(self) -> Dict[str, str]:
        """Get current color settings."""
        return {k: v for k, v in COLORS.items() if k.startswith('atom_')}
    
    def preview_scheme(self, scheme_name: str) -> Dict[str, str]:
        """Preview a color scheme without applying it."""
        scheme = self.scheme_manager.get_scheme(scheme_name)
        return scheme.colors.copy()


def create_simple_color_dialog():
    """Create a simple color dialog interface."""
    customizer = SimpleColorCustomizer()
    
    print("\n" + "="*50)
    print("SMILES Molecular Toolkit - Color Customization")
    print("="*50)
    
    schemes = customizer.get_available_schemes()
    current = customizer.get_current_scheme()
    
    print(f"Current scheme: {current}")
    print(f"Available schemes: {', '.join(schemes)}")
    print()
    
    print("Color Schemes:")
    print("-" * 30)
    
    for i, scheme in enumerate(schemes):
        marker = "→" if scheme == current else " "
        description = customizer.get_scheme_description(scheme)
        print(f"{marker} {i+1}. {scheme}: {description}")
    
    print()
    print("Usage:")
    print("- Click 'Colors' button to cycle through schemes")
    print("- Each click applies the next scheme in the list")
    print("- Current colors are updated in real-time")
    
    return customizer


def get_color_scheme_info() -> str:
    """Get information about current color scheme."""
    customizer = SimpleColorCustomizer()
    current = customizer.get_current_scheme()
    description = customizer.get_scheme_description(current)
    colors = customizer.get_current_colors()
    
    info = f"Current Scheme: {current}\n"
    info += f"Description: {description}\n"
    info += f"Colors: {len(colors)} atom colors defined\n"
    
    return info


def apply_next_color_scheme() -> Dict[str, str]:
    """Apply the next color scheme in the cycle."""
    customizer = SimpleColorCustomizer()
    return customizer.cycle_scheme()


def get_color_scheme_list() -> List[str]:
    """Get list of all available color schemes."""
    customizer = SimpleColorCustomizer()
    return customizer.get_available_schemes()
