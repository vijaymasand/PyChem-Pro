from src.shared.ui.theme import COLORS

DISPLAY_RADIUS = {
    'H': 0.25, 'He': 0.31, 'C': 0.40, 'N': 0.38, 'O': 0.36, 'F': 0.32,
    'P': 0.44, 'S': 0.42, 'Cl': 0.39, 'Br': 0.41, 'I': 0.44, 'B': 0.38,
    'Si': 0.44, 'Se': 0.42, 'Na': 0.50, 'K': 0.55, 'Ca': 0.48, 'Fe': 0.44,
}

def _hex_to_rgb(hex_color: str):
    """Convert hex colour string to (r, g, b) tuple (0-255)."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))

def _hex_to_rgb_float(hex_color: str):
    """Convert hex colour string to (r, g, b) tuple in 0-1 range."""
    r, g, b = _hex_to_rgb(hex_color)
    return (r / 255.0, g / 255.0, b / 255.0)

def _element_color(symbol: str) -> tuple:
    """Return (r, g, b) in 0-255 for *symbol* using the Element colour table.

    Falls back to grey if the element lookup fails.
    """
    try:
        from src.core.domain.models.elements import get_element
        elem = get_element(symbol)
        if elem and elem.color:
            return _hex_to_rgb(elem.color)
    except Exception:
        pass
    return (180, 180, 180)

def _element_color_float(symbol: str) -> tuple:
    """Return (r, g, b) in 0.0-1.0 for *symbol*."""
    r, g, b = _element_color(symbol)
    return (r / 255.0, g / 255.0, b / 255.0)

def _display_radius(symbol: str) -> float:
    """Return display radius for *symbol* (Angstroms, aesthetic scale)."""
    return DISPLAY_RADIUS.get(symbol, 0.35)

