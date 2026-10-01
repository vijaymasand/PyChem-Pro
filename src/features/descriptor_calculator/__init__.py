"""
Molecular Descriptor Calculator Module

This module provides comprehensive molecular descriptor calculation capabilities
including constitutional, topological, geometric, electronic, and custom descriptors.
"""

from .descriptor_engine import DescriptorEngine
from .descriptor_types import DescriptorCategory

# Optional GUI imports
try:
    from .gui import DescriptorCalculatorDialog
    _GUI_AVAILABLE = True
except ImportError:
    DescriptorCalculatorDialog = None
    _GUI_AVAILABLE = False

__all__ = [
    'DescriptorEngine',
    'DescriptorCategory', 
    'DescriptorCalculatorDialog',
    '_GUI_AVAILABLE'
]
