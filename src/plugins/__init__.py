"""
SMILES Molecular Toolkit Plugin System

A comprehensive plugin system that allows users to extend the application
with custom functionality through drop-in Python modules.

Features:
- Tab-based plugin integration
- Automatic plugin discovery
- Safe plugin sandboxing
- Hot loading/unloading
- Comprehensive API for plugin development
"""

__version__ = "1.0.0"
__author__ = "SMILES Development Team"

# Plugin system version for compatibility checking
PLUGIN_API_VERSION = "1.0.0"

# Export main components
from .base_plugin import BasePlugin, PluginType, PluginInfo
from .plugin_manager import PluginManager
from .plugin_types import *

__all__ = [
    'BasePlugin',
    'PluginType', 
    'PluginInfo',
    'PluginManager',
    'PLUGIN_API_VERSION'
]
