"""
Plugin Type Definitions and Enums

Defines the types of plugins supported by the SMILES Molecular Toolkit
and related data structures for plugin management.
"""

from enum import Enum
from dataclasses import dataclass
from typing import Dict, Any, List, Optional


class PluginType(Enum):
    """Enumeration of supported plugin types."""
    
    ANALYSIS = "analysis"
    """Plugins that perform molecular analysis (e.g., descriptor calculation)."""
    
    VISUALIZATION = "visualization"
    """Plugins that provide new molecular visualization methods."""
    
    FILE_FORMAT = "file_format"
    """Plugins that add support for new file formats (import/export)."""
    
    CALCULATION = "calculation"
    """Plugins that calculate molecular properties or perform computations."""
    
    UTILITY = "utility"
    """General utility plugins that don't fit in other categories."""
    
    INTERFACE = "interface"
    """Plugins that modify or extend the user interface."""


class PluginStatus(Enum):
    """Enumeration of plugin status states."""
    
    INACTIVE = "inactive"
    """Plugin is not loaded."""
    
    LOADING = "loading"
    """Plugin is currently being loaded."""
    
    ACTIVE = "active"
    """Plugin is loaded and functioning."""
    
    ERROR = "error"
    """Plugin encountered an error during loading or execution."""
    
    DISABLED = "disabled"
    """Plugin is manually disabled by user."""


@dataclass
class PluginInfo:
    """Information about a plugin."""
    
    name: str
    """Plugin name (must be unique)."""
    
    version: str
    """Plugin version (semantic versioning recommended)."""
    
    description: str
    """Brief description of what the plugin does."""
    
    author: str
    """Plugin author name."""
    
    plugin_type: PluginType
    """Type of plugin (analysis, visualization, etc.)."""
    
    min_api_version: str = "1.0.0"
    """Minimum plugin API version required."""
    
    max_api_version: Optional[str] = None
    """Maximum plugin API version supported (None for no upper limit)."""
    
    dependencies: List[str] = None
    """List of required dependencies (pip packages)."""
    
    homepage: Optional[str] = None
    """Plugin homepage URL."""
    
    license: Optional[str] = None
    """Plugin license."""
    
    keywords: List[str] = None
    """Keywords for plugin search and categorization."""
    
    def __post_init__(self):
        """Initialize default values."""
        if self.dependencies is None:
            self.dependencies = []
        if self.keywords is None:
            self.keywords = []
    
    def is_compatible_with_api(self, api_version: str) -> bool:
        """Check if plugin is compatible with given API version."""
        from packaging import version
        
        # Check minimum version
        if version.parse(api_version) < version.parse(self.min_api_version):
            return False
        
        # Check maximum version if specified
        if self.max_api_version and version.parse(api_version) > version.parse(self.max_api_version):
            return False
        
        return True


@dataclass
class PluginMetadata:
    """Metadata for a loaded plugin instance."""
    
    info: PluginInfo
    """Plugin information."""
    
    module_path: str
    """Path to the plugin module file."""
    
    class_name: str
    """Name of the plugin class."""
    
    instance: Optional['BasePlugin'] = None
    """Plugin instance (None if not loaded)."""
    
    source: Optional[str] = None
    """Source of the plugin ('bundled' or 'user')."""
    
    is_package: bool = False
    """True if the plugin is a package (directory), False if it's a single file."""
    
    status: PluginStatus = PluginStatus.INACTIVE
    """Current plugin status."""
    
    error_message: Optional[str] = None
    """Error message if plugin failed to load."""
    
    load_time: Optional[float] = None
    """Time taken to load the plugin (seconds)."""
    
    last_used: Optional[float] = None
    """Timestamp when plugin was last used."""


class PluginEvent:
    """Base class for plugin events."""
    
    def __init__(self, event_type: str, data: Dict[str, Any] = None):
        self.event_type = event_type
        self.data = data or {}
        self.timestamp = None


class PluginLoadedEvent(PluginEvent):
    """Event fired when a plugin is loaded."""
    
    def __init__(self, plugin_name: str, plugin_info: PluginInfo):
        super().__init__("plugin_loaded", {
            "plugin_name": plugin_name,
            "plugin_info": plugin_info
        })


class PluginUnloadedEvent(PluginEvent):
    """Event fired when a plugin is unloaded."""
    
    def __init__(self, plugin_name: str):
        super().__init__("plugin_unloaded", {
            "plugin_name": plugin_name
        })


class PluginErrorEvent(PluginEvent):
    """Event fired when a plugin encounters an error."""
    
    def __init__(self, plugin_name: str, error_message: str):
        super().__init__("plugin_error", {
            "plugin_name": plugin_name,
            "error_message": error_message
        })


class MoleculeChangedEvent(PluginEvent):
    """Event fired when the current molecule changes."""
    
    def __init__(self, molecule):
        super().__init__("molecule_changed", {
            "molecule": molecule
        })
