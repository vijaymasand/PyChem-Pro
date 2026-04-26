"""
Plugin Manager

Handles plugin discovery, loading, unloading, and management.
Provides the central interface for the plugin system.
"""

import os
import sys
import time
import logging
import importlib.util
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass, field

from .plugin_types import (
    PluginInfo, PluginType, PluginStatus, PluginMetadata,
    PluginLoadedEvent, PluginUnloadedEvent, PluginErrorEvent,
    MoleculeChangedEvent
)
from .base_plugin import BasePlugin
from .utils.safety import PluginSandbox, PluginValidator
from .utils.validation import validate_plugin_file
from . import PLUGIN_API_VERSION


class PluginManager:
    """
    Central plugin management system.
    
    This class handles:
    - Plugin discovery and loading
    - Plugin lifecycle management
    - Plugin safety and validation
    - Event handling and communication
    - Plugin settings and configuration
    """
    
    def __init__(self, plugins_directory: str = None):
        """
        Initialize the plugin manager.
        
        Args:
            plugins_directory: Path to the plugins directory
        """
        self.logger = logging.getLogger("plugin.manager")
        
        # Set plugins directory
        if plugins_directory is None:
            # In a Nuitka-compiled app, data files sit next to the executable
            exe_dir = Path(sys.executable).parent
            if (exe_dir / "plugins").exists():
                self.plugins_directory = exe_dir / "plugins"
            else:
                project_root = Path(__file__).parent.parent.parent
                self.plugins_directory = project_root / "plugins"
        else:
            self.plugins_directory = Path(plugins_directory)
        
        # Ensure plugins directory exists
        self.plugins_directory.mkdir(exist_ok=True)
        
        # Plugin storage
        self._plugins: Dict[str, PluginMetadata] = {}
        self._loaded_plugins: Dict[str, BasePlugin] = {}
        self._plugin_widgets: Dict[str, Any] = {}
        
        # Safety and validation
        self.sandbox = PluginSandbox()
        self.validator = PluginValidator()
        
        # Event handlers
        self._event_handlers: Dict[str, List[Callable]] = {}
        
        # Plugin API for integration
        self._api = None
        
        # Current molecule for plugin initialization
        self._current_molecule = None
        
        # Statistics
        self._stats = {
            'total_discovered': 0,
            'total_loaded': 0,
            'load_time': 0.0,
            'last_scan_time': 0.0
        }
        
        self.logger.info(f"Plugin manager initialized with directory: {self.plugins_directory}")
    
    def set_api(self, api):
        """Set the plugin API for integration."""
        self._api = api
    
    # -------------------------------------------------------------------------
    # Plugin Discovery
    # -------------------------------------------------------------------------
    
    def discover_plugins(self) -> Dict[str, PluginMetadata]:
        """
        Discover all plugins in the plugins directory.
        
        Returns:
            Dictionary of discovered plugins
        """
        start_time = time.time()
        discovered = {}
        
        try:
            self.logger.info(f"Discovering plugins in: {self.plugins_directory}")
            
            # Look for Python files in plugins directory
            for file_path in self.plugins_directory.glob("*.py"):
                if file_path.name.startswith("__"):
                    continue  # Skip __init__.py and other special files
                
                try:
                    plugin_metadata = self._analyze_plugin_file(file_path)
                    if plugin_metadata:
                        discovered[plugin_metadata.info.name] = plugin_metadata
                        self.logger.info(f"Discovered plugin: {plugin_metadata.info.name}")
                
                except Exception as e:
                    self.logger.error(f"Error analyzing plugin {file_path}: {e}")
                    self._emit_event(PluginErrorEvent("discovery", str(e)))
            
            # Update statistics
            self._stats['total_discovered'] = len(discovered)
            self._stats['last_scan_time'] = time.time() - start_time
            
            # Store discovered plugins for later loading
            self._plugins.update(discovered)
            
            self.logger.info(f"Discovery complete: {len(discovered)} plugins found")
            
        except Exception as e:
            self.logger.error(f"Error during plugin discovery: {e}")
            self._emit_event(PluginErrorEvent("discovery", str(e)))
        
        return discovered
    
    def _analyze_plugin_file(self, file_path: Path) -> Optional[PluginMetadata]:
        """
        Analyze a plugin file and extract metadata.
        
        Args:
            file_path: Path to the plugin file
            
        Returns:
            PluginMetadata or None if analysis failed
        """
        try:
            # Validate file safety
            is_safe, safety_msg = self.sandbox.validate_plugin_file(str(file_path))
            if not is_safe:
                self.logger.warning(f"Plugin safety check failed for {file_path}: {safety_msg}")
                return None
            
            # Import the plugin module
            module_name = file_path.stem
            spec = importlib.util.spec_from_file_location(module_name, file_path)
            if spec is None or spec.loader is None:
                self.logger.error(f"Could not load spec for plugin: {file_path}")
                return None
            
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # Find plugin classes
            plugin_classes = []
            for name in dir(module):
                obj = getattr(module, name)
                if (isinstance(obj, type) and 
                    hasattr(obj, '__bases__') and 
                    any('BasePlugin' in str(base) for base in obj.__bases__)):
                    plugin_classes.append(obj)
            
            if not plugin_classes:
                self.logger.warning(f"No plugin classes found in: {file_path}")
                return None
            
            # Use the first plugin class found
            plugin_class = plugin_classes[0]
            
            # Validate plugin
            is_valid, issues = self.validator.validate_plugin(plugin_class, str(file_path))
            if not is_valid:
                self.logger.error(f"Plugin validation failed for {file_path}: {issues}")
                return None
            
            # Create metadata
            try:
                # Create temporary instance to get info
                temp_instance = plugin_class()
                info = temp_instance.info
                
                # Check API compatibility
                if not info.is_compatible_with_api(PLUGIN_API_VERSION):
                    self.logger.warning(f"Plugin {info.name} not compatible with API version {PLUGIN_API_VERSION}")
                    return None
                
                metadata = PluginMetadata(
                    info=info,
                    module_path=str(file_path),
                    class_name=plugin_class.__name__
                )
                
                return metadata
                
            except Exception as e:
                self.logger.error(f"Error creating plugin metadata for {file_path}: {e}")
                return None
                
        except Exception as e:
            self.logger.error(f"Error analyzing plugin file {file_path}: {e}")
            return None
    
    # -------------------------------------------------------------------------
    # Plugin Loading and Unloading
    # -------------------------------------------------------------------------
    
    def load_plugin(self, plugin_name: str) -> bool:
        """
        Load a specific plugin.
        
        Args:
            plugin_name: Name of the plugin to load
            
        Returns:
            True if loading successful, False otherwise
        """
        if plugin_name in self._loaded_plugins:
            self.logger.warning(f"Plugin {plugin_name} is already loaded")
            return True
        
        if plugin_name not in self._plugins:
            self.logger.error(f"Plugin {plugin_name} not found")
            return False
        
        metadata = self._plugins[plugin_name]
        start_time = time.time()
        
        try:
            self.logger.info(f"Loading plugin: {plugin_name}")
            metadata.status = PluginStatus.LOADING
            
            # Import the plugin module
            spec = importlib.util.spec_from_file_location(
                metadata.class_name, metadata.module_path
            )
            if spec is None or spec.loader is None:
                raise Exception("Could not load plugin spec")
            
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # Get the plugin class
            plugin_class = getattr(module, metadata.class_name)
            
            # Create plugin instance
            plugin_instance = plugin_class()
            
            # Initialize plugin - check signature and call appropriately
            try:
                import inspect
                init_sig = inspect.signature(plugin_instance.initialize)
                init_params = list(init_sig.parameters.keys())
                
                if len(init_params) == 0:  # Bound method with no additional args
                    # No additional arguments required
                    plugin_instance.initialize()
                elif len(init_params) >= 2 and self._api:  # Bound: main_window, api
                    plugin_instance.initialize(self._api.main_window, self._api)
                elif len(init_params) == 1:  # Bound: main_window
                    plugin_instance.initialize(self._api.main_window if self._api else None)
                else:
                    # Try without arguments as fallback
                    plugin_instance.initialize()
            except Exception as init_error:
                self.logger.warning(f"Plugin {plugin_name} initialization error: {init_error}")
                # Continue even if initialize fails - plugin may work without it
            
            # Store loaded plugin
            self._loaded_plugins[plugin_name] = plugin_instance
            metadata.instance = plugin_instance
            metadata.status = PluginStatus.ACTIVE
            metadata.load_time = time.time() - start_time
            
            # Update statistics
            self._stats['total_loaded'] += 1
            self._stats['load_time'] += metadata.load_time
            
            self.logger.info(f"Plugin {plugin_name} loaded successfully")
            self._emit_event(PluginLoadedEvent(plugin_name, metadata.info))
            
            # Notify plugin about current molecule if one exists
            if self._current_molecule is not None:
                try:
                    plugin_instance.on_molecule_changed(self._current_molecule)
                    self.logger.info(f"Notified {plugin_name} about current molecule")
                except Exception as e:
                    self.logger.error(f"Error notifying {plugin_name} about current molecule: {e}")
            
            return True
            
        except Exception as e:
            metadata.status = PluginStatus.ERROR
            metadata.error_message = str(e)
            self.logger.error(f"Error loading plugin {plugin_name}: {e}")
            self._emit_event(PluginErrorEvent(plugin_name, str(e)))
            return False
    
    def unload_plugin(self, plugin_name: str) -> bool:
        """
        Unload a specific plugin.
        
        Args:
            plugin_name: Name of the plugin to unload
            
        Returns:
            True if unloading successful, False otherwise
        """
        if plugin_name not in self._loaded_plugins:
            self.logger.warning(f"Plugin {plugin_name} is not loaded")
            return True
        
        try:
            self.logger.info(f"Unloading plugin: {plugin_name}")
            
            plugin = self._loaded_plugins[plugin_name]
            
            # Clean up plugin
            plugin.cleanup()
            
            # Remove from loaded plugins
            del self._loaded_plugins[plugin_name]
            
            # Remove widget if exists
            if plugin_name in self._plugin_widgets:
                widget = self._plugin_widgets[plugin_name]
                if hasattr(widget, 'cleanup'):
                    widget.cleanup()
                del self._plugin_widgets[plugin_name]
            
            # Update metadata
            if plugin_name in self._plugins:
                self._plugins[plugin_name].status = PluginStatus.INACTIVE
                self._plugins[plugin_name].instance = None
            
            # Update statistics
            self._stats['total_loaded'] = max(0, self._stats['total_loaded'] - 1)
            
            self.logger.info(f"Plugin {plugin_name} unloaded successfully")
            self._emit_event(PluginUnloadedEvent(plugin_name))
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error unloading plugin {plugin_name}: {e}")
            self._emit_event(PluginErrorEvent(plugin_name, str(e)))
            return False
    
    def reload_plugin(self, plugin_name: str) -> bool:
        """
        Reload a specific plugin.
        
        Args:
            plugin_name: Name of the plugin to reload
            
        Returns:
            True if reloading successful, False otherwise
        """
        self.logger.info(f"Reloading plugin: {plugin_name}")
        
        # Unload first
        if not self.unload_plugin(plugin_name):
            return False
        
        # Load again
        return self.load_plugin(plugin_name)
    
    # -------------------------------------------------------------------------
    # Plugin Management
    # -------------------------------------------------------------------------
    
    def get_plugin(self, plugin_name: str) -> Optional[BasePlugin]:
        """
        Get a loaded plugin instance.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            Plugin instance or None
        """
        return self._loaded_plugins.get(plugin_name)
    
    def get_plugin_widget(self, plugin_name: str) -> Optional[Any]:
        """
        Get the widget for a plugin.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            Plugin widget or None
        """
        if plugin_name not in self._plugin_widgets:
            plugin = self.get_plugin(plugin_name)
            if plugin:
                try:
                    widget = plugin.create_widget()
                    self._plugin_widgets[plugin_name] = widget
                except Exception as e:
                    self.logger.error(f"Error creating widget for plugin {plugin_name}: {e}")
                    return None
        
        return self._plugin_widgets.get(plugin_name)
    
    def get_all_plugins(self) -> Dict[str, PluginMetadata]:
        """Get all discovered plugins."""
        return self._plugins.copy()
    
    def get_loaded_plugins(self) -> Dict[str, BasePlugin]:
        """Get all loaded plugins."""
        return self._loaded_plugins.copy()
    
    def get_plugins_by_type(self, plugin_type: PluginType) -> List[str]:
        """
        Get plugins of a specific type.
        
        Args:
            plugin_type: Type of plugins to get
            
        Returns:
            List of plugin names
        """
        return [
            name for name, metadata in self._plugins.items()
            if metadata.info.plugin_type == plugin_type
        ]
    
    # -------------------------------------------------------------------------
    # Event System
    # -------------------------------------------------------------------------
    
    def add_event_handler(self, event_type: str, handler: Callable):
        """
        Add an event handler.
        
        Args:
            event_type: Type of event
            handler: Handler function
        """
        if event_type not in self._event_handlers:
            self._event_handlers[event_type] = []
        self._event_handlers[event_type].append(handler)
    
    def remove_event_handler(self, event_type: str, handler: Callable):
        """
        Remove an event handler.
        
        Args:
            event_type: Type of event
            handler: Handler function
        """
        if event_type in self._event_handlers:
            try:
                self._event_handlers[event_type].remove(handler)
            except ValueError:
                pass
    
    def _emit_event(self, event):
        """Emit an event to all handlers."""
        event_type = event.event_type
        if event_type in self._event_handlers:
            for handler in self._event_handlers[event_type]:
                try:
                    handler(event)
                except Exception as e:
                    self.logger.error(f"Error in event handler: {e}")
    
    # -------------------------------------------------------------------------
    # Plugin Lifecycle Events
    # -------------------------------------------------------------------------
    
    def set_current_molecule(self, molecule):
        """
        Set the current molecule and notify all loaded plugins.
        
        Args:
            molecule: Current molecule
        """
        self._current_molecule = molecule
        self.on_molecule_changed(molecule)
    
    def on_molecule_changed(self, molecule):
        """
        Handle molecule changes by notifying all loaded plugins.
        
        Args:
            molecule: New current molecule
        """
        # Store current molecule
        self._current_molecule = molecule
        
        event = MoleculeChangedEvent(molecule)
        self._emit_event(event)
        
        # Also call plugin methods directly
        for plugin in self._loaded_plugins.values():
            try:
                plugin.on_molecule_changed(molecule)
            except Exception as e:
                self.logger.error(f"Error in plugin {plugin.info.name} on_molecule_changed: {e}")
    
    def on_plugin_activated(self, plugin_name: str):
        """
        Handle plugin activation.
        
        Args:
            plugin_name: Name of the activated plugin
        """
        plugin = self.get_plugin(plugin_name)
        if plugin:
            try:
                plugin.on_plugin_activated()
            except Exception as e:
                self.logger.error(f"Error in plugin {plugin_name} on_plugin_activated: {e}")
    
    def on_plugin_deactivated(self, plugin_name: str):
        """
        Handle plugin deactivation.
        
        Args:
            plugin_name: Name of the deactivated plugin
        """
        plugin = self.get_plugin(plugin_name)
        if plugin:
            try:
                plugin.on_plugin_deactivated()
            except Exception as e:
                self.logger.error(f"Error in plugin {plugin_name} on_plugin_deactivated: {e}")
    
    # -------------------------------------------------------------------------
    # Statistics and Information
    # -------------------------------------------------------------------------
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get plugin manager statistics."""
        return self._stats.copy()
    
    def get_plugin_info(self, plugin_name: str) -> Optional[PluginMetadata]:
        """Get information about a specific plugin."""
        return self._plugins.get(plugin_name)
    
    def is_plugin_loaded(self, plugin_name: str) -> bool:
        """Check if a plugin is loaded."""
        return plugin_name in self._loaded_plugins
    
    def get_load_order(self) -> List[str]:
        """Get the recommended load order for plugins."""
        # Sort by type and dependencies
        plugins = list(self._plugins.values())
        
        # Define type priority
        type_priority = {
            PluginType.UTILITY: 0,
            PluginType.FILE_FORMAT: 1,
            PluginType.CALCULATION: 2,
            PluginType.ANALYSIS: 3,
            PluginType.VISUALIZATION: 4,
            PluginType.INTERFACE: 5
        }
        
        plugins.sort(key=lambda p: type_priority.get(p.info.plugin_type, 99))
        
        return [p.info.name for p in plugins]
