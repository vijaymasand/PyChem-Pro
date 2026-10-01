"""
Base Plugin Class and API

This module provides the base class that all plugins must inherit from,
along with the API that plugins can use to interact with the main application.

The plugin system is designed to be:
- Easy to use for plugin developers
- Safe and stable for the main application
- Flexible for different types of plugins
- Well-documented with examples
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, Callable
from dataclasses import dataclass
import logging

from .plugin_types import PluginInfo, PluginType, PluginEvent, MoleculeChangedEvent
from .utils.integration import PluginIntegrationAPI


# Plugin logger for developers
plugin_logger = logging.getLogger("plugin")


class BasePlugin(ABC):
    """
    Base class that all SMILES plugins must inherit from.
    
    This class provides the foundation for creating plugins that can
    integrate with the SMILES Molecular Toolkit. It handles the
    plugin lifecycle and provides access to the main application API.
    
    Plugin developers should inherit from this class and implement
    the required abstract methods.
    
    Example:
        class MyPlugin(BasePlugin):
            def __init__(self):
                super().__init__(
                    name="My Plugin",
                    version="1.0.0",
                    description="A sample plugin",
                    author="Your Name",
                    plugin_type=PluginType.ANALYSIS
                )
            
            def create_widget(self):
                # Create and return the plugin's main widget
                return MyWidget()
            
            def on_molecule_changed(self, molecule):
                # Handle molecule changes
                self.update_display(molecule)
    """
    
    def __init__(self, info: PluginInfo):
        """
        Initialize the plugin.
        
        Args:
            info: Plugin information including name, version, etc.
        """
        self._info = info
        self._api: Optional[PluginIntegrationAPI] = None
        self._main_window = None
        self._is_initialized = False
        self._settings = {}
        
        # Plugin logger (named after plugin)
        self.logger = logging.getLogger(f"plugin.{info.name.lower().replace(' ', '_')}")
        
        plugin_logger.info(f"Plugin initialized: {info.name} v{info.version}")
    
    @property
    def info(self) -> PluginInfo:
        """Get plugin information."""
        return self._info
    
    @property
    def api(self) -> Optional[PluginIntegrationAPI]:
        """Get the plugin API interface."""
        return self._api
    
    @property
    def main_window(self):
        """Get the main application window."""
        return self._main_window
    
    @property
    def is_initialized(self) -> bool:
        """Check if plugin is initialized."""
        return self._is_initialized
    
    @property
    def settings(self) -> Dict[str, Any]:
        """Get plugin settings."""
        return self._settings
    
    def log_info(self, message: str):
        """Log an informational message."""
        self.logger.info(message)
        
    def log_error(self, message: str):
        """Log an error message."""
        self.logger.error(message)
    
    # -------------------------------------------------------------------------
    # Abstract Methods (must be implemented by plugins)
    # -------------------------------------------------------------------------
    
    @abstractmethod
    def create_widget(self):
        """
        Create the plugin's main widget.
        
        This method is called when the plugin tab is created.
        The returned widget will be displayed in the plugin's tab.
        
        Returns:
            QWidget: The plugin's main widget
            
        Example:
            def create_widget(self):
                widget = QWidget()
                layout = QVBoxLayout(widget)
                
                # Add plugin UI elements
                self.label = QLabel("My Plugin")
                layout.addWidget(self.label)
                
                return widget
        """
        pass
    
    # -------------------------------------------------------------------------
    # Optional Override Methods (can be implemented by plugins)
    # -------------------------------------------------------------------------
    
    def initialize(self, main_window, api: PluginIntegrationAPI) -> bool:
        """
        Initialize the plugin with main application access.
        
        This method is called when the plugin is loaded and provides
        access to the main application and API.
        
        Args:
            main_window: The main application window
            api: Plugin API interface for interacting with the application
            
        Returns:
            bool: True if initialization successful, False otherwise
            
        Example:
            def initialize(self, main_window, api):
                self.main_window = main_window
                self.api = api
                
                # Setup plugin-specific initialization
                self.setup_ui()
                self.load_settings()
                
                return True
        """
        self._main_window = main_window
        self._api = api
        self._is_initialized = True
        
        self.logger.info(f"Plugin {self._info.name} initialized successfully")
        return True
    
    def cleanup(self):
        """
        Clean up plugin resources.
        
        This method is called when the plugin is unloaded.
        Plugins should override this to clean up any resources
        (threads, files, connections, etc.).
        
        Example:
            def cleanup(self):
                # Stop any running threads
                if hasattr(self, 'worker_thread'):
                    self.worker_thread.stop()
                
                # Save settings
                self.save_settings()
                
                # Clean up resources
                self.clear_data()
        """
        self.logger.info(f"Plugin {self._info.name} cleaned up")
        self._is_initialized = False
    
    def on_molecule_changed(self, molecule):
        """
        Handle molecule changes in the main application.
        
        This method is called whenever the current molecule changes
        in the main application. Plugins can override this to
        update their display or perform calculations.
        
        Args:
            molecule: The new current molecule (None if no molecule)
            
        Example:
            def on_molecule_changed(self, molecule):
                if molecule:
                    self.update_molecule_display(molecule)
                    self.calculate_properties(molecule)
                else:
                    self.clear_display()
        """
        pass
    
    def on_plugin_activated(self):
        """
        Handle plugin activation (when tab becomes visible).
        
        This method is called when the plugin's tab becomes active.
        Plugins can override this to refresh data or start operations.
        
        Example:
            def on_plugin_activated(self):
                # Refresh data when plugin becomes active
                self.refresh_data()
                self.update_status()
        """
        pass
    
    def on_plugin_deactivated(self):
        """
        Handle plugin deactivation (when tab becomes hidden).
        
        This method is called when the plugin's tab becomes inactive.
        Plugins can override this to pause operations or save state.
        
        Example:
            def on_plugin_deactivated(self):
                # Pause operations when plugin becomes inactive
                self.pause_calculations()
                self.save_current_state()
        """
        pass
    
    def get_settings(self) -> Dict[str, Any]:
        """
        Get plugin-specific settings.
        
        Returns:
            Dict containing plugin settings
        """
        return self._settings
    
    def set_settings(self, settings: Dict[str, Any]):
        """
        Set plugin-specific settings.
        
        Args:
            settings: Dictionary of plugin settings
        """
        self._settings.update(settings)
    
    # -------------------------------------------------------------------------
    # API Helper Methods (convenient access to common functionality)
    # -------------------------------------------------------------------------
    
    def get_current_molecule(self):
        """Get the current molecule from the main application."""
        if self._api:
            return self._api.get_current_molecule()
        return None
    
    def get_molecule_atoms(self):
        """Get atoms from the current molecule."""
        molecule = self.get_current_molecule()
        return molecule.atoms if molecule else []
    
    def get_molecule_bonds(self):
        """Get bonds from the current molecule."""
        molecule = self.get_current_molecule()
        return molecule.bonds if molecule else []
    
    def show_status_message(self, message: str, timeout: int = 3000):
        """Show a message in the main application status bar."""
        if self._api:
            self._api.show_status_message(message, timeout)
        else:
            self.logger.info(f"Status: {message}")
    
    def show_error_message(self, title: str, message: str):
        """Show an error message dialog."""
        if self._api:
            self._api.show_error_message(title, message)
        else:
            self.logger.error(f"{title}: {message}")
    
    def show_info_message(self, title: str, message: str):
        """Show an information message dialog."""
        if self._api:
            self._api.show_info_message(title, message)
        else:
            self.logger.info(f"{title}: {message}")
    
    def add_menu_item(self, menu_path: str, text: str, callback: Callable, shortcut: str = None):
        """
        Add a menu item to the main application.
        
        Args:
            menu_path: Path to menu (e.g., "Tools.My Plugin")
            text: Menu item text
            callback: Function to call when item is clicked
            shortcut: Keyboard shortcut (optional)
            
        Example:
            def add_menu_item(self):
                self.add_menu_item(
                    "Tools.My Plugin",
                    "Calculate My Property",
                    self.calculate_property,
                    "Ctrl+Shift+M"
                )
        """
        if self._api:
            self._api.add_menu_item(menu_path, text, callback, shortcut)
        else:
            self.logger.warning(f"Cannot add menu item: API not available")
    
    def get_viewer_2d(self):
        """Get the 2D molecular viewer."""
        if self._api:
            return self._api.get_viewer_2d()
        return None
    
    def get_viewer_3d(self):
        """Get the 3D molecular viewer."""
        if self._api:
            return self._api.get_viewer_3d()
        return None
    
    def update_viewer_2d(self):
        """Update the 2D molecular viewer."""
        if self._api:
            self._api.update_viewer_2d()
    
    def update_viewer_3d(self):
        """Update the 3D molecular viewer."""
        if self._api:
            self._api.update_viewer_3d()
    
    # -------------------------------------------------------------------------
    # Utility Methods
    # -------------------------------------------------------------------------
    
    def log_info(self, message: str):
        """Log an info message."""
        self.logger.info(message)
    
    def log_warning(self, message: str):
        """Log a warning message."""
        self.logger.warning(message)
    
    def log_error(self, message: str):
        """Log an error message."""
        self.logger.error(message)
    
    def log_debug(self, message: str):
        """Log a debug message."""
        self.logger.debug(message)


class PluginWidget:
    """
    Base class for plugin widgets.
    
    This class provides common functionality for plugin widgets
    and can be inherited by plugin developers for convenience.
    """
    
    def __init__(self, plugin: BasePlugin):
        """
        Initialize the plugin widget.
        
        Args:
            plugin: The plugin instance
        """
        self.plugin = plugin
        self.widget = None
    
    def get_widget(self):
        """Get the Qt widget."""
        return self.widget
    
    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        pass
    
    def cleanup(self):
        """Clean up widget resources."""
        pass
