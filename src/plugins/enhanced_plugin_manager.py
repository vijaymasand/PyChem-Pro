"""
Enhanced Plugin Manager for Runtime Plugin Support

This enhanced plugin manager supports both bundled (compiled) and external plugins,
enabling runtime extensibility even after the application is compiled into an executable.

Key Features:
- Hybrid plugin discovery (bundled + external)
- Runtime plugin installation and removal
- Plugin marketplace integration
- Safe plugin sandboxing
- Version compatibility checking
- Plugin dependency management
"""

import os
import sys
import json
import time
import shutil
import logging
import zipfile
import tempfile
import platform
import importlib.util
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable, Tuple
from dataclasses import dataclass, field
from urllib.parse import urlparse
from urllib.request import urlopen, Request

from src.shared.qt_compat import QStandardPaths

from .plugin_types import (
    PluginInfo, PluginType, PluginStatus, PluginMetadata,
    PluginLoadedEvent, PluginUnloadedEvent, PluginErrorEvent,
    MoleculeChangedEvent
)
from .base_plugin import BasePlugin
from .utils.safety import PluginSandbox, PluginValidator
from .utils.validation import validate_plugin_file
from . import PLUGIN_API_VERSION


@dataclass
class PluginSource:
    """Represents a plugin source location."""
    name: str
    url: str
    description: str
    trusted: bool = False


@dataclass
class PluginRepository:
    """Represents a plugin repository."""
    name: str
    url: str
    description: str = ""
    plugins: List[Dict[str, Any]] = field(default_factory=list)
    last_updated: float = 0.0


class EnhancedPluginManager:
    """
    Enhanced plugin manager with runtime extensibility support.
    
    This manager handles:
    - Bundled plugins (compiled with the app)
    - External plugins (user-installed)
    - Plugin repositories and marketplace
    - Safe plugin installation and removal
    - Plugin version and dependency management
    """
    
    def __init__(self, bundled_plugins_directory: str = None, user_plugins_directory: str = None):
        """
        Initialize the enhanced plugin manager.
        
        Args:
            bundled_plugins_directory: Path to bundled plugins (compiled with app)
            user_plugins_directory: Path to user plugins (external)
        """
        self.logger = logging.getLogger("enhanced_plugin.manager")
        
        # Set up plugin directories
        self._setup_plugin_directories(bundled_plugins_directory, user_plugins_directory)
        
        # Plugin storage
        self._bundled_plugins: Dict[str, PluginMetadata] = {}
        self._user_plugins: Dict[str, PluginMetadata] = {}
        self._loaded_plugins: Dict[str, BasePlugin] = {}
        self._plugin_widgets: Dict[str, Any] = {}
        
        # Plugin repositories
        self._repositories: Dict[str, PluginRepository] = {}
        self._plugin_sources: Dict[str, PluginSource] = {}
        
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
            'total_bundled': 0,
            'total_user': 0,
            'total_loaded': 0,
            'load_time': 0.0,
            'last_scan_time': 0.0
        }
        
        # Initialize default repositories
        self._initialize_default_repositories()
        
        self.logger.info(f"Enhanced plugin manager initialized")
        self.logger.info(f"Bundled plugins directory: {self.bundled_plugins_directory}")
        self.logger.info(f"User plugins directory: {self.user_plugins_directory}")
    
    def _setup_plugin_directories(self, bundled_dir: str = None, user_dir: str = None):
        """Set up plugin directories with proper fallbacks."""
        # Bundled plugins directory (compiled with app)
        if bundled_dir is None:
            # Check for Nuitka-specific indicators or relative paths
            if getattr(sys, 'frozen', False):
                # We are running in a bundled executable
                # Try multiple paths for Nuitka onefile/standalone compatibility
                candidates = []
                
                # 1. Nuitka onefile: data files extracted relative to module __file__
                module_dir = Path(__file__).parent.parent.parent
                candidates.append(module_dir / "plugins")
                
                # 2. Nuitka standalone: data dirs relative to exe
                exe_dir = Path(sys.executable).parent
                candidates.append(exe_dir / "plugins")
                
                # 3. macOS bundle structure
                if platform.system() == "Darwin" and ".app/Contents/MacOS" in str(exe_dir):
                    candidates.append(exe_dir.parent / "Resources" / "plugins")
                
                # Use the first candidate that exists
                bundled_path = None
                for cand in candidates:
                    if cand.exists():
                        bundled_path = cand
                        self.logger.info(f"Found bundled plugins at: {cand}")
                        break
                
                if bundled_path is None:
                    bundled_path = candidates[0]  # Default fallback
                    self.logger.warning(f"No bundled plugins directory found. Tried: {[str(c) for c in candidates]}")
                
                self.bundled_plugins_directory = bundled_path
            else:
                # We are running in development mode
                project_root = Path(__file__).parent.parent.parent
                self.bundled_plugins_directory = project_root / "plugins"
        else:
            self.bundled_plugins_directory = Path(bundled_dir)
        
        # User plugins directory (external, writable)
        if user_dir is None:
            # Use QStandardPaths for reliable cross-platform user data directory
            app_data = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
            if not app_data:
                # Fallback if QStandardPaths fails
                if platform.system() == "Windows":
                    app_data = os.environ.get("APPDATA", "")
                    self.user_plugins_directory = Path(app_data) / "PyChem" / "plugins"
                else:
                    self.user_plugins_directory = Path.home() / ".pychem" / "plugins"
            else:
                self.user_plugins_directory = Path(app_data) / "plugins"
        else:
            self.user_plugins_directory = Path(user_dir)
        
        # Ensure user directory exists (bundled may be read-only, so we just check)
        try:
            self.user_plugins_directory.mkdir(parents=True, exist_ok=True)
            if not self.bundled_plugins_directory.exists():
                self.logger.warning(f"Bundled plugins directory does not exist: {self.bundled_plugins_directory}")
        except Exception as e:
            self.logger.error(f"Error creating plugin directories: {e}")

    
    def _initialize_default_repositories(self):
        """Initialize default plugin repositories."""
        # Official PyChem repository
        self._repositories["official"] = PluginRepository(
            name="PyChem Official",
            url="https://plugins.pychem.org/api/plugins",
            description="Official PyChem plugin repository"
        )
        
        # Community repository
        self._repositories["community"] = PluginRepository(
            name="PyChem Community",
            url="https://community.pychem.org/plugins.json",
            description="Community-contributed plugins"
        )
    
    def set_api(self, api):
        """Set the plugin API for integration."""
        self._api = api
    
    # -------------------------------------------------------------------------
    # Plugin Discovery (Hybrid)
    # -------------------------------------------------------------------------
    
    def discover_all_plugins(self) -> Dict[str, PluginMetadata]:
        """
        Discover all plugins from both bundled and user directories.
        
        Returns:
            Dictionary of all discovered plugins
        """
        start_time = time.time()
        
        # Discover bundled plugins
        self._bundled_plugins = self.discover_bundled_plugins()
        
        # Discover user plugins
        self._user_plugins = self.discover_user_plugins()
        
        # Merge results (user plugins take precedence)
        all_plugins = {**self._bundled_plugins, **self._user_plugins}
        
        # Update statistics
        self._stats['total_bundled'] = len(self._bundled_plugins)
        self._stats['total_user'] = len(self._user_plugins)
        self._stats['last_scan_time'] = time.time() - start_time
        
        self.logger.info(f"Discovery complete: {len(self._bundled_plugins)} bundled, {len(self._user_plugins)} user plugins")
        
        return all_plugins
    
    def discover_bundled_plugins(self) -> Dict[str, PluginMetadata]:
        """Discover bundled plugins (compiled with the app)."""
        return self._discover_plugins_in_directory(self.bundled_plugins_directory, "bundled")
    
    def discover_user_plugins(self) -> Dict[str, PluginMetadata]:
        """Discover user plugins (externally installed)."""
        return self._discover_plugins_in_directory(self.user_plugins_directory, "user")
    
    def _discover_plugins_in_directory(self, directory: Path, plugin_type: str) -> Dict[str, PluginMetadata]:
        """Discover plugins in a specific directory."""
        discovered = {}
        
        try:
            self.logger.info(f"Discovering {plugin_type} plugins in: {directory}")
            if not directory.exists():
                self.logger.error(f"Directory does not exist: {directory}")
                return discovered
            
            # Look for Python files
            files = list(directory.glob("*.py"))
            self.logger.info(f"Found {len(files)} .py files in {directory}")
            
            for file_path in files:
                if file_path.name.startswith("__"):
                    continue
                
                try:
                    plugin_metadata = self._analyze_plugin_file(file_path, plugin_type)
                    if plugin_metadata:
                        discovered[plugin_metadata.info.name] = plugin_metadata
                        self.logger.info(f"Discovered {plugin_type} plugin: {plugin_metadata.info.name}")
                
                except Exception as e:
                    self.logger.error(f"Error analyzing {plugin_type} plugin {file_path}: {e}")
                    self._emit_event(PluginErrorEvent("discovery", str(e)))
            
            # Also look for plugin packages (directories with __init__.py)
            for package_path in directory.iterdir():
                if package_path.is_dir() and (package_path / "__init__.py").exists():
                    try:
                        plugin_metadata = self._analyze_plugin_package(package_path, plugin_type)
                        if plugin_metadata:
                            discovered[plugin_metadata.info.name] = plugin_metadata
                            self.logger.info(f"Discovered {plugin_type} plugin package: {plugin_metadata.info.name}")
                    
                    except Exception as e:
                        self.logger.error(f"Error analyzing {plugin_type} plugin package {package_path}: {e}")
        
        except Exception as e:
            self.logger.error(f"Error during {plugin_type} plugin discovery: {e}")
        
        return discovered
    
    def _analyze_plugin_file(self, file_path: Path, plugin_type: str) -> Optional[PluginMetadata]:
        """Analyze a plugin file and extract metadata."""
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
                temp_instance = plugin_class()
                info = temp_instance.info
                
                # Check API compatibility
                if not info.is_compatible_with_api(PLUGIN_API_VERSION):
                    self.logger.warning(f"Plugin {info.name} not compatible with API version {PLUGIN_API_VERSION}")
                    return None
                
                metadata = PluginMetadata(
                    info=info,
                    module_path=str(file_path),
                    class_name=plugin_class.__name__,
                    source=plugin_type
                )
                
                return metadata
                
            except Exception as e:
                self.logger.error(f"Error creating plugin metadata for {file_path}: {e}")
                return None
                
        except Exception as e:
            self.logger.error(f"Error analyzing plugin file {file_path}: {e}")
            return None
    
    def _analyze_plugin_package(self, package_path: Path, plugin_type: str) -> Optional[PluginMetadata]:
        """Analyze a plugin package (directory with __init__.py)."""
        try:
            init_file = package_path / "__init__.py"
            if not init_file.exists():
                return None
            
            # Import the package
            module_name = package_path.name
            spec = importlib.util.spec_from_file_location(module_name, init_file)
            if spec is None or spec.loader is None:
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
                return None
            
            # Use the first plugin class found
            plugin_class = plugin_classes[0]
            
            # Create metadata
            temp_instance = plugin_class()
            info = temp_instance.info
            
            metadata = PluginMetadata(
                info=info,
                module_path=str(package_path),
                class_name=plugin_class.__name__,
                source=plugin_type,
                is_package=True
            )
            
            return metadata
            
        except Exception as e:
            self.logger.error(f"Error analyzing plugin package {package_path}: {e}")
            return None
    
    # -------------------------------------------------------------------------
    # Plugin Installation and Removal
    # -------------------------------------------------------------------------
    
    def install_plugin_from_file(self, file_path: str) -> Tuple[bool, str]:
        """
        Install a plugin from a file.
        
        Args:
            file_path: Path to plugin file (.py or .zip)
            
        Returns:
            Tuple of (success, message)
        """
        try:
            file_path = Path(file_path)
            
            if file_path.suffix == '.py':
                return self._install_python_plugin(file_path)
            elif file_path.suffix == '.zip':
                return self._install_zip_plugin(file_path)
            else:
                return False, f"Unsupported file type: {file_path.suffix}"
                
        except Exception as e:
            error_msg = f"Error installing plugin: {e}"
            self.logger.error(error_msg)
            return False, error_msg
    
    def _install_python_plugin(self, file_path: Path) -> Tuple[bool, str]:
        """Install a Python plugin file."""
        try:
            # Validate the plugin first
            temp_metadata = self._analyze_plugin_file(file_path, "user")
            if not temp_metadata:
                return False, "Plugin validation failed"
            
            # Check for conflicts
            if temp_metadata.info.name in self._user_plugins:
                return False, f"Plugin '{temp_metadata.info.name}' is already installed"
            
            # Copy to user plugins directory
            dest_path = self.user_plugins_directory / file_path.name
            shutil.copy2(file_path, dest_path)
            
            # Add to user plugins
            self._user_plugins[temp_metadata.info.name] = temp_metadata
            
            success_msg = f"Plugin '{temp_metadata.info.name}' installed successfully"
            self.logger.info(success_msg)
            return True, success_msg
            
        except Exception as e:
            return False, f"Error installing Python plugin: {e}"
    
    def _install_zip_plugin(self, file_path: Path) -> Tuple[bool, str]:
        """Install a plugin from a ZIP file."""
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                # Extract ZIP
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    zip_ref.extractall(temp_dir)
                
                # Find plugin files
                temp_path = Path(temp_dir)
                plugin_files = list(temp_path.glob("*.py"))
                
                if not plugin_files:
                    return False, "No Python plugin files found in ZIP"
                
                # Install each plugin file
                installed_plugins = []
                for plugin_file in plugin_files:
                    success, msg = self._install_python_plugin(plugin_file)
                    if success:
                        installed_plugins.append(plugin_file.stem)
                    else:
                        return False, f"Failed to install {plugin_file.name}: {msg}"
                
                success_msg = f"Installed {len(installed_plugins)} plugins from ZIP"
                self.logger.info(success_msg)
                return True, success_msg
                
        except Exception as e:
            return False, f"Error installing ZIP plugin: {e}"
    
    def uninstall_plugin(self, plugin_name: str) -> Tuple[bool, str]:
        """
        Uninstall a user plugin.
        
        Args:
            plugin_name: Name of the plugin to uninstall
            
        Returns:
            Tuple of (success, message)
        """
        try:
            # Check if plugin exists and is user plugin
            if plugin_name not in self._user_plugins:
                return False, f"Plugin '{plugin_name}' not found or is a bundled plugin"
            
            # Unload plugin if loaded
            if plugin_name in self._loaded_plugins:
                self.unload_plugin(plugin_name)
            
            # Remove plugin file
            metadata = self._user_plugins[plugin_name]
            plugin_path = Path(metadata.module_path)
            
            if plugin_path.is_file():
                plugin_path.unlink()
            elif plugin_path.is_dir():
                shutil.rmtree(plugin_path)
            
            # Remove from user plugins
            del self._user_plugins[plugin_name]
            
            success_msg = f"Plugin '{plugin_name}' uninstalled successfully"
            self.logger.info(success_msg)
            return True, success_msg
            
        except Exception as e:
            error_msg = f"Error uninstalling plugin: {e}"
            self.logger.error(error_msg)
            return False, error_msg
    
    # -------------------------------------------------------------------------
    # Plugin Loading and Unloading
    # -------------------------------------------------------------------------
    
    def load_plugin(self, plugin_name: str) -> bool:
        """Load a specific plugin."""
        if plugin_name in self._loaded_plugins:
            self.logger.warning(f"Plugin {plugin_name} is already loaded")
            return True
        
        # Find plugin metadata
        metadata = self._user_plugins.get(plugin_name) or self._bundled_plugins.get(plugin_name)
        if not metadata:
            self.logger.error(f"Plugin {plugin_name} not found")
            return False
        
        start_time = time.time()
        
        try:
            self.logger.info(f"Loading plugin: {plugin_name}")
            metadata.status = PluginStatus.LOADING
            
            # Import the plugin module
            if metadata.is_package:
                spec = importlib.util.spec_from_file_location(
                    metadata.class_name, metadata.module_path + "/__init__.py"
                )
            else:
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
            
            # Initialize plugin
            try:
                import inspect
                init_sig = inspect.signature(plugin_instance.initialize)
                init_params = list(init_sig.parameters.keys())
                
                if len(init_params) == 0:
                    plugin_instance.initialize()
                elif len(init_params) >= 2 and self._api:
                    plugin_instance.initialize(self._api.main_window, self._api)
                elif len(init_params) == 1:
                    plugin_instance.initialize(self._api.main_window if self._api else None)
                else:
                    plugin_instance.initialize()
            except Exception as init_error:
                self.logger.warning(f"Plugin {plugin_name} initialization error: {init_error}")
            
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
        """Unload a specific plugin."""
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
            metadata = self._user_plugins.get(plugin_name) or self._bundled_plugins.get(plugin_name)
            if metadata:
                metadata.status = PluginStatus.INACTIVE
                metadata.instance = None
            
            # Update statistics
            self._stats['total_loaded'] = max(0, self._stats['total_loaded'] - 1)
            
            self.logger.info(f"Plugin {plugin_name} unloaded successfully")
            self._emit_event(PluginUnloadedEvent(plugin_name))
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error unloading plugin {plugin_name}: {e}")
            self._emit_event(PluginErrorEvent(plugin_name, str(e)))
            return False
    
    # -------------------------------------------------------------------------
    # Plugin Repository Management
    # -------------------------------------------------------------------------
    
    def refresh_repository(self, repo_name: str) -> Tuple[bool, str]:
        """Refresh a plugin repository."""
        try:
            if repo_name not in self._repositories:
                return False, f"Repository '{repo_name}' not found"
            
            repo = self._repositories[repo_name]
            
            # Fetch repository data
            req = Request(repo.url)
            req.add_header('User-Agent', 'PyChem-PluginManager/1.0')
            
            with urlopen(req) as response:
                data = json.loads(response.read().decode('utf-8'))
            
            # Update repository
            repo.plugins = data.get('plugins', [])
            repo.last_updated = time.time()
            
            success_msg = f"Repository '{repo_name}' refreshed successfully"
            self.logger.info(success_msg)
            return True, success_msg
            
        except Exception as e:
            error_msg = f"Error refreshing repository '{repo_name}': {e}"
            self.logger.error(error_msg)
            return False, error_msg
    
    def get_available_plugins(self, repo_name: str = None) -> List[Dict[str, Any]]:
        """Get available plugins from repositories."""
        if repo_name:
            if repo_name not in self._repositories:
                return []
            return self._repositories[repo_name].plugins
        else:
            all_plugins = []
            for repo in self._repositories.values():
                all_plugins.extend(repo.plugins)
            return all_plugins
    
    def install_plugin_from_repository(self, repo_name: str, plugin_id: str) -> Tuple[bool, str]:
        """Install a plugin from a repository."""
        try:
            if repo_name not in self._repositories:
                return False, f"Repository '{repo_name}' not found"
            
            repo = self._repositories[repo_name]
            
            # Find plugin in repository
            plugin_data = None
            for plugin in repo.plugins:
                if plugin.get('id') == plugin_id:
                    plugin_data = plugin
                    break
            
            if not plugin_data:
                return False, f"Plugin '{plugin_id}' not found in repository"
            
            # Download plugin
            download_url = plugin_data.get('download_url')
            if not download_url:
                return False, "Plugin has no download URL"
            
            # Download to temporary file
            with tempfile.NamedTemporaryFile(suffix='.py', delete=False) as temp_file:
                req = Request(download_url)
                req.add_header('User-Agent', 'PyChem-PluginManager/1.0')
                
                with urlopen(req) as response:
                    temp_file.write(response.read())
                
                temp_path = temp_file.name
            
            # Install from temporary file
            try:
                success, msg = self.install_plugin_from_file(temp_path)
                return success, msg
            finally:
                # Clean up temporary file
                try:
                    os.unlink(temp_path)
                except:
                    pass
            
        except Exception as e:
            error_msg = f"Error installing plugin from repository: {e}"
            self.logger.error(error_msg)
            return False, error_msg
    
    # -------------------------------------------------------------------------
    # Plugin Management Interface
    # -------------------------------------------------------------------------
    
    def get_plugin(self, plugin_name: str) -> Optional[BasePlugin]:
        """Get a loaded plugin instance."""
        return self._loaded_plugins.get(plugin_name)
    
    def get_plugin_widget(self, plugin_name: str) -> Optional[Any]:
        """Get the widget for a plugin."""
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
        """Get all discovered plugins (bundled + user)."""
        all_plugins = {**self._bundled_plugins, **self._user_plugins}
        return all_plugins.copy()
    
    def get_bundled_plugins(self) -> Dict[str, PluginMetadata]:
        """Get bundled plugins."""
        return self._bundled_plugins.copy()
    
    def get_user_plugins(self) -> Dict[str, PluginMetadata]:
        """Get user plugins."""
        return self._user_plugins.copy()
    
    def get_loaded_plugins(self) -> Dict[str, BasePlugin]:
        """Get all loaded plugins."""
        return self._loaded_plugins.copy()
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get plugin system statistics."""
        return self._stats.copy()

    def is_plugin_loaded(self, plugin_name: str) -> bool:
        """Check if a plugin is loaded."""
        return plugin_name in self._loaded_plugins
    
    def is_user_plugin(self, plugin_name: str) -> bool:
        """Check if a plugin is a user plugin."""
        return plugin_name in self._user_plugins
    
    def is_bundled_plugin(self, plugin_name: str) -> bool:
        """Check if a plugin is a bundled plugin."""
        return plugin_name in self._bundled_plugins
    
    # -------------------------------------------------------------------------
    # Event System
    # -------------------------------------------------------------------------
    
    def add_event_handler(self, event_type: str, handler: Callable):
        """Add an event handler."""
        if event_type not in self._event_handlers:
            self._event_handlers[event_type] = []
        self._event_handlers[event_type].append(handler)
    
    def remove_event_handler(self, event_type: str, handler: Callable):
        """Remove an event handler."""
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
        """Set the current molecule and notify all loaded plugins."""
        self._current_molecule = molecule
        self.on_molecule_changed(molecule)
    
    def on_molecule_changed(self, molecule):
        """Handle molecule changes by notifying all loaded plugins."""
        self._current_molecule = molecule
        
        event = MoleculeChangedEvent(molecule)
        self._emit_event(event)
        
        for plugin in self._loaded_plugins.values():
            try:
                plugin.on_molecule_changed(molecule)
            except Exception as e:
                self.logger.error(f"Error in plugin on_molecule_changed: {e}")
    
    # -------------------------------------------------------------------------
    # Statistics and Information
    # -------------------------------------------------------------------------
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get plugin manager statistics."""
        return self._stats.copy()
    
    def get_plugin_info(self, plugin_name: str) -> Optional[PluginMetadata]:
        """Get information about a specific plugin."""
        return self._user_plugins.get(plugin_name) or self._bundled_plugins.get(plugin_name)
    
    def get_repositories(self) -> Dict[str, PluginRepository]:
        """Get all repositories."""
        return self._repositories.copy()
    
    def add_repository(self, name: str, url: str, description: str = "") -> Tuple[bool, str]:
        """Add a new repository."""
        try:
            if name in self._repositories:
                return False, f"Repository '{name}' already exists"
            
            self._repositories[name] = PluginRepository(
                name=name,
                url=url,
                description=description
            )
            
            success_msg = f"Repository '{name}' added successfully"
            self.logger.info(success_msg)
            return True, success_msg
            
        except Exception as e:
            error_msg = f"Error adding repository: {e}"
            self.logger.error(error_msg)
            return False, error_msg
    
    def remove_repository(self, name: str) -> Tuple[bool, str]:
        """Remove a repository."""
        try:
            if name not in self._repositories:
                return False, f"Repository '{name}' not found"
            
            if name in ["official", "community"]:
                return False, "Cannot remove default repositories"
            
            del self._repositories[name]
            
            success_msg = f"Repository '{name}' removed successfully"
            self.logger.info(success_msg)
            return True, success_msg
            
        except Exception as e:
            error_msg = f"Error removing repository: {e}"
            self.logger.error(error_msg)
            return False, error_msg
