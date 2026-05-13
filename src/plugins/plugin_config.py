"""
Plugin Configuration System

Manages plugin settings, dependencies, and configuration files.
Provides persistent storage for plugin preferences and enables
proper dependency resolution.
"""

import os
import json
import logging
import platform
import time
from pathlib import Path
from typing import Dict, List, Optional, Any, Set, Tuple
from dataclasses import dataclass, field, asdict

from src.shared.qt_compat import QStandardPaths

# Simple version comparison without external dependency
def parse_version(version_string: str):
    """Parse a version string into comparable parts."""
    parts = version_string.replace('v', '').split('.')
    try:
        return tuple(int(p) for p in parts)
    except ValueError:
        return (0, 0, 0)

def compare_versions(v1: str, v2: str) -> int:
    """Compare two version strings. Returns -1 if v1 < v2, 0 if v1 == v2, 1 if v1 > v2."""
    v1_parts = parse_version(v1)
    v2_parts = parse_version(v2)
    
    for a, b in zip(v1_parts, v2_parts):
        if a < b:
            return -1
        elif a > b:
            return 1
    
    if len(v1_parts) < len(v2_parts):
        return -1
    elif len(v1_parts) > len(v2_parts):
        return 1
    
    return 0


@dataclass
class PluginDependency:
    """Represents a plugin dependency."""
    name: str
    version_requirement: str  # e.g., ">=1.0.0", "~2.1.0", "exact:1.2.3"
    optional: bool = False
    description: str = ""


@dataclass
class PluginConfig:
    """Configuration for a single plugin."""
    name: str
    enabled: bool = True
    auto_load: bool = True
    settings: Dict[str, Any] = field(default_factory=dict)
    dependencies: List[PluginDependency] = field(default_factory=list)
    load_order: int = 0
    last_updated: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'name': self.name,
            'enabled': self.enabled,
            'auto_load': self.auto_load,
            'settings': self.settings,
            'dependencies': [asdict(dep) for dep in self.dependencies],
            'load_order': self.load_order,
            'last_updated': self.last_updated
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PluginConfig':
        """Create from dictionary."""
        dependencies = []
        for dep_data in data.get('dependencies', []):
            dependencies.append(PluginDependency(**dep_data))
        
        return cls(
            name=data['name'],
            enabled=data.get('enabled', True),
            auto_load=data.get('auto_load', True),
            settings=data.get('settings', {}),
            dependencies=dependencies,
            load_order=data.get('load_order', 0),
            last_updated=data.get('last_updated', 0.0)
        )


@dataclass
class RepositoryConfig:
    """Configuration for a plugin repository."""
    name: str
    url: str
    description: str = ""
    trusted: bool = False
    auto_refresh: bool = False
    refresh_interval: int = 86400  # 24 hours in seconds
    last_refreshed: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'RepositoryConfig':
        """Create from dictionary."""
        return cls(**data)


class PluginConfigManager:
    """
    Manages plugin configuration files and settings.
    
    Handles:
    - Plugin settings persistence
    - Dependency management
    - Repository configuration
    - Plugin load order
    - Configuration validation
    """
    
    def __init__(self, config_dir: str = None):
        """
        Initialize the configuration manager.
        
        Args:
            config_dir: Directory to store configuration files
        """
        self.logger = logging.getLogger("plugin.config")
        
        # Set up configuration directory
        if config_dir is None:
            # Use QStandardPaths for reliable cross-platform user config directory
            app_data = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
            if not app_data:
                # Fallback if QStandardPaths fails
                if platform.system() == "Windows":
                    app_data = os.environ.get("APPDATA", "")
                    self.config_dir = Path(app_data) / "PyChem" / "config"
                else:
                    self.config_dir = Path.home() / ".pychem" / "config"
            else:
                self.config_dir = Path(app_data) / "config"
        else:
            self.config_dir = Path(config_dir)
        
        # Ensure directory exists
        try:
            self.config_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            self.logger.error(f"Error creating config directory: {e}")

        
        # Configuration file paths
        self.plugins_config_file = self.config_dir / "plugins.json"
        self.repositories_config_file = self.config_dir / "repositories.json"
        self.global_config_file = self.config_dir / "global.json"
        
        # In-memory configuration
        self._plugin_configs: Dict[str, PluginConfig] = {}
        self._repository_configs: Dict[str, RepositoryConfig] = {}
        self._global_config: Dict[str, Any] = {}
        
        # Load existing configurations
        self.load_configurations()
        
        self.logger.info(f"Plugin config manager initialized with directory: {self.config_dir}")
    
    def load_configurations(self):
        """Load all configuration files."""
        try:
            self._load_plugin_configs()
            self._load_repository_configs()
            self._load_global_config()
            self.logger.info("Configurations loaded successfully")
        except Exception as e:
            self.logger.error(f"Error loading configurations: {e}")
    
    def _load_plugin_configs(self):
        """Load plugin configurations."""
        if not self.plugins_config_file.exists():
            return
        
        try:
            with open(self.plugins_config_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            self._plugin_configs = {}
            for plugin_name, plugin_data in data.items():
                self._plugin_configs[plugin_name] = PluginConfig.from_dict(plugin_data)
                
        except Exception as e:
            self.logger.error(f"Error loading plugin configs: {e}")
    
    def _load_repository_configs(self):
        """Load repository configurations."""
        if not self.repositories_config_file.exists():
            # Create default repository configs
            self._create_default_repository_configs()
            return
        
        try:
            with open(self.repositories_config_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            self._repository_configs = {}
            for repo_name, repo_data in data.items():
                self._repository_configs[repo_name] = RepositoryConfig.from_dict(repo_data)
                
        except Exception as e:
            self.logger.error(f"Error loading repository configs: {e}")
            self._create_default_repository_configs()
    
    def _create_default_repository_configs(self):
        """Create default repository configurations."""
        self._repository_configs = {
            "official": RepositoryConfig(
                name="PyChem Official",
                url="https://plugins.pychem.org/api/plugins",
                description="Official PyChem plugin repository",
                trusted=True,
                auto_refresh=True
            ),
            "community": RepositoryConfig(
                name="PyChem Community",
                url="https://community.pychem.org/plugins.json",
                description="Community-contributed plugins",
                trusted=False,
                auto_refresh=False
            )
        }
        self.save_repository_configs()
    
    def _load_global_config(self):
        """Load global configuration."""
        if not self.global_config_file.exists():
            self._create_default_global_config()
            return
        
        try:
            with open(self.global_config_file, 'r', encoding='utf-8') as f:
                self._global_config = json.load(f)
                
        except Exception as e:
            self.logger.error(f"Error loading global config: {e}")
            self._create_default_global_config()
    
    def _create_default_global_config(self):
        """Create default global configuration."""
        self._global_config = {
            "auto_load_plugins": True,
            "check_updates": True,
            "update_check_interval": 86400,  # 24 hours
            "enable_plugin_sandbox": True,
            "max_plugin_load_time": 30.0,  # seconds
            "log_level": "INFO",
            "plugin_cache_size": 100,
            "last_update_check": 0.0
        }
        self.save_global_config()
    
    # -------------------------------------------------------------------------
    # Plugin Configuration Management
    # -------------------------------------------------------------------------
    
    def get_plugin_config(self, plugin_name: str) -> Optional[PluginConfig]:
        """Get configuration for a specific plugin."""
        return self._plugin_configs.get(plugin_name)
    
    def set_plugin_config(self, config: PluginConfig):
        """Set configuration for a plugin."""
        self._plugin_configs[config.name] = config
        self.save_plugin_configs()
    
    def update_plugin_setting(self, plugin_name: str, key: str, value: Any):
        """Update a specific setting for a plugin."""
        if plugin_name not in self._plugin_configs:
            self._plugin_configs[plugin_name] = PluginConfig(name=plugin_name)
        
        self._plugin_configs[plugin_name].settings[key] = value
        self._plugin_configs[plugin_name].last_updated = time.time()
        self.save_plugin_configs()
    
    def get_plugin_setting(self, plugin_name: str, key: str, default: Any = None) -> Any:
        """Get a specific setting for a plugin."""
        config = self._plugin_configs.get(plugin_name)
        if config:
            return config.settings.get(key, default)
        return default
    
    def enable_plugin(self, plugin_name: str):
        """Enable a plugin."""
        if plugin_name not in self._plugin_configs:
            self._plugin_configs[plugin_name] = PluginConfig(name=plugin_name)
        
        self._plugin_configs[plugin_name].enabled = True
        self.save_plugin_configs()
    
    def disable_plugin(self, plugin_name: str):
        """Disable a plugin."""
        if plugin_name not in self._plugin_configs:
            self._plugin_configs[plugin_name] = PluginConfig(name=plugin_name)
        
        self._plugin_configs[plugin_name].enabled = False
        self.save_plugin_configs()
    
    def is_plugin_enabled(self, plugin_name: str) -> bool:
        """Check if a plugin is enabled."""
        config = self._plugin_configs.get(plugin_name)
        return config.enabled if config else True  # Default to enabled
    
    def set_auto_load(self, plugin_name: str, auto_load: bool):
        """Set whether a plugin should auto-load."""
        if plugin_name not in self._plugin_configs:
            self._plugin_configs[plugin_name] = PluginConfig(name=plugin_name)
        
        self._plugin_configs[plugin_name].auto_load = auto_load
        self.save_plugin_configs()
    
    def should_auto_load(self, plugin_name: str) -> bool:
        """Check if a plugin should auto-load."""
        config = self._plugin_configs.get(plugin_name)
        return config.auto_load if config else True  # Default to auto-load
    
    def get_all_plugin_configs(self) -> Dict[str, PluginConfig]:
        """Get all plugin configurations."""
        return self._plugin_configs.copy()
    
    def remove_plugin_config(self, plugin_name: str):
        """Remove configuration for a plugin."""
        if plugin_name in self._plugin_configs:
            del self._plugin_configs[plugin_name]
            self.save_plugin_configs()
    
    def save_plugin_configs(self):
        """Save plugin configurations to file."""
        try:
            data = {name: config.to_dict() for name, config in self._plugin_configs.items()}
            with open(self.plugins_config_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            self.logger.error(f"Error saving plugin configs: {e}")
    
    # -------------------------------------------------------------------------
    # Dependency Management
    # -------------------------------------------------------------------------
    
    def add_plugin_dependency(self, plugin_name: str, dependency: PluginDependency):
        """Add a dependency to a plugin."""
        if plugin_name not in self._plugin_configs:
            self._plugin_configs[plugin_name] = PluginConfig(name=plugin_name)
        
        # Check if dependency already exists
        for existing_dep in self._plugin_configs[plugin_name].dependencies:
            if existing_dep.name == dependency.name:
                # Update existing dependency
                existing_dep.version_requirement = dependency.version_requirement
                existing_dep.optional = dependency.optional
                existing_dep.description = dependency.description
                self.save_plugin_configs()
                return
        
        # Add new dependency
        self._plugin_configs[plugin_name].dependencies.append(dependency)
        self.save_plugin_configs()
    
    def remove_plugin_dependency(self, plugin_name: str, dependency_name: str):
        """Remove a dependency from a plugin."""
        if plugin_name not in self._plugin_configs:
            return
        
        dependencies = self._plugin_configs[plugin_name].dependencies
        self._plugin_configs[plugin_name].dependencies = [
            dep for dep in dependencies if dep.name != dependency_name
        ]
        self.save_plugin_configs()
    
    def get_plugin_dependencies(self, plugin_name: str) -> List[PluginDependency]:
        """Get all dependencies for a plugin."""
        config = self._plugin_configs.get(plugin_name)
        return config.dependencies.copy() if config else []
    
    def resolve_dependencies(self, plugin_name: str) -> Tuple[List[str], List[str]]:
        """
        Resolve plugin dependencies.
        
        Returns:
            Tuple of (satisfied_dependencies, missing_dependencies)
        """
        config = self._plugin_configs.get(plugin_name)
        if not config:
            return [], []
        
        satisfied = []
        missing = []
        
        for dep in config.dependencies:
            if dep.name in self._plugin_configs:
                # Check version requirement
                dep_config = self._plugin_configs[dep.name]
                if self._check_version_requirement(dep_config.settings.get('version', '1.0.0'), 
                                                   dep.version_requirement):
                    satisfied.append(dep.name)
                else:
                    missing.append(f"{dep.name} (version {dep.version_requirement})")
            else:
                if not dep.optional:
                    missing.append(dep.name)
        
        return satisfied, missing
    
    def _check_version_requirement(self, version: str, requirement: str) -> bool:
        """
        Check if a version satisfies a requirement.
        
        Args:
            version: Current version string
            requirement: Version requirement string (e.g., ">=1.0.0", "~2.1.0")
            
        Returns:
            True if version satisfies requirement
        """
        try:
            v_parts = parse_version(version)
            
            if requirement.startswith("exact:"):
                required_version = requirement[6:]
                required_parts = parse_version(required_version)
                return v_parts == required_parts
            elif requirement.startswith(">="):
                required_version = requirement[2:]
                return compare_versions(version, required_version) >= 0
            elif requirement.startswith("<="):
                required_version = requirement[2:]
                return compare_versions(version, required_version) <= 0
            elif requirement.startswith(">"):
                required_version = requirement[1:]
                return compare_versions(version, required_version) > 0
            elif requirement.startswith("<"):
                required_version = requirement[1:]
                return compare_versions(version, required_version) < 0
            elif requirement.startswith("~"):
                required_version = requirement[1:]
                required_parts = parse_version(required_version)
                return v_parts[0] == required_parts[0] and v_parts[1] == required_parts[1]
            elif requirement.startswith("^"):
                required_version = requirement[1:]
                required_parts = parse_version(required_version)
                return v_parts[0] == required_parts[0]
            else:
                required_version = requirement
                return compare_versions(version, required_version) == 0
        except Exception as e:
            self.logger.error(f"Error checking version requirement: {e}")
            return False
    
    def _create_default_global_config(self):
        """Create default global configuration."""
        self._global_config = {
            "auto_load_plugins": False,
            "check_updates": True,
            "update_check_interval": 86400,  # 24 hours
            "enable_plugin_sandbox": True,
            "max_plugin_load_time": 30.0,  # seconds
            "log_level": "INFO",
            "plugin_cache_size": 100,
            "last_update_check": 0.0
        }
        self.save_global_config()

    # -------------------------------------------------------------------------
    # Plugin Configuration Management
    # -------------------------------------------------------------------------

    def get_plugin_config(self, plugin_name: str) -> Optional[PluginConfig]:
        """Get configuration for a specific plugin."""
        return self._plugin_configs.get(plugin_name)

    def set_plugin_config(self, config: PluginConfig):
        """Set configuration for a plugin."""
        self._plugin_configs[config.name] = config
        self.save_plugin_configs()

    def update_plugin_setting(self, plugin_name: str, key: str, value: Any):
        """Update a specific setting for a plugin."""
        if plugin_name not in self._plugin_configs:
            self._plugin_configs[plugin_name] = PluginConfig(name=plugin_name)
        self._plugin_configs[plugin_name].settings[key] = value
        self.save_plugin_configs()

    def get_load_order(self, plugin_names: List[str]) -> List[str]:
        """
        Get the correct load order for plugins based on dependencies.
        
        Args:
            plugin_names: List of plugin names to order
            
        Returns:
            List of plugin names in correct load order
        """
        # Topological sort based on dependencies
        in_degree = {name: 0 for name in plugin_names}
        graph = {name: [] for name in plugin_names}
        
        # Build dependency graph
        for plugin_name in plugin_names:
            dependencies = self.get_plugin_dependencies(plugin_name)
            for dep in dependencies:
                if dep.name in plugin_names and not dep.optional:
                    graph[dep.name].append(plugin_name)
                    in_degree[plugin_name] += 1
        
        # Topological sort
        queue = [name for name, degree in in_degree.items() if degree == 0]
        result = []
        
        while queue:
            current = queue.pop(0)
            result.append(current)
            
            for neighbor in graph[current]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)
        
        # Add any remaining plugins (circular dependencies or missing deps)
        for plugin_name in plugin_names:
            if plugin_name not in result:
                result.append(plugin_name)
        
        return result
    
    # -------------------------------------------------------------------------
    # Repository Configuration Management
    # -------------------------------------------------------------------------
    
    def get_repository_config(self, repo_name: str) -> Optional[RepositoryConfig]:
        """Get configuration for a repository."""
        return self._repository_configs.get(repo_name)
    
    def set_repository_config(self, config: RepositoryConfig):
        """Set configuration for a repository."""
        self._repository_configs[config.name] = config
        self.save_repository_configs()
    
    def get_all_repository_configs(self) -> Dict[str, RepositoryConfig]:
        """Get all repository configurations."""
        return self._repository_configs.copy()
    
    def remove_repository_config(self, repo_name: str):
        """Remove configuration for a repository."""
        if repo_name in self._repository_configs:
            del self._repository_configs[repo_name]
            self.save_repository_configs()
    
    def save_repository_configs(self):
        """Save repository configurations to file."""
        try:
            data = {name: config.to_dict() for name, config in self._repository_configs.items()}
            with open(self.repositories_config_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            self.logger.error(f"Error saving repository configs: {e}")
    
    # -------------------------------------------------------------------------
    # Global Configuration Management
    # -------------------------------------------------------------------------
    
    def get_global_setting(self, key: str, default: Any = None) -> Any:
        """Get a global setting."""
        return self._global_config.get(key, default)
    
    def set_global_setting(self, key: str, value: Any):
        """Set a global setting."""
        self._global_config[key] = value
        self.save_global_config()
    
    def get_all_global_settings(self) -> Dict[str, Any]:
        """Get all global settings."""
        return self._global_config.copy()
    
    def save_global_config(self):
        """Save global configuration to file."""
        try:
            with open(self.global_config_file, 'w', encoding='utf-8') as f:
                json.dump(self._global_config, f, indent=2, ensure_ascii=False)
        except Exception as e:
            self.logger.error(f"Error saving global config: {e}")
    
    # -------------------------------------------------------------------------
    # Configuration Validation and Migration
    # -------------------------------------------------------------------------
    
    def validate_configurations(self) -> List[str]:
        """Validate all configurations and return list of issues."""
        issues = []
        
        # Validate plugin configs
        for plugin_name, config in self._plugin_configs.items():
            # Check dependencies
            satisfied, missing = self.resolve_dependencies(plugin_name)
            if missing:
                issues.append(f"Plugin '{plugin_name}' has missing dependencies: {', '.join(missing)}")
            
            # Check settings
            if not isinstance(config.settings, dict):
                issues.append(f"Plugin '{plugin_name}' has invalid settings format")
        
        # Validate repository configs
        for repo_name, config in self._repository_configs.items():
            if not config.url:
                issues.append(f"Repository '{repo_name}' has no URL")
            elif not config.url.startswith(('http://', 'https://')):
                issues.append(f"Repository '{repo_name}' has invalid URL: {config.url}")
        
        return issues
    
    def migrate_configuration(self, old_version: str, new_version: str):
        """Migrate configuration from old version to new version."""
        try:
            old_v_parts = parse_version(old_version)
            new_v_parts = parse_version(new_version)
            
            # Add migration logic here based on version changes
            if compare_versions(old_version, new_version) < 0:
                self.logger.info(f"Migrating configuration from {old_version} to {new_version}")
                
                # Example: Add new default settings
                if compare_versions(old_version, "1.1.0") < 0:
                    if "enable_plugin_sandbox" not in self._global_config:
                        self._global_config["enable_plugin_sandbox"] = True
                
                # Save migrated configuration
                self.save_global_config()
                self.save_plugin_configs()
                self.save_repository_configs()
                
                self.logger.info("Configuration migration completed")
                
        except Exception as e:
            self.logger.error(f"Error during configuration migration: {e}")

# End of file
