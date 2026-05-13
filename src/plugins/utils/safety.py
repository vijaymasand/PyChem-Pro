"""
Plugin Safety Mechanisms

Provides safety and sandboxing mechanisms for plugins to ensure
they cannot harm the main application or system.
"""

import sys
import os
import logging
import importlib.util
from typing import Dict, Any, List, Optional
from pathlib import Path


class PluginSandbox:
    """
    Provides a safe environment for plugin execution.
    
    This class implements various safety mechanisms to ensure
    that plugins cannot access sensitive system resources or
    harm the main application.
    """
    
    def __init__(self):
        """Initialize the plugin sandbox."""
        self.logger = logging.getLogger("plugin.sandbox")
        
        # Restricted modules that plugins cannot import
        self._restricted_modules = {
            'os.system', 'subprocess', 'shutil',
            'socket', 'urllib.request', 'http.client',
            'ftplib', 'smtplib', 'telnetlib', 'pickle',
            'marshal', 'code', 'compile', 'eval', 'exec'
        }
        
        # Restricted built-in functions
        self._restricted_builtins = {
            'open', 'file', 'input', 'raw_input', 'eval',
            'exec', 'compile', '__import__', 'reload'
        }
        
        # Allowed file extensions for plugin files
        self._allowed_extensions = {'.py'}
        
        # Maximum plugin file size (1MB)
        self._max_file_size = 1024 * 1024
        
        self.logger.info("Plugin sandbox initialized")
    
    def validate_plugin_file(self, file_path: str) -> tuple[bool, str]:
        """
        Validate a plugin file for safety.
        
        Args:
            file_path: Path to the plugin file
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            path = Path(file_path)
            
            # Check if file exists
            if not path.exists():
                return False, "Plugin file does not exist"
            
            # Check file extension
            if path.suffix not in self._allowed_extensions:
                return False, f"Invalid file extension: {path.suffix}"
            
            # Check file size
            if path.stat().st_size > self._max_file_size:
                return False, f"Plugin file too large: {path.stat().st_size} bytes"
            
            # Read file content
            with open(path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Check for restricted imports
            for restricted in self._restricted_modules:
                if restricted in content:
                    return False, f"Restricted module found: {restricted}"
            
            # Check for restricted built-ins
            for restricted in self._restricted_builtins:
                if restricted in content and not self._is_safe_usage(content, restricted):
                    return False, f"Restricted function found: {restricted}"
            
            return True, "Plugin file is safe"
            
        except Exception as e:
            return False, f"Error from .validation import PluginCodeValidator: {e}"
    
    def _is_safe_usage(self, content: str, restricted_func: str) -> bool:
        """
        Check if a restricted function is used safely.
        
        Args:
            content: File content
            restricted_func: Restricted function name
            
        Returns:
            True if usage is safe, False otherwise
        """
        # Allow 'open' for file I/O operations (needed for I/O plugins and export functionality)
        if restricted_func == 'open':
            # Check for common safe patterns:
            # - 'with open(' for proper file handling
            # - 'open(filename, 'w', ...)' or 'open(filename, 'r', ...)'
            import re
            
            # Pattern: with open(...) as ...:
            safe_patterns = [
                r'with\s+open\s*\([^)]+\)\s+as\s+\w+',
                r'open\s*\([^)]+,[\s]*["\']w',
                r'open\s*\([^)]+,[\s]*["\']r',
                r'open\s*\(\s*filename',
                r'open\s*\(\s*filepath',
                r'open\s*\(\s*file_path',
                r'open\s*\(\s*output',
                r'open\s*\(\s*dest_path',
            ]
            
            for pattern in safe_patterns:
                if re.search(pattern, content):
                    return True
            
            # If open is used in a safe context (inside export/import functions)
            if 'export' in content.lower() or 'import' in content.lower() or 'save' in content.lower():
                return True
                
        # Allow 'file' when used in file path contexts (not the old Python 2 file type)
        if restricted_func == 'file':
            # Only allow 'file' as part of variable names like 'filename', 'filepath', etc.
            # not as a standalone function call like 'file(...)'
            import re
            if not re.search(r'\bfile\s*\(', content):
                return True
                
        # For other restricted functions, be strict
        return False
    
    def safe_import(self, module_name: str, file_path: str = None):
        """
        Safely import a plugin module.
        
        Args:
            module_name: Name of the module to import
            file_path: Path to the plugin file
            
        Returns:
            Imported module or None if import failed
        """
        try:
            if file_path:
                # Import from specific file
                spec = importlib.util.spec_from_file_location(module_name, file_path)
                if spec is None:
                    self.logger.error(f"Could not load spec for {module_name}")
                    return None
                
                module = importlib.util.module_from_spec(spec)
                
                # Add safety attributes to the module
                self._add_safety_attributes(module)
                
                spec.loader.exec_module(module)
                return module
            else:
                # Import from package
                self._add_safety_attributes(sys.modules[module_name])
                return sys.modules[module_name]
                
        except Exception as e:
            self.logger.error(f"Error importing plugin module {module_name}: {e}")
            return None
    
    def _add_safety_attributes(self, module):
        """Add safety attributes to a module."""
        # This could be used to add monitoring or restrictions
        # For now, it's a placeholder for future enhancements
        pass
    
    def check_dependencies(self, dependencies: List[str]) -> tuple[bool, str]:
        """
        Check if plugin dependencies are available and safe.
        
        Args:
            dependencies: List of required dependencies
            
        Returns:
            Tuple of (all_available, missing_or_unsafe)
        """
        missing = []
        unsafe = []
        
        for dep in dependencies:
            try:
                # Check if dependency is available
                __import__(dep)
                
                # Check if dependency is safe (basic check)
                if dep in self._restricted_modules:
                    unsafe.append(dep)
                    
            except ImportError:
                missing.append(dep)
        
        if unsafe:
            return False, f"Unsafe dependencies: {unsafe}"
        elif missing:
            return False, f"Missing dependencies: {missing}"
        else:
            return True, "All dependencies available and safe"
    
    def get_safe_environment(self) -> Dict[str, Any]:
        """
        Get a safe environment for plugin execution.
        
        Returns:
            Safe environment dictionary
        """
        # Start with basic safe built-ins
        safe_env = {
            '__builtins__': {
                'abs': abs,
                'all': all,
                'any': any,
                'bool': bool,
                'dict': dict,
                'enumerate': enumerate,
                'float': float,
                'int': int,
                'len': len,
                'list': list,
                'max': max,
                'min': min,
                'pow': pow,
                'range': range,
                'reversed': reversed,
                'round': round,
                'sorted': sorted,
                'str': str,
                'sum': sum,
                'tuple': tuple,
                'type': type,
                'zip': zip,
            }
        }
        
        # Add safe modules
        safe_modules = {
            'math': __import__('math'),
            'random': __import__('random'),
            'datetime': __import__('datetime'),
            'json': __import__('json'),
            'csv': __import__('csv'),
            're': __import__('re'),
            'collections': __import__('collections'),
            'itertools': __import__('itertools'),
            'functools': __import__('functools'),
            'operator': __import__('operator'),
        }
        
        safe_env.update(safe_modules)
        
        return safe_env


class PluginValidator:
    """
    Validates plugins for compliance and safety.
    """
    
    def __init__(self):
        """Initialize the plugin validator."""
        self.logger = logging.getLogger("plugin.validator")
        self.sandbox = PluginSandbox()
    
    def validate_plugin(self, plugin_class, file_path: str) -> tuple[bool, List[str]]:
        """
        Validate a plugin class.
        
        Args:
            plugin_class: Plugin class to validate
            file_path: Path to the plugin file
            
        Returns:
            Tuple of (is_valid, list_of_issues)
        """
        issues = []
        
        try:
            # Check if plugin inherits from BasePlugin
            # Using string comparison of class names to avoid issues with different import paths
            if not any(base.__name__ == 'BasePlugin' for base in plugin_class.__mro__):
                issues.append("Plugin must inherit from BasePlugin")
            
            # Check required methods
            required_methods = ['create_widget']
            for method in required_methods:
                if not hasattr(plugin_class, method):
                    issues.append(f"Missing required method: {method}")
            
            # Check plugin info
            if hasattr(plugin_class, '__init__'):
                try:
                    # Try to create an instance to check info
                    instance = plugin_class()
                    if not hasattr(instance, 'info'):
                        issues.append("Plugin must have info attribute")
                    else:
                        info = instance.info
                        if not hasattr(info, 'name') or not info.name:
                            issues.append("Plugin info must have a name")
                        if not hasattr(info, 'version') or not info.version:
                            issues.append("Plugin info must have a version")
                        if not hasattr(info, 'plugin_type'):
                            issues.append("Plugin info must have a plugin_type")
                except Exception as e:
                    issues.append(f"Error creating plugin instance: {e}")
            
            # Validate file safety
            is_safe, safety_msg = self.sandbox.validate_plugin_file(file_path)
            if not is_safe:
                issues.append(f"Safety issue: {safety_msg}")
            
            return len(issues) == 0, issues
            
        except Exception as e:
            return False, [f"Validation error: {e}"]
    
    def validate_plugin_dependencies(self, dependencies: List[str]) -> tuple[bool, str]:
        """
        Validate plugin dependencies.
        
        Args:
            dependencies: List of required dependencies
            
        Returns:
            Tuple of (are_valid, message)
        """
        return self.sandbox.check_dependencies(dependencies)
