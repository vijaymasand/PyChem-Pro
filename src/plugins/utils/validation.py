"""
Plugin Validation Utilities

Provides validation functions for plugins to ensure they meet
requirements and are safe to use.
"""

import ast
import re
from typing import List, Dict, Any, Optional
from pathlib import Path


class PluginCodeValidator:
    """
    Validates plugin code structure and content.
    """
    
    def __init__(self):
        """Initialize the code validator."""
        self.required_methods = ['create_widget']
        self.recommended_methods = [
            'initialize', 'cleanup', 'on_molecule_changed',
            'on_plugin_activated', 'on_plugin_deactivated'
        ]
    
    def validate_plugin_class(self, file_path: str) -> Dict[str, Any]:
        """
        Validate plugin class structure.
        
        Args:
            file_path: Path to the plugin file
            
        Returns:
            Validation result dictionary
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Parse AST
            tree = ast.parse(content)
            
            # Find plugin classes
            plugin_classes = []
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    # Check if it looks like a plugin class
                    if self._is_plugin_class(node, content):
                        plugin_classes.append(node)
            
            if not plugin_classes:
                return {
                    'valid': False,
                    'errors': ['No valid plugin class found'],
                    'warnings': [],
                    'suggestions': []
                }
            
            # Validate each plugin class
            results = []
            for plugin_class in plugin_classes:
                result = self._validate_class_node(plugin_class, content)
                results.append(result)
            
            # Return the best result (first valid class)
            for result in results:
                if result['valid']:
                    return result
            
            # If no valid classes, return the first result
            return results[0] if results else {
                'valid': False,
                'errors': ['No plugin classes found'],
                'warnings': [],
                'suggestions': []
            }
            
        except Exception as e:
            return {
                'valid': False,
                'errors': [f'Code parsing error: {e}'],
                'warnings': [],
                'suggestions': []
            }
    
    def _is_plugin_class(self, class_node: ast.ClassDef, content: str) -> bool:
        """Check if a class node represents a plugin class."""
        # Check if it inherits from BasePlugin
        if class_node.bases:
            for base in class_node.bases:
                if isinstance(base, ast.Name) and base.id == 'BasePlugin':
                    return True
                elif isinstance(base, ast.Attribute):
                    if hasattr(base, 'attr') and base.attr == 'BasePlugin':
                        return True
        
        # Check if it has plugin-like methods
        for node in class_node.body:
            if isinstance(node, ast.FunctionDef):
                if node.name in self.required_methods:
                    return True
        
        return False
    
    def _validate_class_node(self, class_node: ast.ClassDef, content: str) -> Dict[str, Any]:
        """Validate a specific plugin class."""
        errors = []
        warnings = []
        suggestions = []
        
        # Check inheritance
        has_base_plugin = False
        for base in class_node.bases:
            if isinstance(base, ast.Name) and base.id == 'BasePlugin':
                has_base_plugin = True
                break
        
        if not has_base_plugin:
            errors.append("Plugin class must inherit from BasePlugin")
        
        # Check methods
        found_methods = []
        for node in class_node.body:
            if isinstance(node, ast.FunctionDef):
                found_methods.append(node.name)
        
        # Check required methods
        for method in self.required_methods:
            if method not in found_methods:
                errors.append(f"Missing required method: {method}")
        
        # Check recommended methods
        for method in self.recommended_methods:
            if method not in found_methods:
                suggestions.append(f"Consider implementing: {method}")
        
        # Check __init__ method
        if '__init__' in found_methods:
            init_node = next(node for node in class_node.body 
                           if isinstance(node, ast.FunctionDef) and node.name == '__init__')
            
            # Check if __init__ calls super().__init__
            calls_super = False
            for node in ast.walk(init_node):
                if isinstance(node, ast.Call):
                    if isinstance(node.func, ast.Attribute):
                        if node.func.attr == '__init__':
                            calls_super = True
                            break
            
            if not calls_super:
                warnings.append("__init__ should call super().__init__()")
        
        # Check docstring
        if not ast.get_docstring(class_node):
            suggestions.append("Add a docstring to the plugin class")
        
        return {
            'valid': len(errors) == 0,
            'errors': errors,
            'warnings': warnings,
            'suggestions': suggestions,
            'class_name': class_node.name,
            'found_methods': found_methods
        }


class PluginMetadataValidator:
    """
    Validates plugin metadata and information.
    """
    
    def __init__(self):
        """Initialize the metadata validator."""
        pass
    
    def validate_plugin_info(self, info) -> Dict[str, Any]:
        """
        Validate plugin information.
        
        Args:
            info: PluginInfo object
            
        Returns:
            Validation result dictionary
        """
        errors = []
        warnings = []
        suggestions = []
        
        # Check required fields
        required_fields = ['name', 'version', 'description', 'author', 'plugin_type']
        for field in required_fields:
            if not hasattr(info, field) or not getattr(info, field):
                errors.append(f"Missing required field: {field}")
        
        # Check name format
        if hasattr(info, 'name') and info.name:
            if not re.match(r'^[a-zA-Z][a-zA-Z0-9_ ]*$', info.name):
                warnings.append("Plugin name should contain only letters, numbers, spaces, and underscores")
            
            if len(info.name) < 3:
                warnings.append("Plugin name should be at least 3 characters long")
            
            if len(info.name) > 50:
                warnings.append("Plugin name should be less than 50 characters")
        
        # Check version format
        if hasattr(info, 'version') and info.version:
            if not re.match(r'^\d+\.\d+\.\d+', info.version):
                suggestions.append("Version should follow semantic versioning (e.g., 1.0.0)")
        
        # Check description
        if hasattr(info, 'description') and info.description:
            if len(info.description) < 10:
                warnings.append("Description should be at least 10 characters long")
            
            if len(info.description) > 500:
                suggestions.append("Consider keeping description under 500 characters")
        
        # Check dependencies
        if hasattr(info, 'dependencies') and info.dependencies:
            for dep in info.dependencies:
                if not re.match(r'^[a-zA-Z][a-zA-Z0-9_-]*$', dep):
                    warnings.append(f"Invalid dependency name format: {dep}")
        
        return {
            'valid': len(errors) == 0,
            'errors': errors,
            'warnings': warnings,
            'suggestions': suggestions
        }


def validate_plugin_file(file_path: str) -> Dict[str, Any]:
    """
    Comprehensive validation of a plugin file.
    
    Args:
        file_path: Path to the plugin file
        
    Returns:
        Complete validation result
    """
    code_validator = PluginCodeValidator()
    metadata_validator = PluginMetadataValidator()
    
    # Validate code structure
    code_result = code_validator.validate_plugin_class(file_path)
    
    # Try to extract and validate metadata
    metadata_result = {'valid': True, 'errors': [], 'warnings': [], 'suggestions': []}
    
    try:
        # Import the plugin to get metadata
        import importlib.util
        spec = importlib.util.spec_from_file_location("temp_plugin", file_path)
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
        
        if plugin_classes:
            # Create instance to get info - handle both old and new style plugins
            try:
                # Try creating without arguments first (new style)
                try:
                    instance = plugin_classes[0]()
                except TypeError:
                    # Try creating with a dummy info (old style that requires info)
                    from src.plugins.plugin_types import PluginInfo, PluginType
                    dummy_info = PluginInfo(
                        name="Temp",
                        version="1.0.0",
                        description="Temp",
                        author="Temp",
                        plugin_type=PluginType.ANALYSIS
                    )
                    instance = plugin_classes[0](dummy_info)
                
                if hasattr(instance, 'get_info'):
                    metadata_result = metadata_validator.validate_plugin_info(instance.get_info())
            except Exception as e:
                metadata_result['errors'].append(f"Error getting plugin metadata: {e}")
    
    except Exception as e:
        metadata_result['errors'].append(f"Error importing plugin for metadata validation: {e}")
    
    # Combine results
    all_errors = code_result['errors'] + metadata_result['errors']
    all_warnings = code_result['warnings'] + metadata_result['warnings']
    all_suggestions = code_result['suggestions'] + metadata_result['suggestions']
    
    return {
        'valid': len(all_errors) == 0,
        'errors': all_errors,
        'warnings': all_warnings,
        'suggestions': all_suggestions,
        'code_validation': code_result,
        'metadata_validation': metadata_result
    }
