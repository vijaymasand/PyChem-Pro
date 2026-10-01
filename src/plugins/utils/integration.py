"""
Plugin Integration API

Provides a safe API for plugins to interact with the main application.
This API acts as a bridge between plugins and the main application,
ensuring that plugins cannot directly access sensitive internal components.
"""

from typing import Optional, Callable, Dict, Any, List
import logging
from datetime import datetime


class PluginIntegrationAPI:
    """
    API interface for plugin integration with the main application.
    
    This class provides safe access to main application functionality
    while maintaining security and stability. Plugins interact with
    the application through this API rather than direct access.
    """
    
    def __init__(self, main_window):
        """
        Initialize the plugin API.
        
        Args:
            main_window: The main application window
        """
        self.main_window = main_window
        self.logger = logging.getLogger("plugin.api")
        
        # Track plugin usage for debugging
        self._usage_stats = {
            'api_calls': {},
            'plugin_activity': {}
        }
        
        self.logger.info("Plugin API initialized")
    
    # -------------------------------------------------------------------------
    # Molecule Access Methods
    # -------------------------------------------------------------------------
    
    def get_current_molecule(self):
        """
        Get the current molecule from the main application.
        
        Returns:
            Current molecule object or None if no molecule loaded
        """
        try:
            molecule = getattr(self.main_window, 'molecule', None)
            self._track_api_call('get_current_molecule')
            return molecule
        except Exception as e:
            self.logger.error(f"Error getting current molecule: {e}")
            return None
    
    def get_molecule_info(self) -> Dict[str, Any]:
        """
        Get information about the current molecule.
        
        Returns:
            Dictionary with molecule information
        """
        try:
            molecule = self.get_current_molecule()
            if not molecule:
                return {'has_molecule': False}
            
            info = {
                'has_molecule': True,
                'atom_count': len(molecule.atoms) if hasattr(molecule, 'atoms') else 0,
                'bond_count': len(molecule.bonds) if hasattr(molecule, 'bonds') else 0,
                'has_coordinates': hasattr(molecule.atoms[0], 'x') if molecule.atoms else False,
                'formula': getattr(molecule, 'formula', 'Unknown')
            }
            
            self._track_api_call('get_molecule_info')
            return info
            
        except Exception as e:
            self.logger.error(f"Error getting molecule info: {e}")
            return {'has_molecule': False, 'error': str(e)}
    
    def get_atoms_data(self) -> List[Dict[str, Any]]:
        """
        Get atom data from the current molecule.
        
        Returns:
            List of dictionaries with atom information
        """
        try:
            molecule = self.get_current_molecule()
            if not molecule or not hasattr(molecule, 'atoms'):
                return []
            
            atoms_data = []
            for i, atom in enumerate(molecule.atoms):
                atom_data = {
                    'index': i,
                    'symbol': getattr(atom, 'symbol', 'Unknown'),
                    'atomic_number': getattr(atom, 'atomic_number', 0)
                }
                
                # Add coordinates if available
                if hasattr(atom, 'x'):
                    atom_data.update({
                        'x': atom.x,
                        'y': atom.y,
                        'z': atom.z
                    })
                
                atoms_data.append(atom_data)
            
            self._track_api_call('get_atoms_data')
            return atoms_data
            
        except Exception as e:
            self.logger.error(f"Error getting atoms data: {e}")
            return []
    
    def get_bonds_data(self) -> List[Dict[str, Any]]:
        """
        Get bond data from the current molecule.
        
        Returns:
            List of dictionaries with bond information
        """
        try:
            molecule = self.get_current_molecule()
            if not molecule or not hasattr(molecule, 'bonds'):
                return []
            
            bonds_data = []
            for i, bond in enumerate(molecule.bonds):
                bond_data = {
                    'index': i,
                    'atom1_index': getattr(bond, 'atom1_index', 0),
                    'atom2_index': getattr(bond, 'atom2_index', 0),
                    'bond_type': str(getattr(bond, 'bond_type', 'Unknown'))
                }
                bonds_data.append(bond_data)
            
            self._track_api_call('get_bonds_data')
            return bonds_data
            
        except Exception as e:
            self.logger.error(f"Error getting bonds data: {e}")
            return []
    
    # -------------------------------------------------------------------------
    # Viewer Access Methods
    # -------------------------------------------------------------------------
    
    def get_viewer_2d(self):
        """
        Get the 2D molecular viewer.
        
        Returns:
            2D viewer widget or None
        """
        try:
            viewer = getattr(self.main_window, 'viewer_2d', None)
            self._track_api_call('get_viewer_2d')
            return viewer
        except Exception as e:
            self.logger.error(f"Error getting 2D viewer: {e}")
            return None
    
    def get_viewer_3d(self):
        """
        Get the 3D molecular viewer.
        
        Returns:
            3D viewer widget or None
        """
        try:
            viewer = getattr(self.main_window, 'viewer_3d', None)
            self._track_api_call('get_viewer_3d')
            return viewer
        except Exception as e:
            self.logger.error(f"Error getting 3D viewer: {e}")
            return None
    
    def update_viewer_2d(self):
        """Update the 2D molecular viewer."""
        try:
            viewer = self.get_viewer_2d()
            if viewer and hasattr(viewer, 'update'):
                viewer.update()
                self._track_api_call('update_viewer_2d')
            else:
                self.logger.warning("2D viewer not available for update")
        except Exception as e:
            self.logger.error(f"Error updating 2D viewer: {e}")
    
    def update_viewer_3d(self):
        """Update the 3D molecular viewer."""
        try:
            viewer = self.get_viewer_3d()
            if viewer and hasattr(viewer, 'update'):
                viewer.update()
                self._track_api_call('update_viewer_3d')
            else:
                self.logger.warning("3D viewer not available for update")
        except Exception as e:
            self.logger.error(f"Error updating 3D viewer: {e}")
    
    # -------------------------------------------------------------------------
    # UI Integration Methods
    # -------------------------------------------------------------------------
    
    def show_status_message(self, message: str, timeout: int = 3000):
        """
        Show a message in the main application status bar.
        
        Args:
            message: Message to display
            timeout: Message timeout in milliseconds
        """
        try:
            status_bar = getattr(self.main_window, 'status_bar', None)
            if status_bar and hasattr(status_bar, 'showMessage'):
                status_bar.showMessage(message, timeout)
                self._track_api_call('show_status_message')
            else:
                self.logger.info(f"Status: {message}")
        except Exception as e:
            self.logger.error(f"Error showing status message: {e}")
    
    def show_error_message(self, title: str, message: str):
        """
        Show an error message dialog.
        
        Args:
            title: Dialog title
            message: Error message
        """
        try:
            from src.shared.qt_compat import QMessageBox
            QMessageBox.critical(self.main_window, title, message)
            self._track_api_call('show_error_message')
        except Exception as e:
            self.logger.error(f"Error showing error message: {e}")
    
    def show_info_message(self, title: str, message: str):
        """
        Show an information message dialog.
        
        Args:
            title: Dialog title
            message: Information message
        """
        try:
            from src.shared.qt_compat import QMessageBox
            QMessageBox.information(self.main_window, title, message)
            self._track_api_call('show_info_message')
        except Exception as e:
            self.logger.error(f"Error showing info message: {e}")
    
    def show_warning_message(self, title: str, message: str):
        """
        Show a warning message dialog.
        
        Args:
            title: Dialog title
            message: Warning message
        """
        try:
            from src.shared.qt_compat import QMessageBox
            QMessageBox.warning(self.main_window, title, message)
            self._track_api_call('show_warning_message')
        except Exception as e:
            self.logger.error(f"Error showing warning message: {e}")
    
    def add_menu_item(self, menu_path: str, text: str, callback: Callable, shortcut: str = None):
        """
        Add a menu item to the main application.
        
        Args:
            menu_path: Path to menu (e.g., "Tools.My Plugin")
            text: Menu item text
            callback: Function to call when item is clicked
            shortcut: Keyboard shortcut (optional)
        """
        try:
            # This would need to be implemented in the main window
            # For now, just log the request
            self.logger.info(f"Menu item requested: {menu_path}.{text}")
            self._track_api_call('add_menu_item')
            
            # TODO: Implement actual menu item addition
            # This would require modifying the main window to support dynamic menu items
            
        except Exception as e:
            self.logger.error(f"Error adding menu item: {e}")
    
    # -------------------------------------------------------------------------
    # File System Methods (Safe Access)
    # -------------------------------------------------------------------------
    
    def get_project_directory(self) -> Optional[str]:
        """
        Get the current project directory.
        
        Returns:
            Project directory path or None
        """
        try:
            # Get the current working directory or project directory
            import os
            cwd = os.getcwd()
            self._track_api_call('get_project_directory')
            return cwd
        except Exception as e:
            self.logger.error(f"Error getting project directory: {e}")
            return None
    
    def show_file_dialog(self, title: str, file_filter: str = "All Files (*)", 
                        save_mode: bool = False) -> Optional[str]:
        """
        Show a file dialog for opening or saving files.
        
        Args:
            title: Dialog title
            file_filter: File filter string
            save_mode: True for save dialog, False for open dialog
            
        Returns:
            Selected file path or None
        """
        try:
            from src.shared.qt_compat import QFileDialog
            
            if save_mode:
                filepath, _ = QFileDialog.getSaveFileName(
                    self.main_window, title, "", file_filter
                )
            else:
                filepath, _ = QFileDialog.getOpenFileName(
                    self.main_window, title, "", file_filter
                )
            
            self._track_api_call('show_file_dialog')
            return filepath if filepath else None
            
        except Exception as e:
            self.logger.error(f"Error showing file dialog: {e}")
            return None
    
    # -------------------------------------------------------------------------
    # Plugin Management Methods
    # -------------------------------------------------------------------------
    
    def get_plugin_settings(self, plugin_name: str) -> Dict[str, Any]:
        """
        Get settings for a specific plugin.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            Plugin settings dictionary
        """
        try:
            # This would integrate with a settings system
            # For now, return empty settings
            self._track_api_call('get_plugin_settings')
            return {}
        except Exception as e:
            self.logger.error(f"Error getting plugin settings: {e}")
            return {}
    
    def set_plugin_settings(self, plugin_name: str, settings: Dict[str, Any]):
        """
        Set settings for a specific plugin.
        
        Args:
            plugin_name: Name of the plugin
            settings: Settings dictionary
        """
        try:
            # This would integrate with a settings system
            # For now, just log the request
            self.logger.info(f"Settings saved for plugin: {plugin_name}")
            self._track_api_call('set_plugin_settings')
        except Exception as e:
            self.logger.error(f"Error setting plugin settings: {e}")
    
    # -------------------------------------------------------------------------
    # Utility Methods
    # -------------------------------------------------------------------------
    
    def _track_api_call(self, method_name: str):
        """Track API call for debugging and statistics."""
        if method_name not in self._usage_stats['api_calls']:
            self._usage_stats['api_calls'][method_name] = 0
        self._usage_stats['api_calls'][method_name] += 1
    
    def get_usage_stats(self) -> Dict[str, Any]:
        """Get API usage statistics."""
        return self._usage_stats.copy()
    
    def clear_usage_stats(self):
        """Clear API usage statistics."""
        self._usage_stats = {'api_calls': {}, 'plugin_activity': {}}
