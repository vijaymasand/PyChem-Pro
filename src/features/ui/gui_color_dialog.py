"""
True GUI-based color dialog using tkinter.

Works without console input - actual GUI interface.
"""

import tkinter as tk
from tkinter import colorchooser, messagebox
from typing import Dict, Optional


class GUIColorDialog:
    """
    True GUI-based color dialog using tkinter.
    """
    
    def __init__(self):
        self.root = None
        self.selected_colors = {}
        
    def show_color_dialog(self) -> Dict[str, str]:
        """Show GUI color selection dialog."""
        try:
            # Create main window
            self.root = tk.Tk()
            self.root.title("Color Selection")
            self.root.geometry("400x500")
            self.root.resizable(False, False)
            
            # Create interface
            self._create_interface()
            
            # Center window
            self._center_window()
            
            # Run dialog
            self.root.mainloop()
            
            return self.selected_colors
            
        except Exception as e:
            print(f"GUI dialog error: {e}")
            return {}
    
    def _create_interface(self):
        """Create the GUI interface."""
        # Title
        title_label = tk.Label(self.root, text="Color Selection", font=("Arial", 16, "bold"))
        title_label.pack(pady=10)
        
        # Instructions
        instructions = tk.Label(self.root, text="Select colors for atoms, spheres, and bonds", font=("Arial", 10))
        instructions.pack(pady=5)
        
        # Create notebook for tabs
        self.notebook = tk.Frame(self.root)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Tab buttons
        tab_frame = tk.Frame(self.notebook)
        tab_frame.pack(fill="x")
        
        self.current_tab = tk.StringVar(value="atoms")
        
        tk.Radiobutton(tab_frame, text="Atoms", variable=self.current_tab, 
                      value="atoms", command=self._show_atoms_tab).pack(side="left", padx=10)
        tk.Radiobutton(tab_frame, text="Spheres", variable=self.current_tab, 
                      value="spheres", command=self._show_spheres_tab).pack(side="left", padx=10)
        tk.Radiobutton(tab_frame, text="Bonds", variable=self.current_tab, 
                      value="bonds", command=self._show_bonds_tab).pack(side="left", padx=10)
        
        # Content area
        self.content_frame = tk.Frame(self.notebook)
        self.content_frame.pack(fill="both", expand=True, pady=10)
        
        # Show default tab
        self._show_atoms_tab()
        
        # Separator
        separator = tk.Frame(self.root, height=2, bg="gray")
        separator.pack(fill="x", padx=10, pady=10)
        
        # Status label
        self.status_label = tk.Label(self.root, text="Select colors and click 'Apply Changes'", font=("Arial", 9))
        self.status_label.pack(pady=5)
        
        # Action buttons
        button_frame = tk.Frame(self.root)
        button_frame.pack(fill="x", padx=10, pady=10)
        
        self.apply_btn = tk.Button(button_frame, text="Apply Changes", command=self._apply_colors, 
                                 bg="green", fg="white", font=("Arial", 12, "bold"), width=15)
        self.apply_btn.pack(side="left", padx=5)
        
        self.cancel_btn = tk.Button(button_frame, text="Cancel", command=self._cancel, 
                                  bg="red", fg="white", font=("Arial", 12, "bold"), width=10)
        self.cancel_btn.pack(side="left", padx=5)
        
        self.reset_btn = tk.Button(button_frame, text="Reset All", command=self._reset_colors, 
                                 bg="orange", fg="white", font=("Arial", 12, "bold"), width=10)
        self.reset_btn.pack(side="left", padx=5)
    
    def _show_atoms_tab(self):
        """Show atoms tab."""
        # Clear content frame
        for widget in self.content_frame.winfo_children():
            widget.destroy()
        
        tk.Label(self.content_frame, text="Atom Colors", font=("Arial", 14, "bold")).pack(pady=5)
        
        # Atom color buttons
        atoms = ["C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"]
        
        for i, atom in enumerate(atoms):
            frame = tk.Frame(self.content_frame)
            frame.pack(fill="x", padx=20, pady=2)
            
            tk.Label(frame, text=f"{atom}:", width=5, anchor="w", font=("Arial", 10, "bold")).pack(side="left")
            
            # Current color display
            color_label = tk.Label(frame, text="    ", bg="white", relief="sunken", width=10)
            color_label.pack(side="left", padx=5)
            
            # Color button
            btn = tk.Button(frame, text="Choose Color", 
                          command=lambda a=atom, lbl=color_label: self._choose_atom_color(a, lbl))
            btn.pack(side="left", padx=5)
            
            # Store references
            setattr(self, f"atom_{atom}_label", color_label)
    
    def _show_spheres_tab(self):
        """Show spheres tab."""
        # Clear content frame
        for widget in self.content_frame.winfo_children():
            widget.destroy()
        
        tk.Label(self.content_frame, text="Sphere Colors", font=("Arial", 14, "bold")).pack(pady=5)
        
        # Sphere color buttons
        spheres = ["default", "com", "centroid", "custom"]
        
        for sphere in spheres:
            frame = tk.Frame(self.content_frame)
            frame.pack(fill="x", padx=20, pady=2)
            
            tk.Label(frame, text=f"{sphere.capitalize()}:", width=12, anchor="w", font=("Arial", 10, "bold")).pack(side="left")
            
            # Current color display
            color_label = tk.Label(frame, text="    ", bg="white", relief="sunken", width=10)
            color_label.pack(side="left", padx=5)
            
            # Color button
            btn = tk.Button(frame, text="Choose Color", 
                          command=lambda s=sphere, lbl=color_label: self._choose_sphere_color(s, lbl))
            btn.pack(side="left", padx=5)
            
            # Store references
            setattr(self, f"sphere_{sphere}_label", color_label)
    
    def _show_bonds_tab(self):
        """Show bonds tab."""
        # Clear content frame
        for widget in self.content_frame.winfo_children():
            widget.destroy()
        
        tk.Label(self.content_frame, text="Bond Colors", font=("Arial", 14, "bold")).pack(pady=5)
        
        # Bond color buttons
        bonds = ["default", "single", "double", "triple", "selected", "highlight"]
        
        for bond in bonds:
            frame = tk.Frame(self.content_frame)
            frame.pack(fill="x", padx=20, pady=2)
            
            tk.Label(frame, text=f"{bond.capitalize()}:", width=12, anchor="w", font=("Arial", 10, "bold")).pack(side="left")
            
            # Current color display
            color_label = tk.Label(frame, text="    ", bg="white", relief="sunken", width=10)
            color_label.pack(side="left", padx=5)
            
            # Color button
            btn = tk.Button(frame, text="Choose Color", 
                          command=lambda b=bond, lbl=color_label: self._choose_bond_color(b, lbl))
            btn.pack(side="left", padx=5)
            
            # Store references
            setattr(self, f"bond_{bond}_label", color_label)
    
    def _choose_atom_color(self, atom: str, label):
        """Choose color for atom."""
        color = colorchooser.askcolor(title=f"Choose color for {atom}")
        if color[1]:  # color[1] is the hex color
            label.config(bg=color[1])
            self.selected_colors[f'atom_{atom.lower()}'] = color[1]
    
    def _choose_sphere_color(self, sphere: str, label):
        """Choose color for sphere."""
        color = colorchooser.askcolor(title=f"Choose color for {sphere} sphere")
        if color[1]:
            label.config(bg=color[1])
            self.selected_colors[f'sphere_{sphere}'] = color[1]
    
    def _choose_bond_color(self, bond: str, label):
        """Choose color for bond."""
        color = colorchooser.askcolor(title=f"Choose color for {bond} bonds")
        if color[1]:
            label.config(bg=color[1])
            self.selected_colors[f'stick_{bond}'] = color[1]
    
    def _apply_colors(self):
        """Apply selected colors and close dialog."""
        try:
            if self.selected_colors:
                # Update status
                self.status_label.config(text=f"Applying {len(self.selected_colors)} colors...")
                self.root.update()
                
                # Print debug info
                print(f"DEBUG: Applying colors: {self.selected_colors}")
                
                # Close dialog first
                self.root.quit()
            else:
                # Show warning
                self.status_label.config(text="No colors selected! Please choose colors first.")
                self.root.update()
                
                # Brief delay then reset message
                self.root.after(2000, lambda: self.status_label.config(text="Select colors and click 'Apply Changes'"))
        except Exception as e:
            print(f"ERROR in _apply_colors: {e}")
            self.status_label.config(text=f"Error applying colors: {str(e)}")
            self.root.update()
            # Don't crash - just show error and continue
    
    def _cancel(self):
        """Cancel color selection."""
        self.selected_colors = {}
        self.root.quit()
    
    def _reset_colors(self):
        """Reset all colors to default."""
        self.selected_colors = {}
        
        # Reset all color labels to white
        atoms = ["C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"]
        for atom in atoms:
            label = getattr(self, f"atom_{atom}_label", None)
            if label:
                label.config(bg="white")
        
        spheres = ["default", "com", "centroid", "custom"]
        for sphere in spheres:
            label = getattr(self, f"sphere_{sphere}_label", None)
            if label:
                label.config(bg="white")
        
        bonds = ["default", "single", "double", "triple", "selected", "highlight"]
        for bond in bonds:
            label = getattr(self, f"bond_{bond}_label", None)
            if label:
                label.config(bg="white")
        
        # Update status
        self.status_label.config(text="All colors reset to default")
        self.root.update()
        
        # Brief delay then reset message
        self.root.after(2000, lambda: self.status_label.config(text="Select colors and click 'Apply Changes'"))
    
    def _center_window(self):
        """Center the window on screen."""
        self.root.update_idletasks()
        width = self.root.winfo_width()
        height = self.root.winfo_height()
        x = (self.root.winfo_screenwidth() // 2) - (width // 2)
        y = (self.root.winfo_screenheight() // 2) - (height // 2)
        self.root.geometry(f'{width}x{height}+{x}+{y}')


def show_gui_color_dialog() -> Dict[str, str]:
    """Show GUI color dialog and return selected colors."""
    dialog = GUIColorDialog()
    return dialog.show_color_dialog()


def apply_gui_colors(colors: Dict[str, str]):
    """Apply GUI selected colors to theme."""
    try:
        from src.shared.ui.theme import COLORS
        print(f"DEBUG: Applying colors to theme: {colors}")
        COLORS.update(colors)
        print(f"DEBUG: Successfully applied {len(colors)} GUI colors to theme")
    except Exception as e:
        print(f"ERROR applying colors to theme: {e}")
        import traceback
        traceback.print_exc()
        # Don't crash - just print error
