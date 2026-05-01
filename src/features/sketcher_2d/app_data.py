# -*- coding: utf-8 -*-
import os
from src.vendors.oasa.periodic_table import periodic_table
from src.shared.qt_compat import QIcon, QPixmap, QStandardPaths

class App:
    """ Stores application wide data """
    window = None
    paper = None # selected current paper
    tool = None # selected current tool
    clipboard = [] # List of cloned objects
    template_manager = None # created only once
    SRC_DIR = os.path.dirname(__file__)
    # We use a custom data dir for PyChem sketcher
    DATA_DIR = QStandardPaths.writableLocation(QStandardPaths.AppDataLocation) + "/PyChem/Sketcher"
    if not os.path.exists(DATA_DIR):
        try:
            os.makedirs(DATA_DIR)
        except:
            pass
    dark_mode = True # Default to dark mode for PyChem

class Default:
    """ default values for drawing properties. """
    atom_font_size = 12
    bond_length = 40
    bond_width = 1.2
    bond_spacing = 6.0
    electron_dot_size = 2.0
    arrow_line_width = 2.0
    arrow_head_dimensions = (5, 2, 1.5)
    plus_size = 18

class Settings:
    """ settings for some properties """
    basic_scale = 1.0
    render_dpi = 100
    atom_font_name = "Sans Serif"
    coord_head_dimensions = 6, 2.5, 2
    min_arrow_length = 30
    text_size = 14
    focus_color = (0, 255, 0)
    selection_color = (150, 150, 255)
    image_export_dpi = 100
    image_export_margin = 10

# initialize Settings with Default values
for key, val in dict(vars(Default)).items():
    if not key.startswith("__"):
        setattr(Settings, key, val)

# Colors
basic_colors = [
    "#000000", "#404040", "#6b6b6b", "#808080", "#909090", "#ffffff",
    "#790874", "#f209f1", "#09007c", "#000def", "#047f7d", "#05fef8",
    "#7e0107", "#f00211", "#fff90d", "#07e00d", "#067820", "#827d05",
]

def atomic_num_to_symbol(atomic_num):
    for symbol, data in periodic_table.items():
        if data.get('ord') == atomic_num:
            return symbol
    raise ValueError("Invalid atomic number %s" % str(atomic_num))

auto_hydrogen_elements = {"B", "C", "Si", "N", "P", "As", "O", "S", "F", "Cl", "Br", "I"}
