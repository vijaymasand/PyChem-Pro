"""
PyDes Molecular Descriptor Calculator

A standalone, custom module for calculating thousands of molecular descriptors.
Operates on PyChem-Pro's native Molecule class.
"""

from .engine import PyDesEngine

__all__ = ["PyDesEngine"]
