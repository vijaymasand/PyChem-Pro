"""
Descriptor Types and Categories

Defines the structure and categories for molecular descriptors.
"""

from enum import Enum
from typing import Dict, List, Any, Callable
from dataclasses import dataclass

class DescriptorCategory(Enum):
    """Descriptor calculation categories."""
    CONSTITUTIONAL = "Constitutional"
    TOPOLOGICAL = "Topological"
    GEOMETRIC = "Geometric"
    ELECTRONIC = "Electronic"
    QUANTUM = "Quantum"
    FINGERPRINTS = "Fingerprints"
    HYBRID = "Hybrid"
    CUSTOM = "Custom"

@dataclass
class DescriptorInfo:
    """Information about a descriptor."""
    name: str
    description: str
    category: DescriptorCategory
    formula: str = ""
    unit: str = ""
    range_min: float = None
    range_max: float = None
    calculation_function: Callable = None
    
@dataclass
class DescriptorResult:
    """Result of descriptor calculation."""
    descriptor_name: str
    value: Any
    molecule_id: str = ""
    selection_id: str = ""
    calculation_time: float = 0.0

@dataclass
class CalculationProgress:
    """Progress tracking for descriptor calculations."""
    current_descriptor: str = ""
    current_category: str = ""
    completed: int = 0
    total: int = 0
    percentage: float = 0.0
    message: str = ""

class SelectionType(Enum):
    """Types of atom selections."""
    ALL = "all"
    CUSTOM = "custom"
    FRAGMENT = "fragment"
    ENVIRONMENT = "environment"
    RING = "ring"
    CHAIN = "chain"
    AROMATIC = "aromatic"
    ALIPHATIC = "aliphatic"

@dataclass
class AtomSelection:
    """Atom selection for descriptor calculation."""
    selection_type: SelectionType
    atom_indices: List[int]
    description: str = ""
    parameters: Dict[str, Any] = None
