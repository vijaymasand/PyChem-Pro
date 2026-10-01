"""Force field protocol — interface for geometry optimization services."""
from __future__ import annotations
from typing import Protocol, TYPE_CHECKING
from dataclasses import dataclass, field
if TYPE_CHECKING:
    from src.core.domain.models.molecule import Molecule

@dataclass
class OptimizationResult:
    converged: bool
    final_energy: float
    num_steps: int
    energy_trajectory: list[float] = field(default_factory=list)
    rms_gradient: float = 0.0

class IForceField(Protocol):
    def add_hydrogens(self, mol: Molecule) -> int: ...
    def assign_atom_types(self, mol: Molecule) -> None: ...
    def assign_charges(self, mol: Molecule) -> None: ...
    def optimize_geometry(self, mol: Molecule, max_iters: int = 500,
        convergence: float = 1e-4, method: str = 'lbfgs') -> OptimizationResult: ...
    def compute_energy(self, mol: Molecule) -> float: ...
