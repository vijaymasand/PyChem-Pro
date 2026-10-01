"""Coordinate generator protocol — interface for 3D structure generation."""
from __future__ import annotations
from typing import Protocol, TYPE_CHECKING
if TYPE_CHECKING:
    from src.core.domain.models.molecule import Molecule

class ICoordinateGenerator(Protocol):
    def generate_3d(self, mol: Molecule, optimize: bool = True, max_steps: int = 200) -> None: ...
    def generate_conformers(self, mol: Molecule, n: int = 10) -> list[Molecule]: ...
