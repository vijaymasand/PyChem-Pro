"""Descriptor calculator protocol — interface for molecular property computation."""
from __future__ import annotations
from typing import Protocol, TYPE_CHECKING
if TYPE_CHECKING:
    from src.core.domain.models.molecule import Molecule

class IDescriptorCalculator(Protocol):
    def available_descriptors(self) -> list[str]: ...
    def calculate(self, mol: Molecule, descriptor_names: list[str] | None = None) -> dict[str, float]: ...
    def calculate_batch(self, molecules: list[Molecule], descriptor_names: list[str] | None = None) -> list[dict[str, float]]: ...
