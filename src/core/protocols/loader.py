"""Loader protocol — interface for molecular file format readers."""
from __future__ import annotations
from typing import Protocol, TYPE_CHECKING
if TYPE_CHECKING:
    from src.core.domain.models.molecule import Molecule

class ILoader(Protocol):
    def supported_extensions(self) -> list[str]: ...
    def can_load(self, path: str) -> bool: ...
    def load(self, path: str, parallel: bool = True) -> Molecule: ...
