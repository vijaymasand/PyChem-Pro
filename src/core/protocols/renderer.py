"""Renderer protocol — interface for molecular visualization backends."""
from __future__ import annotations
from typing import Protocol, Any, TYPE_CHECKING
from dataclasses import dataclass
if TYPE_CHECKING:
    from src.core.domain.models.molecule import Molecule

@dataclass
class Camera:
    rot_x: float; rot_y: float; pan_x: float; pan_y: float
    zoom: float; width: int; height: int

class RenderMode:
    BALL_AND_STICK = 'ball_and_stick'; SPACEFILL = 'spacefill'
    WIREFRAME = 'wireframe'; CARTOON = 'cartoon'
    RIBBON = 'ribbon'; BACKBONE = 'backbone'

class IRenderer(Protocol):
    @property
    def name(self) -> str: ...
    def supports_molecule(self, mol: Molecule) -> bool: ...
    def render(self, mol: Molecule, camera: Camera, mode: str, target: Any) -> None: ...
    def pick_atom(self, mol: Molecule, camera: Camera, screen_x: int, screen_y: int) -> int | None: ...
    def export_image(self, mol: Molecule, camera: Camera, mode: str, width: int, height: int) -> Any: ...
