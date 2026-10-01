"""Event bus and event types for decoupled communication."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Any, TYPE_CHECKING
from collections import defaultdict
if TYPE_CHECKING:
    from src.core.domain.models.molecule import Molecule
    from src.core.protocols.forcefield import OptimizationResult

@dataclass
class MoleculeChangedEvent:
    molecule: Any
    source: str

@dataclass
class OptimizationCompleteEvent:
    molecule: Any
    result: Any

@dataclass
class RenderRequestEvent:
    molecule: Any
    mode: str

class EventBus:
    def __init__(self):
        self._subscribers: dict[type, list[Callable]] = defaultdict(list)

    def subscribe(self, event_type: type, handler: Callable):
        self._subscribers[event_type].append(handler)

    def unsubscribe(self, event_type: type, handler: Callable):
        self._subscribers[event_type].remove(handler)

    def publish(self, event):
        for handler in self._subscribers[type(event)]:
            handler(event)
