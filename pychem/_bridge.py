"""Internal: lazy singleton for ServiceRegistry."""
from src.core.registry import ServiceRegistry

_registry = None


def get_registry() -> ServiceRegistry:
    """Get or create the global ServiceRegistry singleton."""
    global _registry
    if _registry is None:
        _registry = ServiceRegistry()
    return _registry
