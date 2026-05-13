"""
MMFF94 atom typing. Public-ish: plugin authors and Phase 3 will replace
the default typer via this module's exports.
"""
from src.services.forcefield.atom_typing.base import AtomTyper
from src.services.forcefield.atom_typing.hand_coded import HandCodedAtomTyper


def default_typer() -> AtomTyper:
    """Return the currently-recommended atom typer."""
    return HandCodedAtomTyper()


__all__ = ["AtomTyper", "HandCodedAtomTyper", "default_typer"]
