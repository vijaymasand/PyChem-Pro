"""
Molecule Scene — multi-object container for PyMOL-style overlay.

PyChem historically held a single ``Molecule``.  ``MoleculeScene`` lets the
application keep several loaded structures (any mix of SDF / MOL / MOL2 / PDB /
PDBx) at once.  Each structure is wrapped in a :class:`MoleculeObject` carrying
its own display state (visibility, representation, colour).

The 3D viewer is fed a single *merged* display molecule built by
:meth:`MoleculeScene.build_overlay` so the whole existing render pipeline is
reused unchanged: per-object representation maps onto the renderer's
``custom_atom_modes`` override and per-object colour onto ``custom_atom_colors``
(see ``painter_renderer.py``).  Visibility simply includes/excludes an object
from the merge.

The *active* object (``active_id``) is the one driving the 2D view, chemistry
actions and single-structure export — ``window.molecule`` always points at it,
preserving every existing single-molecule code path.
"""

import copy
import string

from src.shared.qt_compat import QObject, Signal
from src.core.domain.models.molecule import Molecule


# Representation strings the renderer honours *per atom* via custom_atom_modes.
# Cartoon / ribbon / backbone are whole-chain modes and stay global.
PER_ATOM_REPRESENTATIONS = ('ball_and_stick', 'spacefill', 'wireframe')

# Pool of single-character chain IDs used to give each merged object its own
# chain namespace (see build_overlay).  Falls back to multi-char IDs if a scene
# somehow exhausts these 62 — harmless, since the merged molecule is display-only.
_CHAIN_ID_POOL = string.ascii_uppercase + string.ascii_lowercase + string.digits


def _looks_like_protein(mol) -> bool:
    """True if the molecule has a protein backbone.

    Prefers the loader's ``is_protein`` flag, but falls back to detecting
    α-carbons directly (``pdb_name == 'CA'`` on a carbon atom).  The parallel
    PDB loader used for large files doesn't set the flag, so relying on it alone
    would make big proteins render as ball-and-stick instead of a cartoon.
    """
    if mol.properties.get('is_protein'):
        return True
    count = 0
    for a in mol.atoms:
        if (getattr(a, 'pdb_name', '') or '').strip() == 'CA' and a.symbol == 'C':
            count += 1
            if count >= 3:
                return True
    return False


def _alloc_chain_id(used: set) -> str:
    """Return a chain ID not present in *used* (without adding it)."""
    for ch in _CHAIN_ID_POOL:
        if ch not in used:
            return ch
    i = 1
    while f"C{i}" in used:
        i += 1
    return f"C{i}"

# PyMOL "color by chain" (cbc) inspired palette for the auto-colour option.
AUTO_COLOR_PALETTE = [
    (0x33, 0xCC, 0x33),  # green
    (0x33, 0xCC, 0xCC),  # cyan
    (0xCC, 0x66, 0xCC),  # light magenta
    (0xCC, 0xCC, 0x33),  # yellow
    (0xCC, 0x66, 0x66),  # salmon
    (0x66, 0x99, 0xCC),  # slate
    (0xCC, 0x99, 0x33),  # orange
    (0x99, 0xCC, 0x66),  # lime
    (0x33, 0x99, 0x99),  # deep teal
    (0xCC, 0x33, 0x99),  # hot pink
    (0x99, 0x66, 0xCC),  # violet
    (0x66, 0xCC, 0xCC),  # marine
]


class MoleculeObject:
    """One loaded structure plus its per-object display state."""

    __slots__ = (
        'id', 'name', 'molecule', 'visible',
        'representation', 'color', 'source_path', 'source_format',
        '_palette_idx',
    )

    def __init__(self, obj_id, name, molecule, source_path=None,
                 source_format='unknown'):
        self.id = obj_id
        self.name = name
        self.molecule = molecule
        self.visible = True
        # None  -> use the global render mode
        # one of PER_ATOM_REPRESENTATIONS -> applied to just this object
        self.representation = None
        # None      -> element (CPK) colouring
        # 'auto'    -> a distinct palette colour
        # (r, g, b) -> a flat custom colour
        self.color = None
        self.source_path = source_path
        self.source_format = source_format
        self._palette_idx = None  # assigned the first time color == 'auto'

    @property
    def num_atoms(self):
        return len(self.molecule.atoms) if self.molecule else 0


class SceneOverlay:
    """Result of merging all visible objects into one display molecule."""

    __slots__ = (
        'molecule', 'atom_origin', 'object_ranges',
        'custom_modes', 'custom_colors', 'any_protein',
    )

    def __init__(self, molecule, atom_origin, object_ranges,
                 custom_modes, custom_colors, any_protein):
        self.molecule = molecule              # merged Molecule (or None if empty)
        self.atom_origin = atom_origin        # merged_idx -> (obj_id, local_idx)
        self.object_ranges = object_ranges    # obj_id -> (start, end)  end exclusive
        self.custom_modes = custom_modes      # merged_idx -> representation
        self.custom_colors = custom_colors    # merged_idx -> (r, g, b)
        self.any_protein = any_protein


class MoleculeScene(QObject):
    """Ordered collection of :class:`MoleculeObject` with an active object."""

    # Emitted whenever the set of objects, their order, or display state changes.
    scene_changed = Signal()
    # Emitted when the active object changes (carries the new active id, or -1).
    active_changed = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.objects = []          # list[MoleculeObject]
        self.active_id = -1
        self._id_counter = 0
        self._color_counter = 0

    # ── Lookup ───────────────────────────────────────────────────────

    def get(self, obj_id):
        for obj in self.objects:
            if obj.id == obj_id:
                return obj
        return None

    def active_object(self):
        return self.get(self.active_id)

    def is_empty(self):
        return not self.objects

    def visible_objects(self):
        return [o for o in self.objects if o.visible]

    # ── Mutation ─────────────────────────────────────────────────────

    def add_object(self, molecule, name, source_path=None,
                   source_format='unknown', make_active=True):
        """Add a structure as a new object. Returns its id."""
        obj_id = self._id_counter
        self._id_counter += 1
        unique = self._unique_name(name or f"object_{obj_id}")
        obj = MoleculeObject(obj_id, unique, molecule, source_path, source_format)
        self.objects.append(obj)
        if make_active or self.active_id == -1:
            self.active_id = obj_id
            self.active_changed.emit(obj_id)
        self.scene_changed.emit()
        return obj_id

    def remove_object(self, obj_id):
        obj = self.get(obj_id)
        if obj is None:
            return
        self.objects.remove(obj)
        if self.active_id == obj_id:
            self.active_id = self.objects[-1].id if self.objects else -1
            self.active_changed.emit(self.active_id)
        self.scene_changed.emit()

    def clear(self):
        self.objects = []
        self.active_id = -1
        self.active_changed.emit(-1)
        self.scene_changed.emit()

    def set_visible(self, obj_id, visible):
        obj = self.get(obj_id)
        if obj is not None and obj.visible != visible:
            obj.visible = visible
            self.scene_changed.emit()

    def set_all_visible(self, visible):
        changed = False
        for obj in self.objects:
            if obj.visible != visible:
                obj.visible = visible
                changed = True
        if changed:
            self.scene_changed.emit()

    def set_active(self, obj_id):
        if obj_id != self.active_id and self.get(obj_id) is not None:
            self.active_id = obj_id
            self.active_changed.emit(obj_id)
            self.scene_changed.emit()

    def set_representation(self, obj_id, representation):
        obj = self.get(obj_id)
        if obj is not None:
            obj.representation = representation
            self.scene_changed.emit()

    def set_color(self, obj_id, color):
        """color: None (CPK) | 'auto' | (r, g, b)."""
        obj = self.get(obj_id)
        if obj is None:
            return
        if color == 'auto' and obj._palette_idx is None:
            obj._palette_idx = self._color_counter
            self._color_counter += 1
        obj.color = color
        self.scene_changed.emit()

    # ── Overlay construction ─────────────────────────────────────────

    def build_overlay(self):
        """Merge all *visible* objects into one display molecule."""
        merged = Molecule(name="__scene__")
        merged.begin_bulk_load()
        atom_origin = {}
        object_ranges = {}
        custom_modes = {}
        custom_colors = {}
        any_protein = False
        # Chain IDs already claimed by earlier objects.  Each object gets its
        # own chain namespace so two structures that share a PDB chain letter
        # (e.g. both use chain 'A') are never fused by chain-based rendering.
        used_chain_ids = set()

        # Merge the active object FIRST so its local atom indices equal the
        # merged indices — this keeps the active molecule (window.molecule)
        # index-compatible with the 2D viewer and every per-atom feature.
        ordered = list(self.objects)
        active = self.active_object()
        if active is not None and active.visible and active in ordered:
            ordered.remove(active)
            ordered.insert(0, active)

        for obj in ordered:
            if not obj.visible or obj.molecule is None or not obj.molecule.atoms:
                continue
            mol = obj.molecule
            start = len(merged.atoms)
            offset = start
            # Atoms from a non-protein object (mol2/sdf/small-molecule PDB) are
            # marked as het-atoms so the renderers treat them as ligands — shown
            # over a protein cartoon, and kept out of protein-only code paths
            # (the GL packer assumes protein atoms carry an integer res_seq).
            obj_is_protein = _looks_like_protein(mol)

            # Remap this object's chain IDs onto a fresh, collision-free set so
            # the cartoon tracer (which groups residues by chain_id) keeps each
            # structure as its own chain instead of splicing them together.
            chain_remap = {}

            for local_idx, atom in enumerate(mol.atoms):
                na = copy.copy(atom)
                if not obj_is_protein:
                    na.is_hetatm = True
                orig_chain = getattr(atom, 'chain_id', None)
                key = orig_chain if orig_chain is not None else ''
                new_chain = chain_remap.get(key)
                if new_chain is None:
                    # Keep the original letter when it's still free (typical for
                    # the first/active object); otherwise allocate a fresh one.
                    if key and key not in used_chain_ids:
                        new_chain = key
                    else:
                        new_chain = _alloc_chain_id(used_chain_ids)
                    used_chain_ids.add(new_chain)
                    chain_remap[key] = new_chain
                na.chain_id = new_chain
                merged_idx = merged.add_atom(na)
                atom_origin[merged_idx] = (obj.id, local_idx)

            for bond in mol.bonds:
                nb = copy.copy(bond)
                nb.begin_atom_idx = bond.begin_atom_idx + offset
                nb.end_atom_idx = bond.end_atom_idx + offset
                nb.index = len(merged.bonds)
                merged.bonds.append(nb)
                merged._adjacency[nb.begin_atom_idx].append((nb.end_atom_idx, nb.index))
                merged._adjacency[nb.end_atom_idx].append((nb.begin_atom_idx, nb.index))

            end = len(merged.atoms)
            object_ranges[obj.id] = (start, end)

            if obj.representation in PER_ATOM_REPRESENTATIONS:
                for midx in range(start, end):
                    custom_modes[midx] = obj.representation

            rgb = self._resolve_color(obj)
            if rgb is not None:
                for midx in range(start, end):
                    custom_colors[midx] = rgb

            if obj_is_protein:
                any_protein = True

        merged.end_bulk_load()
        if any_protein:
            merged.properties['is_protein'] = True

        merged_dummy_spheres = []
        for obj in ordered:
            if not obj.visible or obj.molecule is None or not obj.molecule.atoms:
                continue
            if hasattr(obj.molecule, 'dummy_spheres') and obj.molecule.dummy_spheres:
                merged_dummy_spheres.extend(obj.molecule.dummy_spheres)
        merged.dummy_spheres = merged_dummy_spheres

        return SceneOverlay(
            merged if merged.atoms else None,
            atom_origin, object_ranges, custom_modes, custom_colors, any_protein,
        )

    # ── Helpers ──────────────────────────────────────────────────────

    def _resolve_color(self, obj):
        if obj.color is None:
            return None
        if obj.color == 'auto':
            idx = obj._palette_idx if obj._palette_idx is not None else 0
            return AUTO_COLOR_PALETTE[idx % len(AUTO_COLOR_PALETTE)]
        return tuple(obj.color)

    def _unique_name(self, name):
        existing = {o.name for o in self.objects}
        if name not in existing:
            return name
        i = 2
        while f"{name}_{i}" in existing:
            i += 1
        return f"{name}_{i}"
