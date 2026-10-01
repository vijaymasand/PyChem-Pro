"""
Bridge between main application's Domain Molecule and OASA's molecule graph representations.
"""

import types

import src.vendors.oasa.config as oasa_config


def _install_deterministic_ring_order(o_mol):
    """Monkey-patch ``get_smallest_independent_cycles`` on *this instance*
    to return rings in a stable order derived from each vertex's
    ``_stable_idx`` (set by ``domain_to_oasa_mol``). Without this, OASA's
    coordinate generator iterates rings in hash-of-object order — which
    depends on Python object ids and therefore changes every time a
    file is re-read. The resulting layouts can differ not just in
    rotation but in *which ring OASA picks as the backbone*, producing
    genuinely different embeddings across re-opens.
    """
    orig = type(o_mol).get_smallest_independent_cycles

    def stable_rings(self_mol):
        rings = orig(self_mol)
        # Each ring is a set of vertices. Sort each ring's vertices by
        # their stable index, then sort the whole list of rings
        # lexicographically by those tuples. For vertices without a
        # stable index (freshly-added implicit H's, etc.), fall back to
        # symbol and degree — still stable, just weaker.
        def ring_key(r):
            ids = []
            for v in r:
                sid = getattr(v, '_stable_idx', None)
                if sid is None:
                    sid = (v.symbol, len(v.neighbor_edges))
                ids.append(sid)
            return (len(r), tuple(sorted(ids, key=lambda x: (str(type(x).__name__), str(x)))))
        return sorted(rings, key=ring_key)

    o_mol.get_smallest_independent_cycles = types.MethodType(stable_rings, o_mol)


def _assign_stable_idx_fallback(o_mol):
    """Ensure every vertex has a ``_stable_idx`` attribute. Vertices
    created by ``domain_to_oasa_mol`` already have one, but vertices
    added later (by ``add_missing_hydrogens`` or kekulization) might
    not. Assign fallback indices based on traversal order."""
    next_idx = 1_000_000  # comfortably above any domain atom index
    for v in o_mol.vertices:
        if getattr(v, '_stable_idx', None) is None:
            v._stable_idx = next_idx
            next_idx += 1


def oasa_prepare_for_layout(o_mol, infer_bond_orders=False):
    """Run the OASA cleanup pipeline that produces relaxed, well-kekulized
    2D layouts. This is the sequence used by the internal
    ``docking_pose_visualizer`` plugin — which the user has confirmed
    produces clean, non-overlapping 2D views — and is the missing piece
    from both ``CoordinateGenerator2D`` and ``generate_golden_layout``.

    The trick is *kekulization before layout*. OASA's coordinate generator
    places atoms using bond-order cues; if aromatic bonds are left as
    "aromatic" (order 1.5) instead of being localized to alternating
    single/double, the generator cannot correctly orient fused-ring
    systems and produces visibly overlapping layouts.

    Sequence:
        1. ``remove_unimportant_hydrogens()`` — strip implicit H's so
           the heavy-atom skeleton is clean for the aromatic analysis.
        2. (optional) ``add_missing_bond_orders()`` — only when the
           caller has no reliable bond-order info. mol2/SMILES already
           carry bond orders, so calling this is destructive: OASA
           reinterprets bond orders from valences and can pick a
           different Kekulé form each run, producing layouts with
           different implicit-hydrogen counts. Off by default.
        3. ``mark_aromatic_bonds()`` — detect aromatic systems.
        4. ``localize_aromatic_bonds()`` — Kekulize (alternate
           single/double) on top of the caller's original bond orders.
        5. ``add_missing_hydrogens()`` — re-add implicit H's with
           correct geometry now that bond orders are final.

    Every step is wrapped in try/except because OASA occasionally chokes
    on unusual inputs; failure of any step still leaves a usable (if
    less well-kekulized) molecule for the downstream layout pass.
    """
    try:
        o_mol.remove_unimportant_hydrogens()
    except Exception:
        pass
    if infer_bond_orders:
        try:
            o_mol.add_missing_bond_orders()
        except Exception:
            pass
    try:
        o_mol.mark_aromatic_bonds()
    except Exception:
        pass
    try:
        o_mol.localize_aromatic_bonds()
    except Exception:
        pass
    try:
        o_mol.add_missing_hydrogens()
    except Exception:
        pass
    # Assign stable indices to any vertices that don't have one (e.g.
    # hydrogens added by the previous step), then monkey-patch ring
    # enumeration on this instance to return rings in a stable order.
    _assign_stable_idx_fallback(o_mol)
    _install_deterministic_ring_order(o_mol)
    return o_mol

def domain_to_oasa_mol(domain_mol):
    """
    Converts internal Domain Molecule to OASA compatible graph.
    Returns: (oasa_molecule, atom_idx_map)
    where atom_idx_map maps domain atom indices to OASA vertices.
    """
    oasa_mol = oasa_config.Config.create_molecule()
    
    atom_map = {}
    for atom in domain_mol.atoms:
        if atom.symbol == 'H':
            continue

        v = oasa_mol.create_vertex()
        # Stable index for deterministic ring ordering (see
        # _install_deterministic_ring_order). Uses the domain atom
        # index which is file-order-stable and independent of Python
        # object ids.
        v._stable_idx = atom.index
        v.symbol = atom.symbol
        v.charge = atom.formal_charge
        v.isotope = atom.isotope
        v.explicit_hydrogens = atom.num_explicit_h if atom.num_explicit_h is not None else 0
        
        # Merge physical explicit neighbors into OASA's single text-node format 
        for n_idx in domain_mol.get_neighbors(atom.index):
            if domain_mol.atoms[n_idx].symbol == 'H':
                v.explicit_hydrogens += 1
        
        if atom.is_aromatic:
            v.properties_['aromatic'] = 1
        
        # Handle stereo mapping here if needed in SMILES generator
        from src.core.domain.models.atom import Chirality
        if atom.chirality == Chirality.COUNTERCLOCKWISE:
            v.properties_['stereo'] = '@'
        elif atom.chirality == Chirality.CLOCKWISE:
            v.properties_['stereo'] = '@@'
            
        # Copy 2D coordinates if they exist
        if hasattr(atom, 'x2d') and atom.x2d is not None:
            v.x = atom.x2d
            v.y = atom.y2d
            v.z = getattr(atom, 'z2d', 0.0)
            
        oasa_mol.add_vertex(v)
        atom_map[atom.index] = v
        
    for bond in domain_mol.bonds:
        if bond.begin_atom_idx not in atom_map or bond.end_atom_idx not in atom_map:
            continue
            
        v1 = atom_map[bond.begin_atom_idx]
        v2 = atom_map[bond.end_atom_idx]
        e = oasa_mol.create_edge()
        
        from src.core.domain.models.bond import BondType, BondStereo
        if bond.bond_type == BondType.DOUBLE:
            e.order = 2
        elif bond.bond_type == BondType.TRIPLE:
            e.order = 3
        elif bond.bond_type == BondType.AROMATIC:
            e.aromatic = True
            e.order = 1
        else:
            e.order = 1
            
        if bond.stereo == BondStereo.UP:
            e.properties_['stereo'] = '/'
        elif bond.stereo == BondStereo.DOWN:
            e.properties_['stereo'] = '\\'
            
        oasa_mol.add_edge(v1, v2, e=e)
        
    oasa_mol.add_missing_hydrogens()
    return oasa_mol, atom_map
