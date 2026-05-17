"""
Pharmacophore Features Plugin for PyChem

Ported from liquid_1.0_for python3.py (PyMOL) to PyChem's architecture.
Computes pharmacophore features (Lipophilic, H-Bond Donor, H-Bond Acceptor,
Aromatic, Positive/Negative Ionizable), clusters them using Local Feature
Density, computes PCA-oriented ellipsoids via NIPALS, and generates
120-dimensional Gaussian pharmacophore descriptors.

Algorithm References:
    - liquid_1.0: Dr. Vijay Masand (PyMOL script, 2020)
    - NIPALS PCA: Wold, 1966
    - Pharmacophore typing: Catalyst/Phase conventions

Features:
    - 7 pharmacophore types (3 original + 4 new)
    - Local Feature Density clustering with configurable radii
    - NIPALS PCA for cluster shape/orientation
    - 120-D Gaussian descriptor export (.liquid format)
    - Interactive Qt visualization with ellipses and legend
    - JSON / CSV / clipboard export

Dependencies:
    - numpy (already required by PyChem)
    - PySide6 (already required by PyChem)

Author:
    Dr. Vijay Masand (original algorithm), PyChem Team (port & enhancements)

Version:
    1.0.0

License:
    MIT
"""

import math
import json
import os
import logging
from typing import Dict, List, Tuple, Optional, Any, Set

import numpy as np

# =============================================================================
# Component 1: Pure-Python Pharmacophore Engine (no Qt)
# =============================================================================

# ── Pharmacophore type constants ─────────────────────────────────────────────

PHARM_LIPOPHILIC = "Lipophilic"
PHARM_DONOR = "H-Bond Donor"
PHARM_ACCEPTOR = "H-Bond Acceptor"
PHARM_DONOR_ACCEPTOR = "Donor+Acceptor"
PHARM_AROMATIC = "Aromatic"
PHARM_POSITIVE = "Positive Ionizable"
PHARM_NEGATIVE = "Negative Ionizable"

#: Canonical order for descriptor pair enumeration.
PHARM_TYPES_ORDERED = [
    PHARM_LIPOPHILIC, PHARM_DONOR, PHARM_ACCEPTOR,
    PHARM_AROMATIC, PHARM_POSITIVE, PHARM_NEGATIVE,
]

#: Color scheme for pharmacophore features (hex).
PHARMACOPHORE_COLORS: Dict[str, str] = {
    PHARM_LIPOPHILIC:     "#2E7D32",   # Green
    PHARM_DONOR:          "#1565C0",   # Blue
    PHARM_ACCEPTOR:       "#C62828",   # Red
    PHARM_DONOR_ACCEPTOR: "#9C27B0",   # Magenta/Purple
    PHARM_AROMATIC:       "#FF8F00",   # Orange
    PHARM_POSITIVE:       "#00ACC1",   # Cyan
    PHARM_NEGATIVE:       "#E91E63",   # Pink
}

#: Default cluster radii (Å) per type, matching liquid_1.0.
DEFAULT_RADII: Dict[str, float] = {
    PHARM_LIPOPHILIC: 2.0,
    PHARM_DONOR:      2.0,
    PHARM_ACCEPTOR:   2.0,
    PHARM_AROMATIC:   2.5,
    PHARM_POSITIVE:   2.5,
    PHARM_NEGATIVE:   2.5,
}


# ── Atom Typing ──────────────────────────────────────────────────────────────

def classify_pharmacophore_features(molecule) -> Dict[str, List[int]]:
    """
    Classify each atom into pharmacophore feature types.

    Ported from liquid_1.0 lines 22-24 with extensions for aromatic
    and ionizable types.

    Args:
        molecule: PyChem Molecule object.

    Returns:
        Dict mapping pharmacophore type name to list of atom indices.
    """
    features: Dict[str, List[int]] = {t: [] for t in PHARMACOPHORE_COLORS}

    # Pre-compute ring membership for aromatic detection
    try:
        rings = molecule.find_rings()
    except Exception:
        rings = []

    aromatic_atoms: Set[int] = set()
    for ring in rings:
        ring_aromatic = all(molecule.atoms[idx].is_aromatic for idx in ring)
        if ring_aromatic:
            aromatic_atoms.update(ring)

    for atom in molecule.atoms:
        if not atom.has_coords:
            continue
        sym = atom.symbol
        neighbors = molecule.get_neighbors(atom.index)
        neighbor_symbols = {molecule.atoms[n].symbol for n in neighbors}

        is_donor = False
        is_acceptor = False

        # ── Lipophilic ───────────────────────────────────────────
        # liquid: elem C not neighbor N/O, or S not neighbor H/N/O,
        #         or halogens (Cl, I, Br, F)
        if sym == 'C' and not (neighbor_symbols & {'N', 'O'}):
            features[PHARM_LIPOPHILIC].append(atom.index)
        elif sym == 'S' and not (neighbor_symbols & {'H', 'N', 'O'}):
            features[PHARM_LIPOPHILIC].append(atom.index)
        elif sym in ('Cl', 'Br', 'I', 'F'):
            features[PHARM_LIPOPHILIC].append(atom.index)

        # ── H-Bond Donor ─────────────────────────────────────────
        # liquid: elem N,O neighbor hydro
        if sym in ('N', 'O') and 'H' in neighbor_symbols:
            is_donor = True

        # ── H-Bond Acceptor ──────────────────────────────────────
        # liquid: elem O, or N not neighbor H and formal_charge == 0
        if sym == 'O':
            is_acceptor = True
        elif sym == 'N' and 'H' not in neighbor_symbols and atom.formal_charge == 0:
            is_acceptor = True

        # ── Classify donor/acceptor/both ─────────────────────────
        if is_donor and is_acceptor:
            features[PHARM_DONOR_ACCEPTOR].append(atom.index)
        elif is_donor:
            features[PHARM_DONOR].append(atom.index)
        elif is_acceptor:
            features[PHARM_ACCEPTOR].append(atom.index)

        # ── Aromatic ─────────────────────────────────────────────
        if atom.index in aromatic_atoms:
            features[PHARM_AROMATIC].append(atom.index)

        # ── Positive Ionizable ───────────────────────────────────
        if sym == 'N' and atom.formal_charge > 0:
            features[PHARM_POSITIVE].append(atom.index)

        # ── Negative Ionizable ───────────────────────────────────
        if sym in ('O', 'S') and atom.formal_charge < 0:
            features[PHARM_NEGATIVE].append(atom.index)

    return features


# ── Distance Helpers ─────────────────────────────────────────────────────────

def _get_coords_array(molecule, indices: List[int]) -> np.ndarray:
    """Extract Nx3 coordinate matrix for given atom indices."""
    coords = []
    for idx in indices:
        a = molecule.atoms[idx]
        coords.append([a.x, a.y, a.z])
    return np.array(coords, dtype=np.float64)


def _pairwise_distances(coords: np.ndarray) -> np.ndarray:
    """Compute NxN pairwise Euclidean distance matrix."""
    # Efficient broadcasting: d(i,j) = ||c_i - c_j||
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    return np.sqrt(np.sum(diff ** 2, axis=2))


# ── Local Feature Density ────────────────────────────────────────────────────

def calculate_lfds(molecule, indices: List[int], rc: float) -> np.ndarray:
    """
    Calculate Local Feature Densities for a set of atoms.

    Ported from liquid_1.0 calculateLFDs (lines 50-61).

    LFD(i) = Σ_j max(0, 1 - d(i,j)/rc)

    Args:
        molecule: PyChem Molecule object.
        indices: Atom indices of the same pharmacophore type.
        rc: Cluster radius (Å).

    Returns:
        1D array of LFD values, one per atom in `indices`.
    """
    if len(indices) == 0:
        return np.array([], dtype=np.float64)
    if len(indices) == 1:
        return np.array([1.0], dtype=np.float64)

    coords = _get_coords_array(molecule, indices)
    dists = _pairwise_distances(coords)

    # LFD: sum of max(0, 1 - d/rc) over all j
    contributions = np.maximum(0.0, 1.0 - dists / rc)
    lfds = np.sum(contributions, axis=1)
    return lfds


# ── Clustering ───────────────────────────────────────────────────────────────

def cluster_features(
    molecule, indices: List[int], lfds: np.ndarray, rc: float
) -> Tuple[List[int], List[List[int]]]:
    """
    Cluster features by merging each atom into the highest-LFD neighbor
    within radius rc. Union-find approach.

    Ported from liquid_1.0 cluster() (lines 66-97).

    Args:
        molecule: PyChem Molecule object.
        indices: Atom indices of the same pharmacophore type.
        lfds: LFD array (from calculate_lfds).
        rc: Cluster radius (Å).

    Returns:
        (find, union) where:
            find[i] = index of cluster leader for atom i
            union[i] = list of member indices for cluster i (empty if merged away)
    """
    n = len(indices)
    if n == 0:
        return [], []

    coords = _get_coords_array(molecule, indices)
    dists = _pairwise_distances(coords)

    # Initialize union-find
    union = [[i] for i in range(n)]
    find = list(range(n))

    for k in range(n):
        for g in range(n):
            if dists[k, g] <= rc:
                if lfds[find[k]] <= lfds[find[g]]:
                    if find[k] != find[g]:
                        old_root = find[k]
                        new_root = find[g]
                        # Merge old into new
                        union[new_root] = union[new_root] + union[old_root]
                        union[old_root] = []
                        # Update find pointers
                        for u in union[new_root]:
                            find[u] = new_root

    return find, union


def count_clusters(find: List[int]) -> int:
    """Count the number of non-empty clusters."""
    return sum(1 for i, f in enumerate(find) if f == i)


# ── Geometric Centers ────────────────────────────────────────────────────────

def calculate_geometric_centers(
    molecule, indices: List[int], union: List[List[int]],
    lfds: Optional[np.ndarray] = None, weighted: bool = False
) -> List[Optional[np.ndarray]]:
    """
    Compute geometric centers for each cluster.

    Ported from liquid_1.0 calculateGC() (lines 147-229).

    Args:
        molecule: PyChem Molecule object.
        indices: Atom indices of the same pharmacophore type.
        union: Cluster membership lists from cluster_features.
        lfds: LFD values (required if weighted=True).
        weighted: If True, use LFD-weighted centroids.

    Returns:
        List of 3D center arrays (None for empty clusters).
    """
    if not indices:
        return []

    all_coords = _get_coords_array(molecule, indices)
    centers: List[Optional[np.ndarray]] = []

    for i, members in enumerate(union):
        if not members:
            centers.append(None)
            continue

        member_coords = all_coords[members]

        if not weighted or lfds is None:
            center = np.mean(member_coords, axis=0)
        else:
            member_lfds = lfds[members]
            total_lfd = np.sum(member_lfds)
            if total_lfd > 0:
                weights = member_lfds / total_lfd
                center = np.sum(member_coords * weights[:, np.newaxis], axis=0)
            else:
                center = np.mean(member_coords, axis=0)

        centers.append(center)

    return centers


# ── PCA via NIPALS ───────────────────────────────────────────────────────────

def compute_pca_nipals(
    molecule, indices: List[int], union: List[List[int]],
    centers: List[Optional[np.ndarray]], default_sigma: float = 0.75
) -> List[Optional[List[np.ndarray]]]:
    """
    Find principal components for each cluster via NIPALS algorithm.

    Ported from liquid_1.0 PCA() (lines 239-443).
    Returns 3 eigenvectors per cluster scaled by their standard deviation.

    Special handling:
        - Singletons: 3 default-length axis-aligned vectors
        - 2 atoms: 1 vector along the atom pair, 2 via cross product
        - 3 atoms: 2 via NIPALS, 1 via cross product
        - ≥4 atoms: full 3 NIPALS PCs

    Args:
        molecule: PyChem Molecule object.
        indices: Atom indices of the pharmacophore type.
        union: Cluster membership lists.
        centers: Geometric centers from calculate_geometric_centers.
        default_sigma: Default axis length for singletons/degenerate cases.

    Returns:
        List of [pc1, pc2, pc3] arrays per cluster (None for empty clusters).
        Each pc is a 3D vector whose length = standard deviation along that axis.
    """
    if not indices:
        return []

    all_coords = _get_coords_array(molecule, indices)
    pcs_list: List[Optional[List[np.ndarray]]] = []

    for i, members in enumerate(union):
        if not members:
            pcs_list.append(None)
            continue

        center = centers[i]
        if center is None:
            pcs_list.append(None)
            continue

        n_members = len(members)

        if n_members == 1:
            # Singleton: default axes
            pcs = [
                np.array([default_sigma, 0.0, 0.0]),
                np.array([0.0, default_sigma, 0.0]),
                np.array([0.0, 0.0, default_sigma]),
            ]
            pcs_list.append(pcs)
            continue

        # Center the data
        X = all_coords[members] - center  # (n_members, 3)

        if n_members == 2:
            # 2 atoms: single direction vector
            diff = X[0] - X[1]
            p1 = diff / 2.0
            length = np.linalg.norm(p1)
            if length < default_sigma:
                p1 = p1 / (length + 1e-12) * default_sigma

            # Find orthogonal vectors
            p2 = _find_orthogonal(p1, default_sigma)
            p3 = np.cross(p1, p2)
            p3_len = np.linalg.norm(p3)
            if p3_len > 0:
                p3 = p3 / p3_len * default_sigma
            else:
                p3 = np.array([0., 0., default_sigma])

            pcs_list.append([p1, p2, p3])
            continue

        # General case: n_members >= 3
        # Use eigendecomposition of the covariance matrix for robustness.
        # This is equivalent to the NIPALS approach in liquid_1.0 but handles
        # degenerate/planar cases (e.g. benzene) without numerical issues.
        cov = X.T @ X  # 3x3 (unnormalized covariance)

        eigvals, eigvecs = np.linalg.eigh(cov)
        # eigh returns eigenvalues in ascending order; reverse for descending
        sort_idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[sort_idx]
        eigvecs = eigvecs[:, sort_idx]

        pcs = []
        for pc_idx in range(3):
            p = eigvecs[:, pc_idx]
            # Standard deviation along this PC (matching liquid_1.0 formula)
            projections = X @ p
            std_val = np.sqrt(
                np.sum(projections ** 2) / max(1, len(X) - 1)
            )
            if std_val < default_sigma:
                std_val = default_sigma
            pcs.append(p * std_val)

        pcs_list.append(pcs)


    return pcs_list


def _find_orthogonal(v: np.ndarray, length: float) -> np.ndarray:
    """Find a vector orthogonal to v with given length."""
    v_norm = np.linalg.norm(v)
    if v_norm < 1e-12:
        return np.array([0., length, 0.])

    vn = v / v_norm
    # Choose the axis least parallel to v
    if abs(vn[0]) <= abs(vn[1]) and abs(vn[0]) <= abs(vn[2]):
        candidate = np.array([1., 0., 0.])
    elif abs(vn[1]) <= abs(vn[2]):
        candidate = np.array([0., 1., 0.])
    else:
        candidate = np.array([0., 0., 1.])

    ortho = np.cross(vn, candidate)
    ortho_len = np.linalg.norm(ortho)
    if ortho_len > 0:
        return ortho / ortho_len * length
    return np.array([0., length, 0.])


def _safe_cross(a: np.ndarray, b: np.ndarray, length: float) -> np.ndarray:
    """Cross product scaled to given length, with fallback."""
    c = np.cross(a, b)
    c_len = np.linalg.norm(c)
    if c_len > 1e-12:
        return c / c_len * length
    return np.array([0., 0., length])


# ── Gaussian Pharmacophore Descriptor ────────────────────────────────────────

def _univariate_gaussian(x: float, sigma: float) -> float:
    """Univariate Gaussian density. Ported from liquid uniGD()."""
    if sigma < 1e-12:
        return 0.0
    return (1.0 / math.sqrt(2 * math.pi * sigma ** 2)) * math.exp(
        -0.5 * (x ** 2 / sigma ** 2)
    )


def _trivariate_gaussian(point: np.ndarray, sigmas: np.ndarray, w: float) -> float:
    """
    Trivariate factored Gaussian. Ported from liquid trivGD().
    Uses absolute coordinates as the original does.
    """
    ax, ay, az = abs(point[0]), abs(point[1]), abs(point[2])
    return w * (
        _univariate_gaussian(ax, sigmas[0])
        * _univariate_gaussian(ay, sigmas[1])
        * _univariate_gaussian(az, sigmas[2])
    )


def compute_pharmacophore_descriptors(
    centers_by_type: Dict[str, List[Optional[np.ndarray]]],
    pcs_by_type: Dict[str, List[Optional[List[np.ndarray]]]],
    n_bins: int = 20
) -> Dict[str, np.ndarray]:
    """
    Compute pairwise Gaussian pharmacophore descriptors.

    Ported from liquid_1.0 writeDescriptorfile() / pair() (lines 880-1135).

    For the 3 original types (L, D, A), this produces 6 pair types × 20 bins
    = 120-dimensional descriptor, matching the original .liquid format.

    Args:
        centers_by_type: Dict of pharmacophore type -> list of cluster centers.
        pcs_by_type: Dict of pharmacophore type -> list of PCA axes per cluster.
        n_bins: Number of bins per pair type (default 20).

    Returns:
        Dict mapping pair name (e.g. "L-L", "L-D") to 1D array of bin values.
    """
    # Build compact lists of non-empty centers and their PCs
    def _compact(type_key):
        ctrs = centers_by_type.get(type_key, [])
        pcas = pcs_by_type.get(type_key, [])
        valid = []
        for i, c in enumerate(ctrs):
            if c is not None and i < len(pcas) and pcas[i] is not None:
                valid.append((c, pcas[i]))
        return valid

    type_keys = [PHARM_LIPOPHILIC, PHARM_DONOR, PHARM_ACCEPTOR]
    type_data = {k: _compact(k) for k in type_keys}

    # Short names for pair labeling
    short = {PHARM_LIPOPHILIC: "L", PHARM_DONOR: "D", PHARM_ACCEPTOR: "A"}

    results = {}
    pair_types = [
        (PHARM_LIPOPHILIC, PHARM_LIPOPHILIC),
        (PHARM_LIPOPHILIC, PHARM_DONOR),
        (PHARM_LIPOPHILIC, PHARM_ACCEPTOR),
        (PHARM_DONOR, PHARM_DONOR),
        (PHARM_DONOR, PHARM_ACCEPTOR),
        (PHARM_ACCEPTOR, PHARM_ACCEPTOR),
    ]

    for t1, t2 in pair_types:
        pair_name = f"{short[t1]}-{short[t2]}"
        data1 = type_data[t1]
        data2 = type_data[t2]
        bins = np.zeros(n_bins, dtype=np.float64)

        if not data1 or not data2:
            results[pair_name] = bins
            continue

        # Count pairs
        if t1 == t2:
            n_pairs = len(data1) * (len(data1) - 1) // 2
        else:
            n_pairs = len(data1) * len(data2)

        if n_pairs == 0:
            results[pair_name] = bins
            continue

        for a_idx, (ca, pca_a) in enumerate(data1):
            for b_idx, (cb, pca_b) in enumerate(data2):
                if t1 == t2 and b_idx <= a_idx:
                    continue
                if np.allclose(ca, cb):
                    continue

                # Sigma vectors (lengths of PCA axes)
                sigma_a = np.array([np.linalg.norm(pc) for pc in pca_a])
                sigma_b = np.array([np.linalg.norm(pc) for pc in pca_b])

                # Ensure minimum sigma
                sigma_a = np.maximum(sigma_a, 0.1)
                sigma_b = np.maximum(sigma_b, 0.1)

                # Rotation matrices (axes normalized)
                rot_a = np.array([pc / (np.linalg.norm(pc) + 1e-12) for pc in pca_a])
                rot_b = np.array([pc / (np.linalg.norm(pc) + 1e-12) for pc in pca_b])

                vec_ab = cb - ca
                vec_ba = ca - cb

                for x in range(n_bins):
                    dist = x + 0.5

                    # Point along the line from a to b
                    vec_ab_norm = np.linalg.norm(vec_ab)
                    if vec_ab_norm < 1e-12:
                        continue

                    point_vl = vec_ab / vec_ab_norm * dist
                    point_vl_world = point_vl + ca

                    point_vr = vec_ba / vec_ab_norm * dist
                    point_vr_world = point_vr + cb

                    # Transform into local frames
                    xb = rot_b @ point_vl_world
                    xa = rot_a @ point_vr_world

                    w = 1.0
                    prb_a = _trivariate_gaussian(xa, sigma_a, w)
                    prb_b = _trivariate_gaussian(xb, sigma_b, w)

                    result = (1.0 / max(1, n_pairs)) * prb_a * prb_b
                    bins[x] += result

        # Normalize bins (surface to 1)
        total = np.sum(bins)
        if total > 0:
            bins = bins / total

        results[pair_name] = bins

    return results


# ── Summary Descriptors ──────────────────────────────────────────────────────

def compute_feature_summary(
    features: Dict[str, List[int]],
    cluster_counts: Dict[str, int],
    centers_by_type: Dict[str, List[Optional[np.ndarray]]],
) -> Dict[str, Any]:
    """
    Compute a human-readable summary of pharmacophore features.

    Returns:
        Dict with feature counts, cluster counts, and inter-center
        distance statistics.
    """
    summary = {
        "feature_counts": {t: len(v) for t, v in features.items()},
        "cluster_counts": dict(cluster_counts),
        "total_features": sum(len(v) for v in features.values()),
        "total_clusters": sum(cluster_counts.values()),
    }

    # Inter-center distances for each type
    dist_stats = {}
    for t, ctrs in centers_by_type.items():
        valid = [c for c in ctrs if c is not None]
        if len(valid) >= 2:
            coords = np.array(valid)
            dists = _pairwise_distances(coords)
            np.fill_diagonal(dists, np.inf)
            dist_stats[t] = {
                "min_dist": float(np.min(dists)),
                "max_dist": float(np.max(dists[dists < np.inf])) if np.any(dists < np.inf) else 0.0,
                "mean_dist": float(np.mean(dists[dists < np.inf])) if np.any(dists < np.inf) else 0.0,
            }
    summary["inter_center_distances"] = dist_stats

    return summary


# ── Full Pipeline ────────────────────────────────────────────────────────────

def run_pharmacophore_pipeline(
    molecule, radii: Optional[Dict[str, float]] = None
) -> Dict[str, Any]:
    """
    Run the complete pharmacophore analysis pipeline.

    Args:
        molecule: PyChem Molecule object with 3D coordinates.
        radii: Optional dict overriding default cluster radii per type.

    Returns:
        Dict containing all results:
            features, lfds, clusters, centers, pcs, descriptors, summary
    """
    if radii is None:
        radii = dict(DEFAULT_RADII)

    # Step 1: Atom typing
    features = classify_pharmacophore_features(molecule)

    # Steps 2-5: Per-type clustering and PCA
    lfds_by_type = {}
    clusters_by_type = {}
    centers_by_type = {}
    pcs_by_type = {}
    cluster_counts = {}

    for pharm_type, atom_indices in features.items():
        if not atom_indices:
            lfds_by_type[pharm_type] = np.array([])
            clusters_by_type[pharm_type] = ([], [])
            centers_by_type[pharm_type] = []
            pcs_by_type[pharm_type] = []
            cluster_counts[pharm_type] = 0
            continue

        rc = radii.get(pharm_type, 2.0)

        lfds = calculate_lfds(molecule, atom_indices, rc)
        find, union = cluster_features(molecule, atom_indices, lfds, rc)
        centers = calculate_geometric_centers(molecule, atom_indices, union)
        pcs = compute_pca_nipals(molecule, atom_indices, union, centers)

        lfds_by_type[pharm_type] = lfds
        clusters_by_type[pharm_type] = (find, union)
        centers_by_type[pharm_type] = centers
        pcs_by_type[pharm_type] = pcs
        cluster_counts[pharm_type] = count_clusters(find)

    # Step 6: Gaussian descriptors (original 3 types only)
    descriptors = compute_pharmacophore_descriptors(centers_by_type, pcs_by_type)

    # Step 7: Summary
    summary = compute_feature_summary(features, cluster_counts, centers_by_type)

    return {
        "features": features,
        "lfds": lfds_by_type,
        "clusters": clusters_by_type,
        "centers": centers_by_type,
        "pcs": pcs_by_type,
        "cluster_counts": cluster_counts,
        "descriptors": descriptors,
        "summary": summary,
    }


# =============================================================================
# Component 2: Qt Graphics Visualization
# =============================================================================

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QGroupBox, QFileDialog, QMessageBox, QFormLayout,
    Qt, QColor, QCheckBox, QDialog, QApplication,
    QGraphicsView, QGraphicsScene, QGraphicsEllipseItem,
    QGraphicsLineItem, QGraphicsTextItem, QPen, QBrush, QFont,
    QPainter, QImage, QPointF, QRectF, QScrollArea,
)
from src.shared.ui.theme import COLORS
from src.plugins.base_plugin import BasePlugin
from src.plugins.plugin_types import PluginInfo, PluginType

try:
    from src.shared.qt_compat import QDoubleSpinBox
except ImportError:
    from PySide6.QtWidgets import QDoubleSpinBox

try:
    from src.shared.qt_compat import QSpinBox
except ImportError:
    from PySide6.QtWidgets import QSpinBox

try:
    from src.shared.qt_compat import QComboBox
except ImportError:
    from PySide6.QtWidgets import QComboBox

try:
    from src.shared.qt_compat import QTextEdit
except ImportError:
    from PySide6.QtWidgets import QTextEdit

try:
    from PySide6.QtGui import QRadialGradient, QLinearGradient, QCursor
except ImportError:
    QRadialGradient = None
    QLinearGradient = None
    QCursor = None

try:
    from PySide6.QtWidgets import QGraphicsRectItem, QToolButton, QSizePolicy
except ImportError:
    QGraphicsRectItem = None
    QToolButton = None
    QSizePolicy = None


# ── Element colors for ball-and-stick ────────────────────────────────────────

_ELEMENT_COLORS = {
    'C': '#55ff7f', 'N': '#3050F8', 'O': '#FF0D0D', 'S': '#FFFF30',
    'P': '#FF8000', 'H': '#d0d0d0', 'F': '#90E050', 'Cl': '#1FF01F',
    'Br': '#A62929', 'I': '#940094', 'Fe': '#E06633', 'Zn': '#7D80B0',
}

_ELEMENT_RADII = {
    'C': 5.0, 'N': 4.5, 'O': 4.5, 'S': 6.0,
    'P': 5.5, 'H': 2.5, 'F': 3.5, 'Cl': 5.5,
    'Br': 6.0, 'I': 6.5,
}


# ── 3D Rotation Matrix Helpers ───────────────────────────────────────────────

def _rotation_x(angle_deg: float) -> np.ndarray:
    """Rotation matrix around X axis."""
    rad = np.radians(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _rotation_y(angle_deg: float) -> np.ndarray:
    """Rotation matrix around Y axis."""
    rad = np.radians(angle_deg)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


class PharmacophoreGraphicsView(QGraphicsView):
    """
    Interactive QGraphicsView with 3D rotation, zooming, and distance
    measurement support.

    Controls:
        - Left-click drag: 3D rotation
        - Middle-click drag / scroll wheel: zoom
        - Right-click drag: pan
        - Ctrl+click: distance measurement between atoms
    """

    def __init__(self, parent_widget=None):
        super().__init__()
        self._parent_widget = parent_widget
        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.setRenderHint(QPainter.Antialiasing)
        self.setRenderHint(QPainter.TextAntialiasing)
        self.setDragMode(QGraphicsView.NoDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorUnderMouse)
        self.setBackgroundBrush(QBrush(Qt.white))
        self.setStyleSheet("border: none; border-radius: 8px;")
        self._zoom = 0
        self._last_mouse_pos = None
        self._is_rotating = False
        self._is_panning = False
        self._is_measuring = False

    def wheelEvent(self, event):
        """Scroll wheel zoom."""
        if event.angleDelta().y() > 0:
            self.scale(1.25, 1.25)
            self._zoom += 1
        else:
            self.scale(0.8, 0.8)
            self._zoom -= 1

    def mousePressEvent(self, event):
        """Handle mouse press for rotation, pan, and measurement."""
        if event.button() == Qt.LeftButton:
            # Check for measurement mode (either via button toggle or Ctrl+click)
            is_measure_active = False
            if self._parent_widget and getattr(self._parent_widget, "_measure_mode", False):
                is_measure_active = True
            
            if is_measure_active or (event.modifiers() & Qt.ControlModifier):
                self._handle_measurement_click(event)
                return

            # Left drag: 3D rotation
            self._is_rotating = True
            self._last_mouse_pos = event.pos()
            if QCursor:
                self.setCursor(QCursor(Qt.ClosedHandCursor))
        elif event.button() == Qt.MiddleButton:
            # Middle drag: zoom via drag
            self._last_mouse_pos = event.pos()
        elif event.button() == Qt.RightButton:
            # Right drag: pan
            self._is_panning = True
            self._last_mouse_pos = event.pos()
            self.setDragMode(QGraphicsView.NoDrag)
            if QCursor:
                self.setCursor(QCursor(Qt.OpenHandCursor))
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        """Handle mouse drag for rotation and zoom."""
        if self._last_mouse_pos is not None:
            delta = event.pos() - self._last_mouse_pos

            if self._is_rotating and self._parent_widget:
                # Rotate the 3D view
                dx = delta.x() * 0.5
                dy = delta.y() * 0.5
                self._parent_widget._rotate_view(dx, dy)
                self._last_mouse_pos = event.pos()
                return

            if event.buttons() & Qt.MiddleButton:
                # Middle mouse drag: zoom
                zoom_delta = delta.y() * -0.005
                factor = 1.0 + zoom_delta
                factor = max(0.5, min(2.0, factor))
                self.scale(factor, factor)
                self._last_mouse_pos = event.pos()
                return

            if self._is_panning:
                # Pan the view
                dx = delta.x()
                dy = delta.y()
                self.horizontalScrollBar().setValue(
                    self.horizontalScrollBar().value() - dx
                )
                self.verticalScrollBar().setValue(
                    self.verticalScrollBar().value() - dy
                )
                self._last_mouse_pos = event.pos()
                return

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        """Reset drag state."""
        self._is_rotating = False
        self._is_panning = False
        self._last_mouse_pos = None
        self.setCursor(Qt.ArrowCursor)
        super().mouseReleaseEvent(event)

    def _handle_measurement_click(self, event):
        """Handle Ctrl+click for distance measurement."""
        if self._parent_widget:
            scene_pos = self.mapToScene(event.pos())
            self._parent_widget._on_measure_click(scene_pos)

    def fit_content(self):
        rect = self.scene.itemsBoundingRect()
        if not rect.isNull():
            self.fitInView(rect.adjusted(-50, -50, 50, 50), Qt.KeepAspectRatio)
            self._zoom = 0


class PharmacophoreVisualizerWidget(QWidget):
    """
    Main widget for the Pharmacophore Features plugin.

    Left panel: controls (radii, filters, file loading, tools).
    Right panel: interactive 3D-projected visualization with ball-and-stick
    rendering, PCA ellipsoids, and distance measurement.
    Bottom: horizontal legend in 2-3 rows.
    """

    def __init__(self, plugin: 'PharmacophoreFeatures'):
        super().__init__()
        self.plugin = plugin
        self.widget = self
        self.molecule = None
        self._results = None

        # 3D rotation state
        self._rot_x = 0.0   # Rotation around X axis (degrees)
        self._rot_y = 0.0   # Rotation around Y axis (degrees)
        self._view_scale = 8.0

        # Distance measurement state
        self._measure_mode = False
        self._measure_atom1 = None  # (scene_x, scene_y, atom_idx)
        self._measure_items = []    # Track drawn measurement items

        # Ball-and-stick toggle
        self._show_ball_stick = True

        self._setup_ui()

    def _setup_ui(self):
        self.setStyleSheet(f"background-color: {COLORS['bg_secondary']};")
        outer_lay = QVBoxLayout(self)
        outer_lay.setContentsMargins(10, 10, 10, 10)
        outer_lay.setSpacing(8)

        # ── Top area: Left panel + viewer ────────────────────────
        top_lay = QHBoxLayout()
        top_lay.setSpacing(12)

        # ── Left Panel with Scrolling ────────────────────────────
        scroll = QScrollArea()
        scroll.setFixedWidth(310)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        scroll.setStyleSheet("background: transparent; border: none;")

        lp = QWidget()
        lp.setStyleSheet(
            f"background-color: {COLORS['bg_tertiary']}; border-radius: 10px;"
        )
        llay = QVBoxLayout(lp)
        llay.setContentsMargins(15, 18, 15, 18)
        llay.setSpacing(10)

        # ... (titles and controls) ...
        title = QLabel("PHARMACOPHORE")
        title.setStyleSheet(
            f"font-size: 20px; font-weight: 800; "
            f"color: {COLORS['accent']}; letter-spacing: 1px;"
        )
        llay.addWidget(title)

        subtitle = QLabel("Feature Detection & Clustering")
        subtitle.setStyleSheet(
            f"font-size: 11px; color: {COLORS['text_muted']};"
        )
        llay.addWidget(subtitle)

        # ── Radii Controls ───────────────────────────────────────
        radii_box = QGroupBox("CLUSTER RADII (\u00c5)")
        radii_box.setStyleSheet(
            f"QGroupBox {{ color: {COLORS['accent2']}; font-weight: bold; "
            f"border: 1px solid {COLORS['border']}; "
            f"margin-top: 12px; padding-top: 12px; }}"
        )
        rlay = QFormLayout()
        self._radii_spins = {}
        from src.shared.qt_compat import QDoubleSpinBox
        for pharm_type, default_rc in DEFAULT_RADII.items():
            spin = QDoubleSpinBox()
            spin.setRange(0.5, 8.0)
            spin.setSingleStep(0.1)
            spin.setValue(default_rc)
            spin.setDecimals(1)
            short_name = pharm_type.split()[0][:5]
            rlay.addRow(f"{short_name}:", spin)
            self._radii_spins[pharm_type] = spin
        radii_box.setLayout(rlay)
        llay.addWidget(radii_box)

        # ── Feature Type Filters ─────────────────────────────────
        filter_box = QGroupBox("FEATURE FILTERS")
        filter_box.setStyleSheet(
            f"QGroupBox {{ color: {COLORS['accent2']}; font-weight: bold; "
            f"border: 1px solid {COLORS['border']}; "
            f"margin-top: 12px; padding-top: 12px; }}"
        )
        flay = QVBoxLayout()
        self._filter_checks: Dict[str, QCheckBox] = {}
        for pharm_type, color in PHARMACOPHORE_COLORS.items():
            cb = QCheckBox(pharm_type)
            cb.setChecked(True)
            cb.setStyleSheet(f"color: {color}; font-weight: bold;")
            cb.stateChanged.connect(self._on_filter_changed)
            flay.addWidget(cb)
            self._filter_checks[pharm_type] = cb
        filter_box.setLayout(flay)
        llay.addWidget(filter_box)

        # ── Display options ──────────────────────────────────────
        disp_box = QGroupBox("DISPLAY")
        disp_box.setStyleSheet(
            f"QGroupBox {{ color: {COLORS['accent2']}; font-weight: bold; "
            f"border: 1px solid {COLORS['border']}; "
            f"margin-top: 12px; padding-top: 12px; }}"
        )
        dlay = QVBoxLayout()
        self._cb_ball_stick = QCheckBox("Ball && Stick")
        self._cb_ball_stick.setChecked(True)
        self._cb_ball_stick.stateChanged.connect(self._on_display_changed)
        dlay.addWidget(self._cb_ball_stick)

        self._cb_ellipses = QCheckBox("Show Ellipses")
        self._cb_ellipses.setChecked(True)
        self._cb_ellipses.stateChanged.connect(self._on_display_changed)
        dlay.addWidget(self._cb_ellipses)

        self._cb_features = QCheckBox("Show Feature Markers")
        self._cb_features.setChecked(True)
        self._cb_features.stateChanged.connect(self._on_display_changed)
        dlay.addWidget(self._cb_features)

        disp_box.setLayout(dlay)
        llay.addWidget(disp_box)

        # ── Action Buttons ───────────────────────────────────────
        btn_style = (
            f"QPushButton {{ background-color: {COLORS['accent']}; "
            f"color: white; font-weight: bold; padding: 10px; "
            f"border-radius: 7px; font-size: 12px; }}"
            f"QPushButton:hover {{ background-color: {COLORS['accent_hover']}; }}"
        )
        sec_style = (
            f"QPushButton {{ background-color: transparent; "
            f"border: 2px solid {COLORS['accent']}; "
            f"color: {COLORS['accent']}; font-weight: bold; "
            f"padding: 8px; border-radius: 7px; font-size: 11px; }}"
            f"QPushButton:hover {{ background-color: {COLORS['accent']}; "
            f"color: white; }}"
        )
        tool_style = (
            f"QPushButton {{ background-color: transparent; "
            f"border: 2px solid {COLORS['accent2']}; "
            f"color: {COLORS['accent2']}; font-weight: bold; "
            f"padding: 8px; border-radius: 7px; font-size: 11px; }}"
            f"QPushButton:hover {{ background-color: {COLORS['accent2']}; "
            f"color: white; }}"
        )

        btn_load_viewer = QPushButton("LOAD FROM VIEWER")
        btn_load_viewer.setStyleSheet(btn_style)
        btn_load_viewer.clicked.connect(self._load_current_molecule)
        llay.addWidget(btn_load_viewer)

        btn_load_file = QPushButton("OPEN FILE (MOL/SDF/MOL2)")
        btn_load_file.setStyleSheet(btn_style)
        btn_load_file.clicked.connect(self._load_from_file)
        llay.addWidget(btn_load_file)

        btn_run = QPushButton("ANALYZE FEATURES")
        btn_run.setStyleSheet(btn_style)
        btn_run.clicked.connect(self._run_analysis)
        llay.addWidget(btn_run)

        # ── Tool Buttons ─────────────────────────────────────────
        self._btn_measure = QPushButton("MEASURE DISTANCE")
        self._btn_measure.setStyleSheet(tool_style)
        self._btn_measure.setCheckable(True)
        self._btn_measure.clicked.connect(self._toggle_measure_mode)
        llay.addWidget(self._btn_measure)

        btn_reset_view = QPushButton("RESET VIEW")
        btn_reset_view.setStyleSheet(tool_style)
        btn_reset_view.clicked.connect(self._reset_view)
        llay.addWidget(btn_reset_view)

        # ── Export buttons row ───────────────────────────────────
        export_row_1 = QHBoxLayout()
        btn_export_json = QPushButton("JSON")
        btn_export_json.setStyleSheet(sec_style)
        btn_export_json.clicked.connect(self._export_json)
        export_row_1.addWidget(btn_export_json)

        btn_export_liquid = QPushButton(".LIQUID")
        btn_export_liquid.setStyleSheet(sec_style)
        btn_export_liquid.clicked.connect(self._export_liquid)
        export_row_1.addWidget(btn_export_liquid)
        llay.addLayout(export_row_1)

        export_row_2 = QHBoxLayout()
        btn_export_csv = QPushButton("CSV")
        btn_export_csv.setStyleSheet(sec_style)
        btn_export_csv.clicked.connect(self._export_csv)
        export_row_2.addWidget(btn_export_csv)

        btn_export_img = QPushButton("IMAGE")
        btn_export_img.setStyleSheet(sec_style)
        btn_export_img.clicked.connect(self._export_image)
        export_row_2.addWidget(btn_export_img)
        llay.addLayout(export_row_2)

        llay.addStretch()

        # ── Stats Label ──────────────────────────────────────────
        self.lbl_stats = QLabel("Ready")
        self.lbl_stats.setStyleSheet(
            f"color: {COLORS['text_muted']}; font-size: 10px;"
        )
        self.lbl_stats.setWordWrap(True)
        llay.addWidget(self.lbl_stats)

        scroll.setWidget(lp)
        top_lay.addWidget(scroll)

        # ── Right Panel: Graphics View ───────────────────────────
        self.viewer = PharmacophoreGraphicsView(parent_widget=self)
        top_lay.addWidget(self.viewer)

        outer_lay.addLayout(top_lay)

        # ── Bottom Legend (horizontal, 2-3 rows) ─────────────────
        self._legend_widget = QWidget()
        self._legend_widget.setStyleSheet(
            f"background-color: {COLORS['bg_tertiary']}; "
            f"border-radius: 8px; padding: 6px;"
        )
        self._legend_widget.setFixedHeight(60)
        self._build_legend()
        outer_lay.addWidget(self._legend_widget)

    def _build_legend(self):
        """Build the bottom horizontal legend with interactive toggles."""
        if hasattr(self, "_legend_layout"):
            # Clear existing legend layout
            while self._legend_layout.count():
                item = self._legend_layout.takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
        else:
            self._legend_layout = QHBoxLayout(self._legend_widget)
            self._legend_layout.setContentsMargins(12, 4, 12, 4)
            self._legend_layout.setSpacing(6)

        items = list(PHARMACOPHORE_COLORS.items())
        # Row-based splitting for readability
        row1_items = items[:3]
        row2_items = items[3:6]
        row3_items = items[6:]

        col_lay = QVBoxLayout()
        col_lay.setSpacing(3)

        for row_items in [row1_items, row2_items, row3_items]:
            if not row_items: continue
            row_lay = QHBoxLayout()
            row_lay.setSpacing(14)
            for pharm_type, color_hex in row_items:
                cb = QCheckBox(pharm_type)
                cb.setChecked(self._filter_checks[pharm_type].isChecked())
                # Sync legend checkbox with sidebar filter
                cb.stateChanged.connect(
                    lambda state, pt=pharm_type: self._sync_filters(pt, state)
                )
                cb.setStyleSheet(f"""
                    QCheckBox {{ color: {COLORS['text_primary']}; font-size: 11px; font-weight: 500; }}
                    QCheckBox::indicator {{ width: 12px; height: 12px; border-radius: 6px; }}
                    QCheckBox::indicator:unchecked {{ background-color: transparent; border: 2px solid {color_hex}; }}
                    QCheckBox::indicator:checked {{ background-color: {color_hex}; border: 2px solid {color_hex}; }}
                """)
                row_lay.addWidget(cb)
            row_lay.addStretch()
            col_lay.addLayout(row_lay)

        self._legend_layout.addLayout(col_lay)

    def _sync_filters(self, pharm_type, state):
        """Synchronize legend and sidebar filter checkboxes."""
        self._filter_checks[pharm_type].blockSignals(True)
        self._filter_checks[pharm_type].setChecked(state == Qt.Checked)
        self._filter_checks[pharm_type].blockSignals(False)
        self._render_visualization()

    # ── Actions ──────────────────────────────────────────────────

    def _load_current_molecule(self):
        """Load molecule from the main PyChem viewer."""
        if not self.plugin or not self.plugin.api:
            QMessageBox.warning(self, "No API", "Plugin API not connected.")
            return
        mol = self.plugin.get_current_molecule()
        if mol:
            self.molecule = mol
            self._rot_x = 0.0
            self._rot_y = 0.0
            self.lbl_stats.setText(f"Loaded: {len(mol.atoms)} atoms")
            self._run_analysis()
        else:
            QMessageBox.information(
                self, "No Molecule",
                "No molecule loaded in the main viewer.\n"
                "Use 'Open File' to load a MOL/SDF/MOL2 file directly."
            )

    def _load_from_file(self):
        """Load molecule from a MOL, SDF, or MOL2 file."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Molecule File", "",
            "Molecule Files (*.mol *.sdf *.mol2);;MOL Files (*.mol);;"
            "SDF Files (*.sdf);;MOL2 Files (*.mol2);;All Files (*)"
        )
        if not path:
            return

        try:
            from src.features.io.loaders.file_reader import (
                read_mol, read_sdf, read_mol2
            )
            ext = os.path.splitext(path)[1].lower()
            if ext == '.mol2':
                mol = read_mol2(path)
            elif ext == '.sdf':
                mol = read_sdf(path)
            elif ext == '.mol':
                mol = read_mol(path)
            else:
                QMessageBox.warning(
                    self, "Unknown Format",
                    f"Unsupported file format: {ext}\n"
                    "Supported: .mol, .sdf, .mol2"
                )
                return

            self.molecule = mol
            self._rot_x = 0.0
            self._rot_y = 0.0
            self.lbl_stats.setText(
                f"Loaded: {mol.name or os.path.basename(path)} "
                f"({len(mol.atoms)} atoms)"
            )
            self._run_analysis()

        except Exception as e:
            logging.error(f"File load error: {e}", exc_info=True)
            QMessageBox.critical(self, "Load Error", str(e))

    def set_molecule(self, molecule):
        self.molecule = molecule
        self._rot_x = 0.0
        self._rot_y = 0.0
        # Disabled automatic analysis to prevent UI hang on large proteins.
        # User must click "Run Analysis" or "Load from Viewer" manually.
        if molecule:
            self.lbl_stats.setText(f"Ready to analyze: {len(molecule.atoms)} atoms")
        else:
            self.lbl_stats.setText("No molecule loaded")

    def _get_radii(self) -> Dict[str, float]:
        return {t: s.value() for t, s in self._radii_spins.items()}

    def _run_analysis(self):
        if not self.molecule:
            QMessageBox.warning(self, "No Molecule", "Load a molecule first.")
            return

        has_3d = any(a.has_coords for a in self.molecule.atoms)
        if not has_3d:
            QMessageBox.warning(
                self, "No 3D Coordinates",
                "The molecule has no 3D coordinates.\n"
                "Load a PDB/MOL2/SDF file with 3D coords."
            )
            return

        try:
            radii = self._get_radii()
            self._results = run_pharmacophore_pipeline(self.molecule, radii)
            self._render_visualization()
            self._update_stats()
        except Exception as e:
            logging.error(f"Pharmacophore analysis error: {e}", exc_info=True)
            QMessageBox.critical(self, "Analysis Error", str(e))

    def _on_filter_changed(self, _state):
        """Re-render when filter checkboxes change."""
        if self._results:
            self._render_visualization()

    def _on_display_changed(self, _state):
        """Re-render when display options change."""
        if self._results:
            self._render_visualization()

    def _toggle_measure_mode(self, checked):
        """Toggle distance measurement mode."""
        self._measure_mode = checked
        self._measure_atom1 = None
        if checked:
            self._btn_measure.setStyleSheet(
                f"QPushButton {{ background-color: {COLORS['accent2']}; "
                f"color: white; font-weight: bold; "
                f"padding: 8px; border-radius: 7px; font-size: 11px; }}"
            )
            self.lbl_stats.setText("Ctrl+click on two atoms to measure distance")
        else:
            self._btn_measure.setStyleSheet(
                f"QPushButton {{ background-color: transparent; "
                f"border: 2px solid {COLORS['accent2']}; "
                f"color: {COLORS['accent2']}; font-weight: bold; "
                f"padding: 8px; border-radius: 7px; font-size: 11px; }}"
                f"QPushButton:hover {{ background-color: {COLORS['accent2']}; "
                f"color: white; }}"
            )

    def _reset_view(self):
        """Reset 3D rotation to default."""
        self._rot_x = 0.0
        self._rot_y = 0.0
        self._clear_measurements()
        if self._results:
            self._render_visualization()
            self.viewer.fit_content()

    def _clear_measurements(self):
        """Remove all distance measurement items."""
        for item in self._measure_items:
            if item.scene():
                self.viewer.scene.removeItem(item)
        self._measure_items.clear()
        self._measure_atom1 = None

    def _rotate_view(self, dx: float, dy: float):
        """Rotate the 3D view by dx/dy degrees and re-render."""
        self._rot_y += dx
        self._rot_x += dy
        if self._results:
            self._render_visualization()

    def _update_stats(self):
        if not self._results:
            return
        s = self._results["summary"]
        fc = s["feature_counts"]
        cc = s["cluster_counts"]
        lines = [
            f"<b>Features:</b> {s['total_features']}  "
            f"<b>Clusters:</b> {s['total_clusters']}",
        ]
        for t in PHARMACOPHORE_COLORS:
            c = PHARMACOPHORE_COLORS[t]
            lines.append(
                f"<span style='color:{c};'>\u25cf</span> "
                f"<b>{t.split()[0][:6]}</b>: "
                f"{fc.get(t, 0)} \u2192 {cc.get(t, 0)} cl"
            )
        self.lbl_stats.setText(
            f"<div style='font-family: monospace; font-size: 10px;'>"
            + "<br/>".join(lines)
            + "</div>"
        )

    # ── 3D Projection ────────────────────────────────────────────

    def _project_3d(self, coords_3d: np.ndarray, reference_centroid: np.ndarray) -> np.ndarray:
        """
        Project 3D coordinates to 2D using a fixed reference centroid.
        """
        centered = coords_3d - reference_centroid

        # Use molecule atoms for PCA alignment if possible
        if not hasattr(self, "_view_matrix"):
            all_coords = np.array([[a.x, a.y, a.z] for a in self.molecule.atoms if a.has_coords])
            mol_centered = all_coords - np.mean(all_coords, axis=0)
            cov = mol_centered.T @ mol_centered
            _, eigvecs = np.linalg.eigh(cov)
            self._view_matrix = eigvecs[:, ::-1] # Save PCs in descending order

        aligned = centered @ self._view_matrix

        # Apply user rotation
        rot = _rotation_x(self._rot_x) @ _rotation_y(self._rot_y)
        rotated = aligned @ rot.T

        # Orthographic projection: take XY
        return rotated[:, :2] * self._view_scale, rotated[:, 2]  # Return project coords and Z-depth

    # ── Measurement ──────────────────────────────────────────────

    def _on_measure_click(self, scene_pos):
        """Handle a measurement click at the given scene position."""
        if not self._measure_mode or not self.molecule:
            return

        # Find nearest atom to click position
        nearest_idx = None
        nearest_dist = float('inf')

        if not hasattr(self, '_proj_map') or not self._proj_map:
            return

        for atom_idx, (px, py) in self._proj_map.items():
            dx = scene_pos.x() - px
            dy = scene_pos.y() - py
            d = math.sqrt(dx * dx + dy * dy)
            if d < nearest_dist and d < 15.0:  # 15px click radius
                nearest_dist = d
                nearest_idx = atom_idx

        if nearest_idx is None:
            return

        if self._measure_atom1 is None:
            # First atom selected
            self._measure_atom1 = nearest_idx
            px, py = self._proj_map[nearest_idx]
            # Draw selection ring
            ring = QGraphicsEllipseItem(-8, -8, 16, 16)
            ring.setPos(px, py)
            ring.setPen(QPen(QColor("#FF6600"), 2.5, Qt.SolidLine))
            ring.setBrush(QBrush(Qt.NoBrush))
            ring.setZValue(50)
            self.viewer.scene.addItem(ring)
            self._measure_items.append(ring)
            self.lbl_stats.setText(
                f"Selected atom {nearest_idx} "
                f"({self.molecule.atoms[nearest_idx].symbol}). "
                f"Ctrl+click second atom."
            )
        else:
            # Second atom — compute and draw distance
            a1_idx = self._measure_atom1
            a2_idx = nearest_idx
            a1 = self.molecule.atoms[a1_idx]
            a2 = self.molecule.atoms[a2_idx]

            # 3D distance
            dist_3d = math.sqrt(
                (a1.x - a2.x) ** 2 +
                (a1.y - a2.y) ** 2 +
                (a1.z - a2.z) ** 2
            )

            p1 = self._proj_map[a1_idx]
            p2 = self._proj_map[a2_idx]

            # Draw selection ring on second atom
            ring2 = QGraphicsEllipseItem(-8, -8, 16, 16)
            ring2.setPos(p2[0], p2[1])
            ring2.setPen(QPen(QColor("#FF6600"), 2.5, Qt.SolidLine))
            ring2.setBrush(QBrush(Qt.NoBrush))
            ring2.setZValue(50)
            self.viewer.scene.addItem(ring2)
            self._measure_items.append(ring2)

            # Draw dashed line
            line = QGraphicsLineItem(p1[0], p1[1], p2[0], p2[1])
            line.setPen(QPen(QColor("#FF6600"), 1.5, Qt.DashDotLine))
            line.setZValue(49)
            self.viewer.scene.addItem(line)
            self._measure_items.append(line)

            # Draw distance label at midpoint
            mx = (p1[0] + p2[0]) / 2
            my = (p1[1] + p2[1]) / 2
            label = QGraphicsTextItem(f"{dist_3d:.2f} \u00c5")
            label.setFont(QFont("Segoe UI", 8, QFont.Bold))
            label.setDefaultTextColor(QColor("#FF6600"))
            label.setPos(mx + 5, my - 12)
            label.setZValue(51)
            self.viewer.scene.addItem(label)
            self._measure_items.append(label)

            self.lbl_stats.setText(
                f"Distance: {a1.symbol}{a1_idx} \u2194 "
                f"{a2.symbol}{a2_idx} = <b>{dist_3d:.3f} \u00c5</b>"
            )
            self._measure_atom1 = None

    # ── Rendering ────────────────────────────────────────────────

    def _render_visualization(self):
        """Render pharmacophore features with 3D rotation and ball-and-stick."""
        if not self._results or not self.molecule:
            return

        # Ensure reference view matrix is reset if molecule changed
        if not hasattr(self, "_molecule_id") or self._molecule_id != id(self.molecule):
            if hasattr(self, "_view_matrix"): delattr(self, "_view_matrix")
            self._molecule_id = id(self.molecule)

        # Preserve measurement items
        preserved_items = list(self._measure_items)
        self.viewer.scene.clear()
        self._measure_items.clear()

        features = self._results["features"]
        centers = self._results["centers"]
        pcs = self._results["pcs"]

        active_types = {
            t for t, cb in self._filter_checks.items() if cb.isChecked()
        }

        show_ball_stick = self._cb_ball_stick.isChecked()
        show_ellipses = self._cb_ellipses.isChecked()
        show_features = self._cb_features.isChecked()

        # Collect all 3D coordinates for reference centroid
        all_coords = []
        atom_list = []
        for a in self.molecule.atoms:
            if a.has_coords:
                # Suppress H atoms in ball-and-stick view
                if show_ball_stick and a.symbol == 'H':
                    continue
                all_coords.append([a.x, a.y, a.z])
                atom_list.append(a)
        
        if not all_coords:
            return
        all_coords = np.array(all_coords)
        ref_centroid = np.mean(all_coords, axis=0)

        # Project all atoms
        projected, depths = self._project_3d(all_coords, ref_centroid)

        self._proj_map = {}
        depth_map = {}
        for i, a in enumerate(atom_list):
            self._proj_map[a.index] = (projected[i, 0], projected[i, 1])
            depth_map[a.index] = depths[i]

        # ── Draw bonds ───────────────────────────────────────────
        if show_ball_stick:
            bond_pen = QPen(QColor("#9E9E9E"), 2.5, Qt.SolidLine, Qt.RoundCap)
            for bond in self.molecule.bonds:
                a1, a2 = bond.begin_atom_idx, bond.end_atom_idx
                if a1 in self._proj_map and a2 in self._proj_map:
                    p1, p2 = self._proj_map[a1], self._proj_map[a2]
                    avg_depth = (depth_map[a1] + depth_map[a2]) / 2
                    line = self.viewer.scene.addLine(p1[0], p1[1], p2[0], p2[1], bond_pen)
                    line.setZValue(10 + avg_depth * 0.1)

            # Draw atoms as balls (sorted back-to-front)
            sorted_atom_indices = sorted(self._proj_map.keys(), key=lambda idx: depth_map[idx])
            for idx in sorted_atom_indices:
                atom = self.molecule.atoms[idx]
                px, py = self._proj_map[idx]
                radius = _ELEMENT_RADII.get(atom.symbol, 4.0)
                color = QColor(_ELEMENT_COLORS.get(atom.symbol, '#808080'))
                
                ball = QGraphicsEllipseItem(-radius, -radius, 2*radius, 2*radius)
                ball.setPos(px, py)
                if QRadialGradient:
                    grad = QRadialGradient(QPointF(-radius*0.3, -radius*0.3), radius*1.5)
                    grad.setColorAt(0.0, color.lighter(170))
                    grad.setColorAt(0.6, color)
                    grad.setColorAt(1.0, color.darker(150))
                    ball.setBrush(QBrush(grad))
                else:
                    ball.setBrush(QBrush(color))
                ball.setPen(QPen(color.darker(130), 0.5))
                ball.setZValue(20 + depth_map[idx]*0.1)
                self.viewer.scene.addItem(ball)
        else:
            # Wireframe mode: thin lines + small dots
            bond_pen = QPen(QColor("#BDBDBD"), 1.0, Qt.SolidLine, Qt.RoundCap)
            for bond in self.molecule.bonds:
                a1, a2 = bond.begin_atom_idx, bond.end_atom_idx
                if a1 in self._proj_map and a2 in self._proj_map:
                    p1, p2 = self._proj_map[a1], self._proj_map[a2]
                    self.viewer.scene.addLine(p1[0], p1[1], p2[0], p2[1], bond_pen)
            
            for atom_idx in self._proj_map:
                px, py = self._proj_map[atom_idx]
                symbol = self.molecule.atoms[atom_idx].symbol
                color = QColor(_ELEMENT_COLORS.get(symbol, '#808080'))
                dot = QGraphicsEllipseItem(-2, -2, 4, 4)
                dot.setPos(px, py)
                dot.setBrush(QBrush(color))
                dot.setPen(QPen(Qt.NoPen))
                dot.setZValue(15)
                self.viewer.scene.addItem(dot)

        # ── Draw cluster ellipses ────────────────────────────────
        if show_ellipses:
            for pharm_type in PHARMACOPHORE_COLORS:
                if pharm_type not in active_types:
                    continue
                type_centers = centers.get(pharm_type, [])
                type_pcs = pcs.get(pharm_type, [])
                type_color = QColor(PHARMACOPHORE_COLORS[pharm_type])

                for ci, ctr_3d in enumerate(type_centers):
                    if ctr_3d is None: continue
                    
                    # Project cluster center using fixed ref_centroid
                    c_2d, c_depth = self._project_3d(np.array([ctr_3d]), ref_centroid)
                    cx, cy, cz = c_2d[0, 0], c_2d[0, 1], c_depth[0]

                    # Project PCA axes
                    pc_axes = type_pcs[ci]
                    axis_endpoints = np.array([ctr_3d + pc for pc in pc_axes])
                    ends_2d, _ = self._project_3d(axis_endpoints, ref_centroid)
                    proj_axes = ends_2d - c_2d

                    a_len = max(np.linalg.norm(proj_axes[0]), 4.0)
                    b_len = max(np.linalg.norm(proj_axes[1]), 4.0) if len(proj_axes) > 1 else a_len
                    angle = math.degrees(math.atan2(proj_axes[0,1], proj_axes[0,0])) if a_len > 0.1 else 0

                    ellipse = QGraphicsEllipseItem(-a_len, -b_len, 2*a_len, 2*b_len)
                    ellipse.setPos(cx, cy)
                    ellipse.setRotation(angle)
                    
                    fill = QColor(type_color)
                    fill.setAlpha(pharm_type == PHARM_LIPOPHILIC and 60 or 35)
                    ellipse.setBrush(QBrush(fill))
                    ellipse.setPen(QPen(type_color, pharm_type == PHARM_LIPOPHILIC and 2.5 or 1.5, 
                                       pharm_type == PHARM_LIPOPHILIC and Qt.SolidLine or Qt.DashLine))
                    ellipse.setZValue(120 + cz*0.1) # Place ellipses in front
                    self.viewer.scene.addItem(ellipse)

        # ── Draw feature atom markers ────────────────────────────
        if show_features:
            for pharm_type, atom_indices in features.items():
                if pharm_type not in active_types: continue
                color = QColor(PHARMACOPHORE_COLORS[pharm_type])
                for idx in atom_indices:
                    if idx not in self._proj_map: continue
                    px, py = self._proj_map[idx]
                    cz = depth_map[idx]
                    
                    marker = QGraphicsEllipseItem(-4, -4, 8, 8)
                    marker.setPos(px, py)
                    marker.setBrush(QBrush(color))
                    marker.setPen(QPen(Qt.white, 0.8))
                    marker.setZValue(150 + cz*0.1)
                    self.viewer.scene.addItem(marker)

        # Re-add measurements
        for item in preserved_items:
            # Re-verify and re-add in case some measure items were lost
            if item.scene() is None:
                self.viewer.scene.addItem(item)
            self._measure_items.append(item)

        if not self.viewer._is_rotating:
            self.viewer.fit_content()

    # ── Export Methods ────────────────────────────────────────────

    def _export_json(self):
        if not self._results:
            QMessageBox.warning(self, "No Data", "Run analysis first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export JSON", "pharmacophore.json",
            "JSON Files (*.json);;All Files (*)"
        )
        if not path:
            return
        try:
            export_data = {
                "summary": self._results["summary"],
                "features": {
                    t: indices for t, indices in self._results["features"].items()
                },
                "clusters": {},
            }
            for t in PHARMACOPHORE_COLORS:
                ctrs = self._results["centers"].get(t, [])
                pc_list = self._results["pcs"].get(t, [])
                type_clusters = []
                for ci, c in enumerate(ctrs):
                    if c is None:
                        continue
                    cluster_info = {
                        "center": c.tolist(),
                        "pca_axes": (
                            [pc.tolist() for pc in pc_list[ci]]
                            if ci < len(pc_list) and pc_list[ci]
                            else []
                        ),
                    }
                    type_clusters.append(cluster_info)
                if type_clusters:
                    export_data["clusters"][t] = type_clusters

            with open(path, "w") as f:
                json.dump(export_data, f, indent=2)
            QMessageBox.information(
                self, "Exported", f"Saved to {os.path.basename(path)}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def _export_liquid(self):
        """Export 120-D descriptor in original .liquid format."""
        if not self._results:
            QMessageBox.warning(self, "No Data", "Run analysis first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Descriptor", "molecule.liquid",
            "LIQUID Files (*.liquid);;All Files (*)"
        )
        if not path:
            return
        try:
            descs = self._results["descriptors"]
            with open(path, "w") as f:
                f.write(os.path.basename(path))
                for pair_name in ["L-L", "L-D", "L-A", "D-D", "D-A", "A-A"]:
                    bins = descs.get(pair_name, np.zeros(20))
                    for val in bins:
                        f.write(f"\t{val}")
                f.write("\n")
            QMessageBox.information(
                self, "Exported", f"Saved to {os.path.basename(path)}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def _export_csv(self):
        """Export feature/cluster summary as CSV."""
        if not self._results:
            QMessageBox.warning(self, "No Data", "Run analysis first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export CSV", "pharmacophore_summary.csv",
            "CSV Files (*.csv);;All Files (*)"
        )
        if not path:
            return
        try:
            s = self._results["summary"]
            fc = s["feature_counts"]
            cc = s["cluster_counts"]
            with open(path, "w") as f:
                f.write("Feature Type,Atom Count,Cluster Count\n")
                for t in PHARMACOPHORE_COLORS:
                    f.write(f"{t},{fc.get(t, 0)},{cc.get(t, 0)}\n")
                f.write(f"\nTotal,{s['total_features']},{s['total_clusters']}\n")
            QMessageBox.information(
                self, "Exported", f"Saved to {os.path.basename(path)}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def _export_image(self):
        """Export the visualization as a high-res PNG."""
        if not self.viewer.scene.items():
            QMessageBox.warning(self, "Empty", "Nothing to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Image", "pharmacophore.png",
            "PNG Images (*.png);;All Files (*)"
        )
        if not path:
            return
        try:
            source_rect = self.viewer.scene.itemsBoundingRect().adjusted(
                -20, -20, 20, 20
            )
            scale = 3.0
            w = int(source_rect.width() * scale)
            h = int(source_rect.height() * scale)
            image = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
            image.fill(Qt.white)
            painter = QPainter(image)
            painter.setRenderHint(QPainter.Antialiasing)
            painter.setRenderHint(QPainter.TextAntialiasing)
            target_rect = QRectF(0, 0, w, h)
            self.viewer.scene.render(painter, target_rect, source_rect)
            painter.end()
            dpm = int(300 / 0.0254)
            image.setDotsPerMeterX(dpm)
            image.setDotsPerMeterY(dpm)
            if image.save(path):
                QMessageBox.information(
                    self, "Exported",
                    f"Saved: {os.path.basename(path)} ({w}\u00d7{h} @ 300 DPI)"
                )
            else:
                raise Exception("QImage.save failed.")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))


# =============================================================================
# Component 3: Plugin Class
# =============================================================================

class PharmacophoreFeatures(BasePlugin):
    """
    PyChem plugin for pharmacophore feature detection, clustering,
    PCA analysis, and descriptor generation.

    Ported from liquid_1.0 (PyMOL) with enhancements for aromatic
    and ionizable feature types, interactive 3D-rotatable visualization,
    ball-and-stick rendering, distance measurement, and multiple export
    formats (JSON, .liquid, CSV, PNG).

    Supported molecule formats: MOL, SDF, MOL2 (via built-in file readers).
    """

    def __init__(self):
        info = PluginInfo(
            name="Pharmacophore-Features",
            version="1.0.0",
            description=(
                "Pharmacophore feature detection with LFD clustering, "
                "PCA ellipsoids, Gaussian descriptors, 3D rotation, "
                "ball-and-stick rendering, and distance measurement. "
                "Supports MOL/SDF/MOL2. Ported from liquid_1.0."
            ),
            author="Dr. Vijay Masand",
            plugin_type=PluginType.ANALYSIS,
            keywords=["pharmacophore", "clustering", "PCA", "descriptor",
                       "ball-and-stick", "3D"],
        )
        super().__init__(info)
        self.main_widget: Optional[PharmacophoreVisualizerWidget] = None

    def create_widget(self) -> PharmacophoreVisualizerWidget:
        self.main_widget = PharmacophoreVisualizerWidget(self)

        # Auto-load molecule if available
        if self.api:
            mol = self.get_current_molecule()
            if mol:
                logging.info(
                    f"PharmacophoreFeatures: Auto-loading "
                    f"{len(mol.atoms)} atoms"
                )
                self.main_widget.set_molecule(mol)
            else:
                self.main_widget.lbl_stats.setText(
                    "Click LOAD FROM VIEWER or OPEN FILE to begin"
                )
        return self.main_widget

    def initialize(self, main_window=None, api=None) -> bool:
        try:
            if main_window is not None and api is not None:
                super().initialize(main_window, api)
            else:
                self._main_window = main_window
                self._api = api
                self._is_initialized = True
            self.logger.info("Pharmacophore-Features v1.1.0 initialized")
            return True
        except Exception as e:
            self.logger.error(f"Initialization failed: {e}")
            return False

    def cleanup(self):
        if self.main_widget:
            self.main_widget.deleteLater()
            self.main_widget = None
        self.logger.info("Pharmacophore-Features plugin cleaned up")

    def on_molecule_changed(self, molecule):
        if self.main_widget:
            self.main_widget.set_molecule(molecule)

