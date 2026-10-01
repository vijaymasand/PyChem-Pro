"""
Mesh Builder — Faithful port of CPPCartoon's createSegmentMesh.

The key insight from CPPCartoon is that mesh generation works on a
4-PeptidePlane sliding window. For each profile vertex, TWO separate
B-splines are computed (one for start profile, one for end profile),
then blended using easing functions.

Includes vertex deduplication via spatial hashing for watertight meshes.
"""

import numpy as np
from typing import List, Tuple, Optional

from src.features.visualization_3d.services.peptide_plane import (
    PeptidePlane, HELIX, STRAND, COIL, transition
)
from src.features.visualization_3d.services.spline_math import (
    spline_for_planes, linear, in_out_quad, out_circ
)
from src.features.visualization_3d.services.profiles import (
    segment_profiles, segment_colors
)


class MeshAccumulator:
    """
    Accumulates vertices, colors, and triangles with vertex deduplication.
    Matches CPPCartoon's Mesh struct with verticesDict.
    """
    
    def __init__(self, tolerance=1e-4):
        self.vertices = []
        self.colors = []
        self.triangles = []
        self._vertex_dict = {}
        self._tolerance = tolerance
    
    def _quantize(self, v):
        """Quantize vertex position for spatial hashing."""
        t = self._tolerance
        return (
            round(v[0] / t) * t,
            round(v[1] / t) * t,
            round(v[2] / t) * t,
        )
    
    def add_vertex(self, pos, color):
        """Add a vertex, returning its index (deduplicated)."""
        key = self._quantize(pos)
        if key in self._vertex_dict:
            return self._vertex_dict[key]
        
        idx = len(self.vertices)
        self.vertices.append(pos.copy())
        self.colors.append(color.copy())
        self._vertex_dict[key] = idx
        return idx
    
    def add_quad(self, p1, p2, p3, p4, c1, c2, c3, c4):
        """
        Add a quad as 2 triangles (matching CPPCartoon triangulateQuad).
        Winding order: p10, p11, p01, p00 from CPPCartoon line 416.
        """
        id1 = self.add_vertex(p1, c1)
        id2 = self.add_vertex(p2, c2)
        id3 = self.add_vertex(p3, c3)
        id4 = self.add_vertex(p4, c4)
        
        self.triangles.append([id1, id2, id3])
        self.triangles.append([id1, id3, id4])
    
    def to_arrays(self):
        """Convert to numpy arrays."""
        if not self.vertices:
            return None, None, None
        return (
            np.array(self.vertices, dtype=np.float64),
            np.array(self.triangles, dtype=np.int32),
            np.array(self.colors, dtype=np.float64),
        )


def create_segment_mesh(mesh, pp1, pp2, pp3, pp4, spline_steps, profile_detail):
    """
    Generate mesh for one backbone segment using 4 peptide planes.
    Matches CPPCartoon cartoon.cpp createSegmentMesh() lines 335-428 exactly.
    
    Args:
        mesh: MeshAccumulator to add geometry to
        pp1..pp4: the 4-plane sliding window
        spline_steps: number of interpolation steps per segment
        profile_detail: number of points in cross-section profile
    """
    type0 = pp2.ss1  # SS of residue before this segment
    type1, type2 = transition(pp2)

    c1, c2 = segment_colors(pp2)

    # Get start and end profiles for this segment
    profile1, profile2 = segment_profiles(pp2, pp3, profile_detail)

    # Select easing function for profile blending
    if type1 == STRAND and type2 != STRAND:
        ease_func = linear  # Arrow tip: linear collapse
    else:
        ease_func = in_out_quad  # Default: smooth ease

    if type0 == STRAND and type1 != STRAND:
        ease_func = out_circ  # Coming out of strand: circular ease

    n = profile_detail

    # Build splines for each profile vertex through all 4 planes
    # splines1[j] = B-spline for profile1 vertex j → shape (spline_steps+1, 3)
    # splines2[j] = B-spline for profile2 vertex j → shape (spline_steps+1, 3)
    splines1 = []
    splines2 = []
    for j in range(n):
        p1_pt = profile1[j]
        p2_pt = profile2[j]
        s1 = spline_for_planes(pp1, pp2, pp3, pp4, spline_steps, p1_pt[0], p1_pt[1])
        s2 = spline_for_planes(pp1, pp2, pp3, pp4, spline_steps, p2_pt[0], p2_pt[1])
        splines1.append(s1)
        splines2.append(s2)

    # Generate quads by stepping along the spline and around the profile
    for i in range(spline_steps):
        t0 = ease_func(float(i) / spline_steps)
        t1 = ease_func(float(i + 1) / spline_steps)

        # Cap the arrow tip at the start of strand-to-coil transition
        if i == 0 and type1 == STRAND and type2 != STRAND:
            p00 = splines1[0][i]
            p10 = splines1[n // 4][i]
            p11 = splines1[2 * n // 4][i]
            p01 = splines1[3 * n // 4][i]
            mesh.add_quad(p00, p01, p11, p10, c1, c1, c1, c1)

        for j in range(n):
            j_next = (j + 1) % n

            # 4 corners of the quad from both profile splines
            p100 = splines1[j][i]
            p101 = splines1[j][i + 1]
            p110 = splines1[j_next][i]
            p111 = splines1[j_next][i + 1]

            p200 = splines2[j][i]
            p201 = splines2[j][i + 1]
            p210 = splines2[j_next][i]
            p211 = splines2[j_next][i + 1]

            # Blend between start and end profiles using easing
            p00 = p100 * (1 - t0) + p200 * t0
            p01 = p101 * (1 - t1) + p201 * t1
            p10 = p110 * (1 - t0) + p210 * t0
            p11 = p111 * (1 - t1) + p211 * t1

            # Interpolate colors
            c00 = c1 * (1 - t0) + c2 * t0
            c01 = c1 * (1 - t1) + c2 * t1
            c10 = c1 * (1 - t0) + c2 * t0
            c11 = c1 * (1 - t1) + c2 * t1

            # Add quad (winding order matches CPPCartoon line 416)
            mesh.add_quad(p10, p11, p01, p00, c10, c11, c01, c00)


def diff_pp(id1, id2):
    """Check if two residue IDs indicate a discontinuity. Matches CPPCartoon diffPP."""
    diff = 1
    if id1 < 0 and id2 > 0:
        diff = 2
    return (id2 - id1) > diff


def discontinuity(pp1, pp2, pp3, pp4):
    """Check if the 4-plane window spans a chain break. Matches CPPCartoon."""
    for pp in [pp1, pp2, pp3, pp4]:
        if diff_pp(pp.res_id1, pp.res_id2) or diff_pp(pp.res_id2, pp.res_id3):
            return True
    return False
