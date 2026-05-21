"""
Cross-Section Profiles — Faithful port of CPPCartoon's profile functions.

Generates the 2D cross-section shapes that are swept along the backbone:
  - ellipseProfile: round tubes for coils
  - rectangleProfile: flat sheets with arrow heads
  - roundedRectangleProfile: flattened ribbons for helices

Also includes segmentProfiles() which selects and configures the correct
profiles for a segment based on SS transitions.

Constants match CPPCartoon cartoon.h exactly.
"""

import numpy as np
from src.features.visualization_3d.services.peptide_plane import (
    HELIX, STRAND, COIL, transition
)

# ─── Constants from CPPCartoon cartoon.h ─────────────────────────────────────

RIBBON_WIDTH      = 2.0
RIBBON_HEIGHT     = 0.125
RIBBON_OFFSET     = 1.5
ARROW_HEAD_WIDTH  = 3.0
ARROW_WIDTH       = 2.0
ARROW_HEIGHT      = 0.5
TUBE_SIZE         = 0.25


# ─── Profile Generators ─────────────────────────────────────────────────────

def ellipse_profile(n, w, h):
    """
    Elliptical cross-section for coils/tubes.
    Matches CPPCartoon cartoon.cpp lines 15-25.
    """
    profile = np.zeros((n, 3))
    for i in range(n):
        t = float(i) / float(n)
        a = t * 2.0 * np.pi + np.pi / 4.0
        profile[i, 0] = np.cos(a) * w / 2.0
        profile[i, 1] = np.sin(a) * h / 2.0
    return profile


def rectangle_profile(n, w, h):
    """
    Rectangular cross-section for beta sheets.
    Matches CPPCartoon cartoon.cpp lines 28-60.
    """
    hw = w / 2.0
    hh = h / 2.0
    segments = [
        (np.array([ hw,  hh, 0.0]), np.array([-hw,  hh, 0.0])),
        (np.array([-hw,  hh, 0.0]), np.array([-hw, -hh, 0.0])),
        (np.array([-hw, -hh, 0.0]), np.array([ hw, -hh, 0.0])),
        (np.array([ hw, -hh, 0.0]), np.array([ hw,  hh, 0.0])),
    ]
    m = n // 4
    profile = np.zeros((n, 3))
    cpt = 0
    for s_start, s_end in segments:
        for i in range(m):
            t = float(i) / float(m)
            profile[cpt] = s_start * (1 - t) + s_end * t
            cpt += 1
    return profile


def rounded_rectangle_profile(n, w, h):
    """
    Rounded rectangle cross-section for helices.
    Matches CPPCartoon cartoon.cpp lines 65-116 exactly.
    Key: r = h/2 (the corner radius equals half the height).
    """
    r = h / 2.0
    hw = w / 2.0 - r
    hh = h / 2.0

    # Segments: top flat, left arc, bottom flat, right arc
    segments_flat = [
        (np.array([ hw,  hh, 0.0]), np.array([-hw,  hh, 0.0])),  # top
        (np.array([-hw, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])),  # left center (used for arc center)
        (np.array([-hw, -hh, 0.0]), np.array([ hw, -hh, 0.0])),  # bottom
        (np.array([ hw, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])),  # right center (used for arc center)
    ]

    m = n // 4
    profile = np.zeros((n, 3))
    cpt = 0

    for s in range(4):
        for i in range(m):
            t = float(i) / float(m)
            if s == 0 or s == 2:
                # Flat segments: linear interpolation
                p = segments_flat[s][0] * (1 - t) + segments_flat[s][1] * t
            elif s == 1:
                # Left semicircle arc
                a = np.pi / 2.0 + np.pi * t
                x = np.cos(a) * r
                y = np.sin(a) * r
                p = segments_flat[s][0] + np.array([x, y, 0.0])
            elif s == 3:
                # Right semicircle arc
                a = 3.0 * np.pi / 2.0 + np.pi * t
                x = np.cos(a) * r
                y = np.sin(a) * r
                p = segments_flat[s][0] + np.array([x, y, 0.0])
            profile[cpt] = p
            cpt += 1

    return profile


def _translate_profile(profile, dx, dy):
    """Translate a profile by (dx, dy). Matches CPPCartoon translateProfile."""
    result = profile.copy()
    result[:, 0] += dx
    result[:, 1] += dy
    return result


def segment_profiles(pp1, pp2, n):
    """
    Determine the start and end profiles for a segment.
    Matches CPPCartoon cartoon.cpp segmentProfiles() lines 130-185.
    
    Args:
        pp1: PeptidePlane at start of segment (pp2 in CPPCartoon's createSegmentMesh)
        pp2: PeptidePlane at end of segment   (pp3 in CPPCartoon's createSegmentMesh)
        n: profile detail (number of profile points)
    
    Returns:
        (profile1, profile2) — both np.ndarray shape (n, 3)
    """
    type0 = pp1.ss1  # SS of the residue BEFORE this segment
    type1, type2 = transition(pp1)

    offset1 = RIBBON_OFFSET
    offset2 = RIBBON_OFFSET

    if pp1.flipped:
        offset1 = -offset1
    if pp2.flipped:
        offset2 = -offset2

    # Profile 1 (start of segment)
    if type1 == HELIX:
        if type0 == STRAND:
            p1 = rounded_rectangle_profile(n, 0.0, 0.0)  # Collapsed
        else:
            p1 = rounded_rectangle_profile(n, RIBBON_WIDTH, RIBBON_HEIGHT)
        p1 = _translate_profile(p1, 0.0, offset1)
    elif type1 == STRAND:
        _, t2_check = transition(pp1)
        if t2_check == STRAND:
            p1 = rectangle_profile(n, ARROW_WIDTH, ARROW_HEIGHT)
        else:
            p1 = rectangle_profile(n, ARROW_HEAD_WIDTH, ARROW_HEIGHT)  # Arrow head
    else:  # Coil
        if type0 == STRAND:
            p1 = ellipse_profile(n, 0.0, 0.0)  # Collapsed at strand end
        else:
            p1 = ellipse_profile(n, TUBE_SIZE, TUBE_SIZE)

    # Profile 2 (end of segment)
    if type2 == HELIX:
        p2 = rounded_rectangle_profile(n, RIBBON_WIDTH, RIBBON_HEIGHT)
        p2 = _translate_profile(p2, 0.0, offset2)
    elif type2 == STRAND:
        p2 = rectangle_profile(n, ARROW_WIDTH, ARROW_HEIGHT)
    else:  # Coil
        p2 = ellipse_profile(n, TUBE_SIZE, TUBE_SIZE)

    # At strand end (arrow tip closes to zero width)
    if type1 == STRAND and type2 != STRAND:
        p2 = rectangle_profile(n, 0.0, ARROW_HEIGHT)

    return p1, p2


def hex_to_rgb_array(hex_str: str) -> np.ndarray:
    """Convert hex string (e.g. '#ff0000') to float rgb array."""
    hex_str = str(hex_str).lstrip('#')
    if len(hex_str) == 6:
        return np.array([int(hex_str[0:2], 16) / 255.0,
                         int(hex_str[2:4], 16) / 255.0,
                         int(hex_str[4:6], 16) / 255.0])
    return np.array([0.2, 0.8, 0.2])

def segment_colors(pp, theme_colors=None):
    """
    Determine the start and end colors for a segment.
    Uses dynamic colors from src.shared.ui.theme if theme_colors not provided.
    """
    if theme_colors is None:
        from src.shared.ui.theme import COLORS
        theme_colors = COLORS
    
    type1, type2 = transition(pp)

    color_map = {
        HELIX:  hex_to_rgb_array(theme_colors.get('ss_helix', '#dc3232')),
        STRAND: hex_to_rgb_array(theme_colors.get('ss_sheet', '#3296dc')),
        COIL:   hex_to_rgb_array(theme_colors.get('ss_coil', '#b4b4b4')),
    }
    
    c1 = color_map.get(type1, color_map.get(COIL)).copy()
    c2 = color_map.get(type2, color_map.get(COIL)).copy()

    # Strand keeps same color through the whole arrow
    if type1 == STRAND:
        c2 = c1.copy()

    return c1, c2
