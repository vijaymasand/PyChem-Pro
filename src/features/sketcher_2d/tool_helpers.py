# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from math import sqrt, sin, cos
from math import pi as PI

from . import geometry as geo
from .app_data import Settings

def get_objs_with_all_children(objs):
    stack = list(objs)
    result = set()
    while len(stack):
        obj = stack.pop()
        result.add(obj)
        stack += obj.children
    return list(result)

def draw_objs_recursively(objs):
    objs = get_objs_with_all_children(objs)
    objs = sorted(objs, key=lambda x : x.redraw_priority)
    for o in objs: o.draw()

def move_objs(objs, dx, dy):
    tr = geo.Transform()
    tr.translate(dx, dy)
    objs = get_objs_with_all_children(objs)
    for o in objs: o.transform(tr)

def find_least_crowded_place_around_atom(atom, distance=10):
    neighbors = atom.neighbors
    if not neighbors:
        return atom.x + distance, atom.y
    angles = [geo.line_get_angle_from_east([atom.x, atom.y, n.x, n.y]) for n in neighbors]
    angles.append(2 * PI + min(angles))
    angles.sort(reverse=True)
    diffs = [angles[i] - angles[i + 1] for i in range(len(angles) - 1)]
    i = diffs.index(max(diffs))
    angle = (angles[i] + angles[i + 1]) / 2
    return atom.x + distance * cos(angle), atom.y + distance * sin(angle)

def calc_average_bond_length(bonds):
    if not bonds: return Settings.bond_length
    lens = [sqrt((b.atom1.x - b.atom2.x)**2 + (b.atom1.y - b.atom2.y)**2) for b in bonds]
    lens.sort()
    return lens[len(lens) // 2]
