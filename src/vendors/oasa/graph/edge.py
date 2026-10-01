#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2008 Beda Kosata <beda@zirael.org>

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file gpl.txt in the
#     main directory of the program

#--------------------------------------------------------------------------



from warnings import warn
import copy

class edge:

  attrs_to_copy = ("disconnected",)

  def __init__( self, vs=[]):
    # Stable, per-graph ordinal assigned by graph.add_edge; used as the hash so
    # that ``set``/``frozenset`` operations over edges (OASA stores edges in a
    # set and ring perception pops/iterates them) are deterministic across runs.
    # See graph/vertex.py for the rationale — id()-based hashing made the 2D
    # coordinate layout differ every time the same molecule was drawn.
    self._ord = None
    self.vertices = set([])
    self.set_vertices( vs)
    self.properties_ = {}
    self.disconnected = False

  def __hash__( self):
    o = self._ord
    return o if o is not None else object.__hash__( self)

  def __str__( self):
    return "edge between %s %s" % tuple( map( str, self.vertices))

  def copy( self):
    other = self.__class__()
    for attr in self.attrs_to_copy:
      setattr( other, attr, copy.copy( getattr( self, attr)))
    return other

  def set_vertices( self, vs=[]):
    if vs and len( vs) == 2:
      self.vertices = set( vs)

  def get_vertices( self):
    return self.vertices

  def get_neighbor_edges( self):
    v1, v2 = self.vertices
    out1 = [e for e in v1.get_neighbor_edges() if e!=self]
    out2 = [e for e in v2.get_neighbor_edges() if e!=self]
    return out1 + out2

  def get_neighbor_edges2( self):
    """returns 2 lists of neighbor edges (one for one side, one for the other)"""
    v1, v2 = self.vertices
    out1 = [e for e in v1.get_neighbor_edges() if e!=self]
    out2 = [e for e in v2.get_neighbor_edges() if e!=self]
    return out1, out2
    
  def _get_disconnected( self):
    return self._disconnected

  def _set_disconnected( self, d):
    self._disconnected = d
   # if d:
   #   warn(a "aaaaaa", UserWarning, 4)


  disconnected = property( _get_disconnected, _set_disconnected)
