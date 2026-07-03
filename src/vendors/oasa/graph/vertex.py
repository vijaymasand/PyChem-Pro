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

"""implements the vertex class"""

#import basic

import copy


class vertex( object):
  """simple vertex class, normaly would not be needed but it can speed up many analytical tasks
  to store data directly in vertex and not get them from the graph connectivity matrix.
  vertex has a value attribute used to store arbitrary object"""

  attrs_to_copy = ("value",)

  def __init__( self):
    self.properties_ = {} # used to store intermediate properties such as distances etc.
    self.value = None  # used to store any object associated with the vertex
    self._neighbors = {} # set of all neighbors
    # Stable, per-graph ordinal assigned by graph.add_vertex.  It is used as
    # the hash so that Python ``set``/``frozenset`` operations over vertices
    # (ring perception, path sorting, coordinate generation) iterate in the
    # SAME order on every run.  Without it, the default id()-based hash made
    # the 2D layout non-deterministic (a different, sometimes overlapping,
    # depiction every time the same molecule was drawn).
    self._ord = None
    self._clean_cache()


  def __hash__( self):
    o = self._ord
    return o if o is not None else object.__hash__( self)


  def __str__( self):
    return ("vertex, value=%s, degree=%d, " % (str( self.value), self.get_degree()) )+str(self.properties_)

  def _clean_cache( self):
    self._cache = {}

  def copy( self):
    other = self.__class__()
    for attr in self.attrs_to_copy:
      setattr( other, attr, copy.copy( getattr( self, attr)))
    return other
  
  def add_neighbor( self, v, e):
    """adds a neighbor connected via e"""
    self._clean_cache()
    self._neighbors[ e] = v

  def remove_neighbor( self, v):
    self._clean_cache()
    to_del = None
    for k, vv in self._neighbors.items():
      if v == vv:
        to_del = k
        break
    if to_del:
      del self._neighbors[ to_del]
    else:
      raise Exception("cannot remove non-existing neighbor")

  def remove_edge_and_neighbor( self, e):
    self._clean_cache()
    if e in self._neighbors.keys():
      del self._neighbors[ e]
    else:
      raise Exception("cannot remove non-existing edge", e)


  def get_neighbors( self):
    return [v for (e,v) in self._neighbors.items() if not e.disconnected]
    #for i in self._neighbors.values():
    #  yield i


  def get_neighbor_connected_via( self, e):
    return self._neighbors[ e]

  def get_edge_leading_to( self, a):
    for b, at in self._neighbors.items():
      if a == at:
        return b
    return None

  def get_degree( self):
    return len( self.neighbors)


  def get_neighbors_with_distance( self, d):
    ret = []
    for v in self.neighbors:
      if 'd' in v.properties_ and v.properties_['d'] == d:
        ret.append( v)
    return ret


  def get_neighbor_edge_pairs( self):
    for e,v in self._neighbors.items():
      if not e.disconnected:
        yield e,v

  def get_neighbor_edges( self):
    return [e for e in self._neighbors.keys() if not e.disconnected]



  # PROPERTIES

  neighbors = property( get_neighbors, None, None, "the neighboring vertices")
  neighbor_edges = property( get_neighbor_edges, None, None, "the neighboring edges")
  degree = property( get_degree, None, None, "the degree of the vertex")
  
