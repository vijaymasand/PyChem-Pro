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


"""this module contains a graph class that provides a minimalistic
graph implementation suitable for analysis of chemical problems"""


from .diedge import diedge
from .vertex import vertex
import warnings
import copy
import operator
from . import graph
from functools import reduce



class digraph( graph.graph):
  """provides a minimalistic graph implementation suitable for analysis of chemical problems,
  even if some care was taken to make the graph work with nonsimple graphs, there are cases where it won't!"""

  vertex_class = vertex
  edge_class = diedge


  def __str__( self):
    str = "digraph G(V,E), |V|=%d, |E|=%d" % ( len( self.vertices), len( self.edges))
    return str

  ## MODIFICATION METHODS


  def add_edge( self, v1, v2, e=None):
    """adds an edge to a graph connecting vertices v1 and v2, if e argument is not given creates a new one.
    returns None if operation fails or the edge instance if successful"""
    i1 = self._get_vertex_index( v1)
    i2 = self._get_vertex_index( v2)
    if i1 == None or i2 == None:
      warnings.warn( "Adding edge to a vertex not present in graph failed (of course)", UserWarning, 3)
      return None
    # to get the vertices if v1 and v2 were indexes
    v1 = self.vertices[ i1]
    v2 = self.vertices[ i2]
    if not e:
      e = self.create_edge()
    e.set_vertices( (v1,v2))
    self.edges.add( e)
    v1.add_neighbor( v2, e)
    return e



  def get_diameter( self):
    diameter = 0
    best = None
    best_path = None
    for v in self.vertices:
      dist = self.mark_vertices_with_distance_from( v)
      if dist > diameter:
        diameter = dist
        best = v
        end = [x for x in self.vertices if 'd' in x.properties_ and x.properties_['d'] == dist][0]
        best_path = self.get_random_longest_path_numbered( v, end)

    print("path")
    best_path.reverse()
    for v in best_path:
      print(v)
    return diameter


  def get_connected_components( self):
    """returns the connected components of graph in a form o list of lists of vertices"""
    comp = set() # just processed component 
    comps = []
    not_processed = set( self.vertices)
    if not_processed:
      recent = set() # [not_processed.pop()])
    processed = set()
      
    while not_processed:
      recent = set( reduce( operator.add, [a.get_neighbors() for a in recent], [])) & not_processed
      processed = recent | comp
      for e in self.edges:
        v1, v2 = e.get_vertices()
        if (v1 in processed and v2 in not_processed):
          recent.add( v2)
        elif (v2 in processed and v1 in not_processed):
          recent.add( v1)
      if not recent:
        if comp:
          yield comp
        recent = set( [not_processed.pop()])
        comp = recent
      else:
        comp |= recent
        not_processed -= recent
    # when there is only one atom in the last piece it is not yielded in the loop
    yield comp

  
  def get_random_longest_path_numbered( self, start, end):
    """vertices have to be freshly marked with distance"""
    now = end
    path = []
    d = end.properties_['d']
    while now:
      d -= 1
      path.append( now)
      ns = [v for v in self.vertices if 'd' in v.properties_ and v.properties_['d'] == d and now in v.neighbors]
      if ns:
        now = ns[0]
      else:
        now = None
    return path


  def get_graphviz_text_dump( self):
    ret = '''digraph "dump" {
    fontpath="/usr/share/fonts/corefonts";
    ratio=compress
    nodesep=0.1;
    ranksep=0.3;
    rankdir=LR;
    node [shape=box,fontsize=10,fontname=Arial,height=0.3];
    '''
    for v in self.vertices:
      for n in v.neighbors:
        ret += '"%s" -> "%s";\n' % (v,n);
      ret += '"%s" [label="%s"];\n' % (v,v.value)
    ret += "}"
    return ret
