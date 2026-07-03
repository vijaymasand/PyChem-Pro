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

from .edge import edge
from .vertex import vertex
import warnings
import copy
pass  # types wildcard removed
import operator
import time
from functools import reduce


class graph( object):
  """provides a minimalistic graph implementation suitable for analysis of chemical problems,
  even if some care was taken to make the graph work with nonsimple graphs, there are cases where it won't!"""

  uses_cache = True


  def __init__( self, vertices = []):
    if vertices:
      self.vertices = vertices
    else:
      self.vertices = []
    self.edges = set()
    self.disconnected_edges = set()
    self._cache = {}



  def __str__( self):
    str = "graph G(V,E), |V|=%d, |E|=%d" % ( len( self.vertices), len( self.edges))
    return str



  def copy( self):
    """provides a really shallow copy, the vertex and edge objects will remain the same,
    only the graph itself is different"""
    c = self.create_graph()
    c.vertices = copy.copy( self.vertices)
    for e in self.edges:
      i, j = e.get_vertices()
      c.add_edge( i, j, e)
    return c

  def deep_copy( self):
    """provides a deep copy of the graph. The result is an isomorphic graph,
    all the used objects are different"""
    c = self.create_graph()
    for v in self.vertices:
      new = v.copy()
      c.add_vertex( new)
    for e in self.edges:
      v1, v2 = e.get_vertices()
      i1 = self.vertices.index( v1)
      i2 = self.vertices.index( v2)
      new_e = e.copy()
      c.add_edge( c.vertices[i1], c.vertices[i2], new_e)
    return c
    
  def create_vertex( self):
    return vertex()

  def create_edge( self):
    return edge()

  def create_graph( self):
    return self.__class__()
  


  ## MODIFICATION METHODS


  def delete_vertex( self, v):
    self.vertices.remove( v)
    self._flush_cache()


  def add_vertex( self, v=None):
    """adds a vertex to a graph, if v argument is not given creates a new one.
    returns None if vertex is already present or the vertex instance if successful"""
    if not v:
      v = self.create_vertex()
    if v not in self.vertices:
      self.vertices.append( v)
      # Stable ordinal (used as the hash → deterministic set order).  Assign it
      # only ONCE per object: OASA re-adds the SAME vertices to induced
      # subgraphs, and mutating the hash of an object that already lives in
      # another graph's sets would corrupt them.
      if v._ord is None:
        v._ord = len( self.vertices) - 1
    else:
      warnings.warn( "Added vertex is already present in graph %s" % str( v), UserWarning, 2)
      return None
    self._flush_cache()
    return v



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
    # Assign the stable ordinal BEFORE the edge is hashed anywhere (added to
    # self.edges — a set — and used as a neighbor-dict key below), else the set
    # would store it under the old id()-hash and later lookups with the new
    # hash would miss it.  Assign only ONCE: the same edge is re-added to
    # induced subgraphs, and changing its hash mid-flight corrupts the sets
    # that already hold it.
    if e._ord is None:
      e._ord = len( self.edges)
    e.set_vertices( (v1,v2))
    self.edges.add( e)
    v1.add_neighbor( v2, e)
    v2.add_neighbor( v1, e)
    self._flush_cache()
    return e



  def insert_a_graph( self, gr):
    """inserts all edges and vertices to the graph"""
    self.vertices.extend( gr.vertices)
    self.edges.update( gr.edges)
    self._flush_cache()
    

  # do we need it?
  def connect_a_graph( self, gr, v1, v2, e=None):
    """gr is a graph, v1 is vertex in self, v2 is vertex in gr, bond is what to use for connection"""
    self.insert_a_graph( gr)
    self.add_edge( v1, v2, e=e)
    


  def disconnect( self, v1, v2):
    """disconnects vertices v1 and v2, on success returns the edge"""
    if v1 != None and v2 != None:
      e = self.get_edge_between( v1, v2)
      if e:
        self.edges.remove( e)
        v1.remove_neighbor( v2)
        v2.remove_neighbor( v1)
      self._flush_cache()
      return e
    else:
      return None



  def disconnect_edge( self, e):
    self.edges.remove( e)
    v1, v2 = e.get_vertices()
    v1.remove_edge_and_neighbor( e)
    if v1 is not v2:
      v2.remove_edge_and_neighbor( e)
    self._flush_cache()
    



  def remove_vertex( self, v):
    for neigh in v.get_neighbors():
      self.disconnect( v, neigh)
    self.delete_vertex( v)
      


  def get_edge_between( self, v1, v2):
    """takes two vertices"""
    for e in v1.get_neighbor_edges():
      if e in v2.get_neighbor_edges():
        return e
    return None


  def temporarily_disconnect_edge( self, e):
    self.edges.remove( e)
    self.disconnected_edges.add( e)
    e.disconnected = True
    self._flush_cache()
    return e



  def reconnect_temporarily_disconnected_edge( self, e):
    assert e in self.disconnected_edges
    self.disconnected_edges.remove( e)
    self.edges.add( e)
    e.disconnected = False
    self._flush_cache()
    

  def reconnect_temporarily_disconnected_edges( self):
    while self.disconnected_edges:
      e = self.disconnected_edges.pop()
      e.disconnected = False
      self.edges.add( e)
    self._flush_cache()


  ## PROPERTIES METHODS
  ## BOOLEAN
    
  def is_connected( self):
    if len( self.edges) < len( self.vertices) -1:
      # in this case it cannot be connected
      return False
    i = 0
    for x in self.get_connected_components():
      i += 1
      if i > 1:
        return False
    return True
  #return len( self.get_connected_components()) == 1

  def is_tree( self):
    if self.is_connected():
      return len( self.vertices)-1 == len( self.edges)
    else:
      return len( self.vertices)-len( list( self.get_connected_components())) == len( self.edges)

  def is_cycle( self):
    for i in self.get_degrees():
      if i != 2:
        return False
    return True

  def is_euler( self):
    for i in self.get_degrees():
      if i % 2:
        return False
    return True

  def contains_cycle( self):
    """this assumes that the graph is connected"""
    return not self.is_tree()


  def is_edge_a_bridge( self, e):
    """tells whether an edge is bridge"""
    start = list( e.vertices)[0]
    # find number of vertices accessible from one of the edge endpoints
    self.mark_vertices_with_distance_from( start)
    c1 = len( [v for v in self.vertices if 'd' in v.properties_])
    # disconnect the eddge
    self.temporarily_disconnect_edge( e)
    # find the number of vertices accessible now
    self.mark_vertices_with_distance_from( start)
    c2 = len( [v for v in self.vertices if 'd' in v.properties_])
    # if they differ, we've got a bridge
    if c1 > c2:
      x = 1
    else:
      x = 0
    self.reconnect_temporarily_disconnected_edge( e)
    return x



  def is_edge_a_bridge_fast_and_dangerous( self, e):
    """should be used only in case of repetitive questions for the same edge in cases
    where no edges are added to the graph between the questions (if brigde==1 the value
    is stored and returned, which is safe only in case no edges are added)"""
    try:
      return e.properties_['bridge']
    except:
      if self.is_edge_a_bridge( e):
        e.properties_['bridge'] = 1
        return 1
      else:
        return 0


  def get_pieces_after_edge_removal( self, e):
    self.temporarily_disconnect_edge( e)
    ps = [i for i in self.get_connected_components()]
    self.reconnect_temporarily_disconnected_edge( e)
    return ps


  def get_size_of_pieces_after_edge_removal( self, e):
    v1, v2 = e.vertices
    self.disconnect( v1, v2)
    ps = [i for i in self.get_connected_components()]
    self.add_edge( v1, v2, e=e)
    return map( len, ps)
    


  ## ANALYSIS

  def get_connected_components( self):
    """returns the connected components of graph in a form o list of lists of vertices"""
    comp = set() # just processed component 
    comps = []
    not_processed = set( self.vertices)
    if not_processed:
      recent = set() # [not_processed.pop()])
    while not_processed:
      recent = set( reduce( operator.add, [a.get_neighbors() for a in recent], [])) & not_processed
      if not recent:
        if comp:
          yield comp
        recent = set( [not_processed.pop()])
        comp = recent
      else:
        comp.update( recent)
        not_processed -= recent
    # when there is only one atom in the last piece it is not yielded in the loop
    yield comp



  def get_disconnected_subgraphs( self):
    """returns the subgraphs of self, it is dangerous as it reuses the original vertices and
    edges, therefore it should be used only when the old self is no longer needed."""
    vss = self.get_connected_components()
    out = []
    for vs in vss:
      out.append( self.get_induced_subgraph_from_vertices( vs))
    return out


  def get_induced_subgraph_from_vertices( self, vs):
    """it creates a new graph, however uses the old vertices and edges!"""
    g = self.create_graph()
    for v in vs:
      g.add_vertex( v)
    for e in self.vertex_subgraph_to_edge_subgraph( vs):
      v1, v2 = e.get_vertices()
      if v1 in vs and v2 in vs:
        g.add_edge( v1, v2, e)  # BUG - it should copy the edge?
    return g


  def get_induced_copy_subgraph_from_vertices_and_edges( self, vertices, edges, add_back_links=False):
    """it creates a new graph, populates it with copies of vertices and edges;
    it only uses the vertices and edges that are supplied as argument;
    add_back_links - add x.properties_['original'] link to the original edges and vertices
    """
    c = self.create_graph()
    old_v_to_new_v = {}
    for v in vertices:
      new = v.copy()
      c.add_vertex( new)
      old_v_to_new_v[v] = new
      if add_back_links:
        new.properties_['original'] = v
    for e in edges:
      v1, v2 = e.get_vertices()
      if (v1 in old_v_to_new_v) and (v2 in old_v_to_new_v):
        # exclude edges to not-copied vertices (this prevents programmers errors and adds the possibility
        # to replace edges in this call by all edges)
        new_e = e.copy()
        if add_back_links:
          new_e.properties_['original'] = e
        c.add_edge( old_v_to_new_v[v1], old_v_to_new_v[v2], new_e)
    return c

  
      
  def get_degrees( self):
    """returns a generator of degrees, this is useful because for many properties
    the whole list is not important"""
    for v in self.vertices:
      yield v.get_degree()
  

  def get_neighbors( self, v):
    """Info - available also trough the vertex.get_neighbors()"""
    for i in self.get_neighbors_indexes( v):
      yield self.vertices[ i]

  def get_neighbors_indexes( self, v):
    i1 = self._get_vertex_index( v)
    for i, x in enumerate( self.vertices[ i1]):
      if x and i!=i1:
        yield i


  def get_smallest_independent_cycles( self):
    """returns a set of smallest possible independent cycles,
    other cycles in graph are guaranteed to be combinations of them"""
    return list(map( self.edge_subgraph_to_vertex_subgraph, self.get_smallest_independent_cycles_e()))
                

  def get_smallest_independent_cycles_dangerous_and_cached( self):
    try:
      return self._cache['cycles']
    except KeyError:
      self._cache['cycles'] = self.get_smallest_independent_cycles()
      return self._cache['cycles']
    



  def get_almost_all_cycles_e( self):
    """returns almost all cycles found in the graph as sets of edges
    this version is not perfect as it sometimes forgets a few rings"""
    all_cycles = []
    self.mark_vertices_with_distance_from( self.vertices[0])
    to_go = set()
    for ps in self._get_all_ring_end_points():
      to_go.update( ps)
    for ps in self._get_all_ring_start_points():
      to_go.update( ps)
    while to_go:
      v = to_go.pop()
      cycles = self._get_cycles_for_vertex( v, to_reach=v)
      all_cycles += cycles
    all_cycles = set( map( frozenset, all_cycles))
    self.clean_distance_from_vertices()
    return all_cycles



  def get_all_cycles_e_old( self):
    """returns all cycles found in the graph as sets of edges;
    this version of the algorithm strips all non-cyclic (bridge) edges
    and then searches for cycles in the rest"""
    if not self.contains_cycle():
      return set()
    self.temporarily_strip_bridge_edges()

    to_go = set( [v for v in self.vertices if v.degree])
    all_cycles = []
    while to_go:
      v = to_go.pop()
      cycles = self._get_cycles_for_vertex( v, to_reach=v)
      all_cycles += cycles
      to_go -= set( [ver for ver in reduce( operator.or_, map( self.edge_subgraph_to_vertex_subgraph, cycles), set()) if ver.degree == 2])
    all_cycles = set( map( frozenset, all_cycles))

    self.reconnect_temporarily_disconnected_edges()
    return all_cycles



  def get_all_cycles_e_oldest( self):
    """returns all cycles found in the graph as sets of edges"""
    to_go = set( self.vertices)
    for v in self.vertices:
      if v.degree == 1:
        to_go.remove( v)
    all_cycles = []
    while to_go:
      v = to_go.pop()
      cycles = self._get_cycles_for_vertex( v, to_reach=v)
      all_cycles += cycles
      to_go -= set( [ver for ver in reduce( operator.or_, map( self.edge_subgraph_to_vertex_subgraph, cycles), set()) if ver.degree == 2])
    all_cycles = set( map( frozenset, all_cycles))
    return all_cycles




  def _get_cycles_for_vertex( self, v, to_reach=None, processed=None):
    if not processed:
      processed = set()
    for e, neigh in v.get_neighbor_edge_pairs():
      if neigh == to_reach and len( processed) > 1:
        return [set( [e])]
    all_cycles = []
    for e, neigh in v.get_neighbor_edge_pairs():
      if neigh not in processed:
        p = set([v]) | processed
        cycles = self._get_cycles_for_vertex( neigh, to_reach=to_reach, processed=p)
        if cycles:
          for cycle in cycles:
            cycle.add( e)
        all_cycles += cycles
    return all_cycles





  def get_smallest_independent_cycles_e( self):
    """returns a set of smallest possible independent cycles as list of Sets of edges,
    other cycles in graph are guaranteed to be combinations of them.
    Gasteiger J. (Editor), Engel T. (Editor), Chemoinformatics : A Textbook, John Wiley & Sons 2001,
    ISBN 3527306811, 174."""
    #assert self.is_connected()
    ncycles = len( self.edges) - len( self.vertices) + len( list( self.get_connected_components()))
    # check if the graph is connected, don't know if we should do it...
    if ncycles < 0:
      warnings.warn( "The number of edges is smaller than number of vertices-1, the molecule must be disconnected, which means there is something wrong with it.", UserWarning, 3)
      ncycles = 0

    # trivial case of tree graph
    if ncycles == 0:
      return set()

    # the code itself
    self.temporarily_strip_bridge_edges()
    cycles = set()

    vs = [v for v in self.vertices if v.degree]
    while vs and len( cycles) < ncycles:
      new_cycles = set()
      vs2 = [v for v in vs if v.degree == 2]
      # disconnect something if there are no vertices of degree 2
      removed_e = None
      if not vs2:
        for v in vs:
          if v.degree == 3:
            removed_e = list( v.neighbor_edges)[0]
            self.temporarily_disconnect_edge( removed_e)
            break
      vs2 = [v for v in vs if v.degree == 2]
      assert len( vs2) > 0
      # get rings for all degree==2 vertices
      for v in vs2:
        gen = self._get_smallest_cycles_for_vertex( v, to_reach=v)
        for x in gen:
          if x:
            new_cycles.update( set( x))
            break
      if removed_e:
        # we removed an edge - we need to check what cycles it would influence
        # we can also assume that there are only two vertices with degree 2
        # after the removal of this edge and these are the end vertices
        # therefore the code to detect longest path of degree 2 vertices is
        # superfluous
        to_disconnect = [list( removed_e.vertices[0].neighbor_edges)[0]]
        # reconnect the edge removed on the top
        if removed_e:
          self.reconnect_temporarily_disconnected_edge( removed_e)
      else:
        # strip the cycles
        to_disconnect = set()
        for cycle in new_cycles:
          # find the longest degree==2 chain in each cycle
          paths = set()
          to_go = set( [v for v in self.edge_subgraph_to_vertex_subgraph( cycle) if v.degree == 2])
          while to_go:
            now = set( [to_go.pop()])            
            path = set( now)
            while now:
              now = reduce( operator.or_, [set( [n for n in v.neighbors if n.degree == 2]) for v in now])
              now &= to_go
              to_go -= now
              path.update( now)
            if path:
              paths.add( frozenset( path))
          l = max( map( len, paths))
          path = [p for p in paths if len( p) == l][0]
          # now mark them for disconnection
          v1 = set( path).pop()
          to_disconnect.add( list( v1.neighbor_edges)[0])
      # disconnect what needs to be disconnected
      [self.temporarily_disconnect_edge( e) for e in to_disconnect]

      # add new_cycles to cycles
      cycles.update( new_cycles)

      # strip the degree==1 vertices
      vs1 = [v for v in self.vertices if v.degree == 1]
      while vs1:
        for v in vs1:
          # we have to ask the degree, because the bond might have been stripped in this run
          if v.degree:
            e = v.neighbor_edges[0]
            self.temporarily_disconnect_edge( e)
        vs1 = [v for v in self.vertices if v.degree == 1]

      self.temporarily_strip_bridge_edges()
      vs = [v for v in self.vertices if v.degree]

    # remove extra cycles in some cases like adamantane
    if len( cycles) - ncycles > 0:
      # sort cycles according to length
      cs = [(len( c), c) for c in cycles]
      cs.sort()
      cs = [c[1] for c in cs]
      # now try to remove the biggest ones
      while len( cs) - ncycles > 0:
        c = set( cs.pop( -1))
        if not c <= reduce( operator.or_, map( set, cs)):
          cs.insert( 0, c)
      cycles = set( [frozenset( _c) for _c in cs])

    # count the cycles and report warnings if their number is wrong
    if len( cycles) < ncycles:
      warnings.warn( "The number of cycles found (%d) is smaller than the theoretical value %d (|E|-|V|+1)" % (len( cycles), ncycles), UserWarning, 2)
    elif len( cycles) > ncycles: 
      warnings.warn( "The number of independent cycles found (%d) is larger than the theoretical value %d (|E|-|V|+1), but I cannot improve it." % (len( cycles), ncycles), UserWarning, 2)

    self.reconnect_temporarily_disconnected_edges()
      
    return cycles



  def _get_smallest_cycle_for_vertex( self, v, to_reach=None, came_from=None, went_through=None):
    """ingenious generator-based breadth-first search (BFS) to find one (random) smallest
    cycle for given vertex. It yields None or cycle for each depth level"""
    for e, neigh in v.get_neighbor_edge_pairs():
      if neigh == to_reach and e != came_from:
        yield set( [came_from, e])
    gens = []
    yield None
    w = went_through and went_through+[v] or [v]
    for e, neigh in v.get_neighbor_edge_pairs():
      if (not went_through or neigh not in went_through) and not e == came_from:
        gens.append( self._get_smallest_cycle_for_vertex( neigh, to_reach=to_reach, came_from=e, went_through=w))
    while 1:
      for i, gen in enumerate( gens):
        ret = next(gen)
        if ret:
          if came_from:
            yield set( [came_from]) | ret
          else:
            yield ret
      yield None


  def _get_smallest_cycles_for_vertex( self, v, to_reach=None, came_from=None, went_through=None):
    """ingenious generator-based breadth-first search (BFS) to find smallest
    cycles for given vertex. It yields None or cycles for each depth level"""
    ret = []
    for e, neigh in v.get_neighbor_edge_pairs():
      if neigh == to_reach and e != came_from:
        ret.append( frozenset( [came_from, e]))
    yield ret

    gens = []
    w = went_through and went_through+[v] or [v]
    for e, neigh in v.get_neighbor_edge_pairs():
      # we dont want to go back, therefore we use went_through
      if (not went_through or neigh not in went_through) and not e == came_from:
        gens.append( self._get_smallest_cycles_for_vertex( neigh, to_reach=to_reach, came_from=e, went_through=w))

    counter = 0
    while 1:
      counter += 1
      if counter > 10000:
        raise Exception( "This exception is triggered when ring-searching code does not find anything in 10000 steps. This means that there is with high probability something wrong with the code and it would be futile to continue in an endless search. Please let me know at beda@zirael.org where this happened")

      all_rets = []
      for gen in gens:
        rets = next(gen)
        new_rets = []
        if rets:
          for ret in rets:
            if came_from:
              ret = set( ret)
              ret.update( set( [came_from]))
              ret = frozenset( ret)
            new_rets.append( ret)
          all_rets.extend( frozenset( new_rets))
      yield all_rets



  def get_all_cycles_old( self):
    """returns all cycles found in the graph as sets of vertices,
    use get_all_cycles_e to get the edge variant, which is better for building new
    graphs as the mapping edges => vertices is unambiguous, while edges=>vertices=>edges might
    include some more edges"""
    return map( self.edge_subgraph_to_vertex_subgraph, self.get_all_cycles_e_old())


  def get_all_cycles_e( self):
    return map( self.vertex_subgraph_to_edge_subgraph, self.get_all_cycles())
  

  def get_all_cycles( self):
    """
    implementation of:
    A New Algorithm for Exhaustive Ring Perception in a Molecular Graph
    Th. Hanser, Ph. Jauffret, and G. Kaufmann
    J. Chem. Inf. Comput. Sci., 1996, 36 (6), 1146-1152 . DOI: 10.1021/ci960322f
    """
    pgraph = self._get_p_graph()
    rings = set()
    for pv in copy.copy( pgraph.vertices):
      rings |= graph._p_graph_remove( pv, pgraph)
    final_rings = set()
    for ring in rings:
      final_ring = frozenset( [v.properties_['original'] for v in ring])
      final_rings.add( final_ring)
    return final_rings
  

  def _get_p_graph( self):
    """helper method for p-graph (path graph) generation"""
    p = self.deep_copy()
    p.temporarily_strip_bridge_edges()
    # we count on order of vertices remaining the same
    for i,v in enumerate( self.vertices):
      p.vertices[i].properties_['original'] = v
    for e in p.edges:
      e.path_ = set( e.vertices)
    to_remove = [v for v in p.vertices if not v.neighbors]
    for v in to_remove:
      p.delete_vertex( v)
    return p

  @staticmethod
  def _p_graph_remove( v, pgraph):
    rings = set()
    neighbor_edge_vertex_pairs = list( v.get_neighbor_edge_pairs())
    new_edges = []
    for i,(ne1,nv1) in enumerate( neighbor_edge_vertex_pairs):
      for ne2,nv2 in neighbor_edge_vertex_pairs[i+1:]:
        if (nv1 is nv2 and (ne1.path_ & ne2.path_ == set( [v,nv2]))) or (ne1.path_ & ne2.path_ == set( [v])):
          new_path = ne1.path_ | ne2.path_
          new_edge = pgraph.create_edge()
          new_edge.path_ = new_path
          pgraph.add_edge( nv1, nv2, new_edge)
          new_edges.append( new_edge)
    for ne,nv in neighbor_edge_vertex_pairs:
      pgraph.disconnect_edge( ne)
    for e in new_edges:
      end1, end2 = e.vertices
      if end1 is end2:
        # we have a loop here
        rings.add( frozenset( e.path_))
        pgraph.disconnect_edge( e)
    pgraph.remove_vertex( v)
    return rings

  def mark_vertices_with_distance_from( self, v):
    """returns the maximum d"""
    self.clean_distance_from_vertices()
    return self._mark_vertices_with_distance_from( v)

  def clean_distance_from_vertices( self):
    for i in self.vertices:
      try:
        del i.properties_['d']
      except KeyError:
        pass

  def mark_edges_with_distance_from( self, e1):
    for e in self.edges:
      try:
        del e.properties_['d']
      except KeyError:
        pass
    marked = set( [e1])
    new = set( [e1])
    dist = 0
    e1.properties_['dist'] = dist
    while new:
      new_new = set()
      dist += 1
      for e in new:
        for ne in e.get_neighbor_edges():
          if not ne in marked:
            ne.properties_['dist'] = dist
            new_new.add( ne)
      new = new_new
      marked.update( new)

  def get_path_between_edges( self, e1, e2):
    self.mark_edges_with_distance_from( e1)
    if not "dist" in e2.properties_:
      return None
    else:
      path = [e2]
      _e = e2
      for i in range( e2.properties_['dist']-1, -1, -1):
        _e = [ee for ee in _e.get_neighbor_edges() if ee.properties_['dist'] == i][0]
        path.append( _e)
      return path
  

  def _gen_diameter_progress( self):
    """this generator iteratively generates graph diameter during its computation,
    the result is the last value, it is only interesting for monitoring of the computation
    as it can be pretty time consuming"""
    diameter = 0
    best = None
    best_path = None
    for v in self.vertices:
      dist = self.mark_vertices_with_distance_from( v)
      if dist > diameter:
        diameter = dist
        best = v
        #end = [x for x in self.vertices if x.properties_['d'] == dist][0]
        #best_path = get_path_down_to( end, v)
        yield diameter
    if diameter == 0:
      yield 0
    #best_path.reverse()


  def _get_width_from_vertex( self, v):
    """returns width of the graph as calculated from vertex v"""
    d = 0
    to_mark = set([v])
    marked = set()
    while to_mark:
      marked_before = marked
      marked = to_mark
      to_mark_next = set( reduce( operator.add, [i.get_neighbors() for i in to_mark], []))
      to_mark = to_mark_next - marked - marked_before
      d += 1
    return d-1


  def get_diameter( self):
    d = self._get_cache( "diameter")
    if not d:
      g = self._gen_diameter_progress()
      for i in g:
        pass
      self._set_cache( "diameter", i)
      return i
    else:
      return d


# does not work because the mark_vertices_with_distance_from is not thread safe,
# it writes to v.properties_['d']. It would have to be switched to v.properties_['dX']
# where X is the number of the thread
##   def get_diameter_multi_thread( self, processes=2, group_by=100):
##     threads = []
##     attrs = list( self.vertices)
##     diameter = 0
    
##     while attrs or threads:
##       if len( threads) < processes and attrs:
##         if len( attrs) > group_by:
##           _as = attrs[0:group_by]
##           del attrs[0:group_by]
##         else:
##           _as = attrs
##           attrs = []

##         t = MyThread( self.mark_vertices_with_distance_from, _as)
##         threads.append( t)
##         t.start()
##       else:
##         time.sleep( 0.05)

##       dead = [t for t in threads if not t.isAlive()]
##       for d in dead:
##         dist = max( d.ret)
##         if dist > diameter:
##           diameter = dist

##       threads = [t for t in threads if t.isAlive()]

##     return diameter


    


  def vertex_subgraph_to_edge_subgraph( self, cycle):
    ret = set()
    for v1 in cycle:
      for (e,n) in v1.get_neighbor_edge_pairs():
        if n in cycle:
          ret.add( e)
    return ret

  def edge_subgraph_to_vertex_subgraph( self, cycle):
    ret = set()
    for e in cycle:
      v1, v2 = e.get_vertices()
      ret.add( v1)
      ret.add( v2)
    return ret

  def get_new_induced_subgraph( self, vertices, edges):
    """returns a induced subgraph that is newly created and can be therefore freely
    changed without worry about the original."""
    sub = self.create_graph()
    vertex_map = {}
    i = 0
    for v in vertices:
      new_v = v.copy()
      sub.add_vertex( new_v)
      vertex_map[v] = i
      i += 1
    for e in edges:
      new_e = e.copy()
      v1, v2 = e.get_vertices()
      sub.add_edge( vertex_map[v1], vertex_map[v2], new_e)
    return sub

  def defines_connected_subgraph_e( self, edges):
    sub = self.get_new_induced_subgraph( self.edge_subgraph_to_vertex_subgraph( edges), edges)
    return sub.is_connected()

  def defines_connected_subgraph_v( self, vertices):
    sub = self.get_new_induced_subgraph( vertices, self.vertex_subgraph_to_edge_subgraph( vertices))
    return sub.is_connected()
  

  def find_path_between( self, start, end, dont_go_through=[]):
    """finds path between two vertices, if dont_go_through is given (as a list of vertices and edges),
    only paths not containing these vertices will be given (or None is returned if such a path
    does not exist"""
    ### DOES NOT WORK WELL WITH RINGS, FOR THIS RECURSIVE DESIGN WILL BE NEEDED
    self.mark_vertices_with_distance_from( start)
    d = end.properties_['d']
    out = [end]
    rend = end
    for i in range( d, 0, -1):
      vs = [v for e,v in rend.get_neighbor_edge_pairs() if v.properties_['d'] == i-1 and (v not in dont_go_through or e not in dont_go_through)]
      if not vs:
        return None
      v = vs[0]
      out.append( v)
      rend = v
    return out

  def sort_vertices_in_path( self, path, start_from=None):
    """returns None if there is no path"""
    rng = copy.copy( path)
    if start_from:
      a = start_from
      rng.remove( a)
    else:
      a = None
      # for acyclic path we need to find one end
      for at in path:
        if len( [v for v in at.get_neighbors() if v in path]) == 1:
          a = at
          rng.remove( at)
          break
      if not a:
        a = rng.pop() # for rings
    out = [a]
    while rng:
      try:
        a = [i for i in a.get_neighbors() if i in rng][0]
      except IndexError:
        return None
      out.append( a)
      rng.remove( a)
    return out


  def temporarily_strip_bridge_edges( self):
    """strip all edges that are a bridge, thus leaving only the cycles connected"""
    bridge_found = True
    while bridge_found:
      vs = [v for v in self.vertices if v.degree == 1]
      while vs:
        for v in vs:
          # we have to ask the degree, because the bond might have been stripped in this run
          if v.degree:
            e = v.neighbor_edges[0]
            self.temporarily_disconnect_edge( e)
        vs = [v for v in self.vertices if v.degree == 1]
      bridge_found = False
      for e in self.edges:
        if self.is_edge_a_bridge( e):
          bridge_found = True
          break
      if bridge_found:
        self.temporarily_disconnect_edge( e)


  def dump_simple_text_file( self, f):
    """dumps itself into a simple text file, one edge per line described as
    indexes of vertices separated with a space, the first line has the vertex
    count.
    f is an file-like object opened for writing"""
    f.write( "%d\n" % len( self.vertices))
    for v in self.vertices:
      i1 = self.vertices.index( v)
      for n in v.neighbors:
        i2 = self.vertices.index( n)
        f.write( "%d %d\n" % (i1, i2))


  def read_simple_text_file( self, f):
    """fills itself from a simple text file, one edge per line described as
    indexes of vertices separated with a space, the first line has the vertex
    count.
    f is an file-like object opened for reading"""
    num = int( f.readline().strip())
    for i in range( num):
      self.add_vertex( self.create_vertex())
    for line in f:
      i1, i2 = map( int, line.split())
      self.add_edge( self.vertices[i1], self.vertices[i2])

  def path_exists( self, a1, a2):
    self.mark_vertices_with_distance_from( a1)
    if 'd' in a2.properties_:
      return True
    else:
      return False

  ## MAXIMUM MATCHING RELATED STUFF

  def get_initial_matching( self):
    """simple method to get initial matching"""
    nrex = len( self.vertices)
    mate = {}.fromkeys( self.vertices, 0)
    #return mate, nrex
    for v in self.vertices[:-1]:
      if mate[v] == 0:
        for n in v.neighbors:
          if mate[n] == 0:
            mate[n] = v
            mate[v] = n
            nrex -= 2
            break
    return mate, nrex

  def get_initial_matching_old( self):
    """simple method to get initial matching"""
    nrex = len( self.vertices)
    mate = nrex*[0]
    for i,v in enumerate( self.vertices[:-1]):
      if mate[i] == 0:
        for n in v.neighbors:
          j = self.vertices.index( n)
          if mate[j] == 0:
            mate[j] = v
            mate[i] = n
            nrex -= 2
            break
    return mate, nrex

  def find_augmenting_path_from( self, start, mate):
    """tries to find augmenting path from start to any other exposed vertex;
    it is not the most sophisticated algorithm and does not treat blossoms very
    effectively, but it should work.
    """
    def add_copy_vertex( _v):
      _new_v = alt_tree.create_vertex()
      alt_tree.add_vertex( _new_v)
      _new_v.properties_['original'] = _v
      return _new_v
    
    alt_tree = self.create_graph() # alternating tree
    copy_root = add_copy_vertex( start)
    new_layer = [copy_root]
    inner = True # we are now going to process inner vertices
    hit = None
    while new_layer and not hit:
      next_new_layer = []
      for copy_parent in new_layer:
        if inner:
          # to prevent backtracking to already walked path
          path_from_root_to_copy_parent_original = [x.properties_['original'] for x in alt_tree.find_path_between( copy_root, copy_parent)]
          for child in copy_parent.properties_['original'].neighbors:
            if child not in path_from_root_to_copy_parent_original:
              new_v = add_copy_vertex( child)
              alt_tree.add_edge( copy_parent, new_v)
              if mate[child] == 0:
                hit = new_v
                break
              next_new_layer.append( new_v)
        else:
          m = mate[copy_parent.properties_['original']]
          new_v = add_copy_vertex( m)
          alt_tree.add_edge( copy_parent, new_v)
          next_new_layer.append( new_v)
      inner = not inner
      new_layer = next_new_layer
    if hit:
      path = alt_tree.find_path_between( hit, copy_root)
      aug_p = [v.properties_['original'] for v in path]
      return aug_p
    return None
  

  def update_matching_using_augmenting_path( self, path, mate):
    assert len( path) % 2 == 0 # there must be an even number of vertices
    for i in range( 0, len( path), 2):
      v1 = path[i]
      v2 = path[i+1]
      mate[v1] = v2
      mate[v2] = v1
    return mate


  def get_maximum_matching( self):
    mate, nrex = self.get_initial_matching()
    while nrex > 1:
      #print "NREEEEEEX", nrex
      #self._print_mate( mate)
      exposed = [v for v,m in mate.items() if m == 0]
      aug = self.find_augmenting_path_from( exposed[0], mate)
      if not aug:
        break
      mate = self.update_matching_using_augmenting_path( aug, mate)
      nrex -= 2
    return mate, nrex
    
  def _print_mate( self, mate):
    print("MATE", end="")
    for k,v in mate.items():
      if v:
        print("%d-%d" % (self.vertices.index(k), self.vertices.index(v)), end="")
    print("END")

  ## STATIC METHODS

##   def get_random_longest_path_numbered( self, start):
##     """vertices have to be freshly marked with distance"""
##     now = start
##     path = []
##     d = 0
##     while now:
##       d += 1
##       path.append( now)
##       ns = [v for v in now.neighbors if 'd' in v.properties_ and v.properties_['d'] == d]
##       if ns:
##         now = ns[0]
##       else:
##         now = None
##     return path

##   get_random_longest_path_numbered = classmethod( get_random_longest_path_numbered)


##   def get_random_longest_path( self, start):
##     self.mark_vertices_with_distance_from( start)
##     return self.get_random_longest_path_numbered( start)

##   get_random_longest_path = classmethod( get_random_longest_path)




  # PRIVATE METHODS

  def _get_vertex_index( self, v):
    """if v is already an index, return v, otherwise return index of v on None"""
    if type( v) == int and v < len( self.vertices):
      return v
    try:
      return self.vertices.index( v)
    except ValueError:
      return None


  def _read_file( self, name="/home/beda/oasa/oasa/mol.graph"):
    self.vertices = []
    self.edges = set()
    f = file( name, 'r')
    vs = f.readline()
    [self.add_vertex() for i in vs.split(" ")]
    for l in f.readlines():
      c, a, b = l.split(' ')
      self.add_edge( self.vertices[int(a)], self.vertices[int(b)])
    f.close()


  def _mark_vertices_with_distance_from( self, v):
    """returns the maximum d"""
    d = 0
    to_mark = set([v])
    while to_mark:
      to_mark_next = set()
      for i in to_mark:
        i.properties_['d'] = d

      for i in to_mark:
        for j in i.get_neighbors():
          if 'd' not in j.properties_:
            to_mark_next.add( j)

      to_mark = to_mark_next
      d += 1

    return d-1


  def _get_some_cycles( self):
    if len( self.vertices) <= 2:
      raise StopIteration
    self.mark_vertices_with_distance_from( self.vertices[0])
    for end in self._get_all_ring_end_points():
      for start in self._get_all_ring_start_points():
        ring = is_there_a_ring_between( start, end)
        if ring:
          yield ring


  def _get_all_ring_end_points( self):
    already_there = []
    for v in self.vertices:
      if v in already_there:
        continue
      vs_ed, vs_ver = is_ring_end_vertex( v)
      if vs_ed:
        yield vs_ed
      if vs_ver:
        already_there.extend( vs_ver)
        yield vs_ver

  def _get_all_ring_start_points( self):
    already_there = []
    for v in self.vertices:
      if v in already_there:
        continue
      vs_ed, vs_ver = is_ring_start_vertex( v)
      if vs_ed:
        yield vs_ed
      if vs_ver:
        already_there.extend( vs_ver)
        yield vs_ver



  def _flush_cache( self):
    self._cache = {}


  def _set_cache( self, name, value):
    if self.uses_cache:
      self._cache[ name] = value


  def _get_cache( self, name):
    return self._cache.get( name, None)
    





def is_ring_end_vertex( v):
# NEEDS NEW COMMENT
#  """tells if a vertex has two neighbors with 'd' one smaller and equal or
#  one neighbor with equal 'd'. These are the conditions for a cycle end.
#  Returns the set([v]), in second case 'v' and the other with same 'd'"""
  ed, ver = None, None
  for x in v.get_neighbors():
    if v.properties_['d'] == x.properties_['d']:
      if len( v.get_neighbors_with_distance( v.properties_['d']-1)) and len( x.get_neighbors_with_distance( v.properties_['d']-1)):
        ed = set([ x, v])
    for y in v.get_neighbors():
      if x != y:
        if (x.properties_['d'] == y.properties_['d']) and (x.properties_['d'] == v.properties_['d']-1):
          ver = set([v])
  return ed, ver

def is_ring_start_vertex( v):
# NEEDS NEW COMMENT
#  """tells if a vertex has two neighbors with 'd' one higher and equal or one neighbor with
#  equal 'd' (then both have neighbors with one higher 'd'). These are the conditions for a cycle start.
#  Returns boolean"""
  ed, ver = None, None
  for x in v.get_neighbors():
    if v.properties_['d'] == x.properties_['d']:
      if len( v.get_neighbors_with_distance( v.properties_['d']+1)) and len( x.get_neighbors_with_distance( v.properties_['d']+1)):
        ed = set([ x, v])
    for y in v.get_neighbors():
      if x != y:
        if (x.properties_['d'] == y.properties_['d']) and (x.properties_['d'] == v.properties_['d']+1):
          ver = set([v])
          break
  return ed, ver

def get_first_closer_by_one( v):
  d = v.properties_['d']
  for x in v.get_neighbors():
    if x.properties_['d'] == d-1:
      return x
  return None

def is_there_a_ring_between( start, end):
  pths = []
  for e in end:
    for s in start:
      ps = get_paths_down_to( e, s)
      if ps:
        ## for end-edge there is only one path from each point important
        if len( end) == 2:
          ps = [ps[0]]
        ## for end-vertex the end vertex must not be part of the set in order to filter it
        if len( end) == 1:
          for p in ps:
            p.remove( e)
        pths.extend( ps)
  ## we use sets for filtering off the path that have common vertices 
  paths = [set( p) for p in pths]
  final = []
  while len( paths) > 0:
    p1 = paths.pop(0)
    to_del= []
    for p2 in paths:
      if p1 & p2:
        to_del.append( p2)
    for p2 in to_del:
      paths.remove( p2)
    final.append( p1)
  if len( final) >= 2:
    final.append( end | start)
    ret = set( reduce( operator.or_, final))
    return ret
  return False
                
def get_paths_down_to( end, start):
  paths = []
  if end == start:
    return None
  for x in end.get_neighbors():
    if x.properties_['d'] == end.properties_['d']-1:
      ps = get_path_down_to( x, start)
      if ps != None:
        ps.append( end)
        paths.append( ps)
  return paths


def get_path_down_to( end, start):
  if end == start:
    return []
  for x in end.get_neighbors():
    if x.properties_['d'] == end.properties_['d']-1:
      ps = get_path_down_to( x, start)
      if ps != None:
        ps.append( end)
        return ps

  return None





def filter_off_dependent_cycles( cycles):
  """this filtres off all cycles that are combinations of smaller cycles."""
  cs = list( cycles)
  level = 2
  cs.sort( lambda a,b: -len( a) + len(b))
  while level < len( cs):
    to_del = set()
    for c in cs:
      combs = gen_variations( [x for x in cs if x != c and x not in to_del and x & c], level)
      for combl in combs:
        comb = set( combl)
        if not comb & to_del and (c <= reduce( operator.or_, comb)) and reduce( operator.and_, comb):
          to_del.add( c)
          break
    [cs.remove( x) for x in to_del]
    level += 1
  return cs



## def filter_off_dependent_cycles( cycles):
##   cs = list( cycles)
##   level = 2
##   cs.sort( lambda a,b: -len( a) + len(b))
##   print len( cs), map( len, cs)
##   while level < len( cs):
##     to_del = set()
##     for c in cs:
##       combs = gen_variations( [x for x in cs if x != c and x & c], level)
##       for combl in combs:
##         comb = set( combl)
##         if not comb & to_del and (c <= reduce( operator.or_, comb)):
##           to_del.add( c)
##           print "remove", len( c), map( len, comb)
##           break
##     [cs.remove( x) for x in to_del]
##     level += 1
##   return cs




def gen_variations(items, n):
    if n==0:
      yield []
    else:
      for i in range( len(items)-n+1):
        for v in gen_variations(items[i+1:],n-1):
          yield [items[i]]+v




from threading import Thread
import time


class MyThread( Thread):

  def __init__( self, function, attrs):
    """function is the function to run, attrs a list of attrs for which to run the
    consecutively the function, additional_attrs and kw_attrs are attrs that are added to
    each call of function"""
    Thread.__init__( self)
    self.function = function
    self.attrs = attrs
    self.ret = []


  def run( self):
    for a in self.attrs:
      self.ret.append( self.function( a))





## import profile

## g = graph()
## g._read_file()
## print g
## print [c for c in g.get_connected_components()]

## for e in g.edges:
##   print g.is_edge_a_bridge( e),
## print

## print len( [i for i in g.get_connected_components()])
## for e in g.edges:
##   d1, d2 = [x.get_degree() for x in e.get_vertices()]
##   if d1 != 2 or d2 != 2:
##     v1, v2 = e.vertices
##     break
## g.disconnect( v1, v2)
## print len( [i for i in g.get_connected_components()])


## g = g.deep_copy()
## print g
## #profile.run( 'g.get_all_cycles()')
## c = g.get_smallest_independent_cycles()
## print 'cycles %d, lengths %s' % (len( c), str( map( len, c)))
## c = g.get_all_cycles()
## print 'cycles %d, lengths %s' % (len( c), str( map( len, c)))




##################################################
# TODO


##################################################
