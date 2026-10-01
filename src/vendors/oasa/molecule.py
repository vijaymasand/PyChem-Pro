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


from . import graph
from .atom import atom
from .query_atom import query_atom
from .bond import bond
import copy
from . import common
import operator
from . import misc
from . import periodic_table as PT
import math
from . import transform3d


class molecule( graph.graph):


  def __init__( self, vertices =[]):
    graph.graph.__init__( self, vertices=vertices)
    # aliases
    self.atoms = self.vertices
    self.bonds = self.edges
    self.stereochemistry = []

  def __str__( self):
    return "molecule, %d atoms, %d bonds" % (len( self.vertices), len( self.edges))

  def create_vertex( self, vertex_class=None):
    if not vertex_class:
      return atom()
    else:
      return vertex_class()

  def create_edge( self):
    return bond()

  def create_graph( self):
    from . import config
    return config.Config.molecule_class()

  def add_stereochemistry( self, stereo):
    self.stereochemistry.append( stereo)

  def remove_stereochemistry( self, stereo):
    if not stereo in self.stereochemistry:
      raise ValueError("cannot remove non-existent stereochemistry information")
    self.stereochemistry.remove( stereo)

  def get_stereochemistry_by_center( self, center):
    for st in self.stereochemistry:
      if st.center == center:
        return st
    return None

  # analytics

  # override of graphs method to add stereochemistry support
  def get_disconnected_subgraphs( self):
    out = graph.graph.get_disconnected_subgraphs( self)
    for part in out:
      for st in self.stereochemistry:
        if set( [ref for ref in st.references if isinstance(ref,atom)]) <= set( part.vertices):
          part.add_stereochemistry( st)
    return out

  def get_mol_weight( self):
    w = 0
    for v in self.vertices:
      w += v.weight
      if v.free_valency > 0:
        w += v.free_valency * PT.periodic_table['H']['weight']
    return w

  weight = property( get_mol_weight, None, None, "molecular weight")


  def get_charge( self):
    return sum( [a.charge for a in self.vertices])

  charge = property( get_charge, None, None, "Net charge of the molecule")


  def get_formula_dict( self):
    """returns a formula dict as defined in the periodic_table.py::formula_dict"""
    comp = PT.formula_dict()
    for a in self.atoms:
      comp += a.get_formula_dict()
    return comp



  def add_missing_hydrogens( self):
    hs = set()
    for v in copy.copy( self.vertices):
      hs |= self.add_missing_hydrogens_to_atom( v)
    return hs


  def add_missing_hydrogens_to_atom( self, v):
    hs = set()
    for i in range( v.free_valency):
      h = self.create_vertex()
      h.symbol = 'H'
      self.add_vertex( h)
      self.add_edge( h, v)
      hs.add( h)
    return hs


  def explicit_hydrogens_to_real_atoms( self, v):
    hs = set()
    for i in range( v.explicit_hydrogens):
      h = self.create_vertex()
      h.symbol = 'H'
      self.add_vertex( h)
      self.add_edge( h, v)
      hs.add( h)
    v.explicit_hydrogens = 0
    return hs
  

  def add_missing_bond_orders( self, retry=False):
    """retry means to try a different approach because the last one was not successful"""
    # we must at first find orphans that will force raising of valence of their neigbors
    # NO2 is a typical example
    [v.raise_valency_to_senseful_value() for v in self.vertices if v.free_valency < 0]

    #print [sum( [xxx.free_valency for xxx in x]) for x in self._gen_free_valency_connected_components()]

##     print [map( str, n) for n in [v.neighbors for v in self.vertices if v.symbol == "N"]]
##     print [v.charge for v in self.vertices if v.symbol == "N"]

    fix = True
    while fix:
      fix = False
      for v in self.vertices:
        if sum( [n.free_valency for n in v.neighbors]) < v.free_valency:
          for n in v.neighbors:
            if n.raise_valency():
              fix = True
              break
        if fix:
          break

    #print [v.valency for v in self.vertices if v.symbol == "N"]

##     if retry:
##       for ring in self.get_smallest_independent_cycles_e():
##         if not filter( lambda x: x<=0, [v.free_valency for v in self.edge_subgraph_to_vertex_subgraph( ring)]):
##           for e in ring:
##             e.order = 4
##       self.localize_aromatic_bonds()

    # now we process it
    processed = [1]
    step = 0
    while processed:
      processed = []
      step += 1
      assert step < 100 # there is something strange if we reach this number
      for b in [bo for bo in self.edges if min( [a.free_valency for a in bo.vertices])]:
        a1, a2 = b.get_vertices()
        as1 = [a for a in a1.get_neighbors() if a.free_valency > 0]
        as2 = [a for a in a2.get_neighbors() if a.free_valency > 0]
        if len( as1) == 1 or len( as2) == 1:
          b.order += min( [a.free_valency for a in b.vertices])
          processed.append( b)

      # we process non-bridge bonds first, so that their fixation does not
      # force unsolvable possitioning of bonds in rings
      if not processed:
        for b in self.edges:
          i = min( [a.free_valency for a in b.vertices])
          if i and not self.is_edge_a_bridge( b):
            processed = [b]
            b.order += i
            break

      if not processed:
        for b in self.edges:
          i = min( [a.free_valency for a in b.vertices])
          if i:
            processed = [b]
            b.order += i
            break


  def _gen_free_valency_connected_components( self):
    done = set()
    vs = set( copy.copy( self.vertices))
    while vs:
      now = set()
      v = vs.pop()
      if v.free_valency:
        now.add( v)
        new = set( now)
        while new:
          added = False
          now_new = set()
          for v in new:
            for n in v.neighbors:
              if n.free_valency and not n in done:
                now.add( n)
                done.add( n)
                now_new.add( n)
          new = now_new
        vs -= now
        yield now
      else:
        done.add( v)




  def mark_aromatic_bonds( self):
    sssr = self.get_smallest_independent_cycles()
    if len( sssr) > 10:
      # turn off processing of all cycles - it would be too slow, just use SSSR
      rings = sssr
    else:
      rings = self.get_all_cycles()   #self.get_smallest_independent_cycles()
    solved = [1]
    # we need to repeat it to mark such things as 'badly' drawn naphtalene (no double bond in the centre) 
    while solved:
      solved = []
      for aring in rings:
        els = [self._get_atoms_possible_aromatic_electrons( a, aring) for a in aring]
        if () in els:
          continue
        for comb in common.gen_combinations_of_series( els):
          if sum( comb) % 4 == 2:
            bring = self.vertex_subgraph_to_edge_subgraph( aring)
            for b in bring:
              b.aromatic = 1
              #b.order = 4
            solved.append( aring)
            break
      for r in solved:
        rings.remove( r)




  def _get_atoms_possible_aromatic_electrons( self, at, ring):
    out = set()
    if at.charge > 0 and at.symbol not in PT.accept_cation:
      out.add( 0)
    elif at.charge < 0:
      out.add( 2)
    if at.symbol in PT.accept_cation and not at.charge > 0:
      out.add( 2)
    if at.symbol in ('B','Al'):
      out.add( 0)
    #if at.multiplicity == 2:
    #  out.add( 1)
    for b, a in at.get_neighbor_edge_pairs():
      if b.order > 1 and (a in ring or b.aromatic):
        if b.order == 4:
          # we localize aromatic bonds
          if at.get_highest_possible_free_valency() and a.get_highest_possible_free_valency():
            # we calculate free valency with all aromatic bonds as single
            out.add( 1)
        else:
          # we detect aromatic bonds
          out.add( 1)
      elif b.order == 2 and not (a in ring):
        out.add( 0)
    return tuple( out)



  def localize_aromatic_bonds( self):
    """localizes aromatic bonds (does not relocalize already localized ones),
    for those that are not aromatic but marked so
    (it is for instance possible to misuse 'cccc' in smiles to create butadiene)
    they will be properly localized but marked as non-aromatic"""
    # at first get all clusters of aromatic bonds
    for b in copy.copy( self.edges):
      if not b.aromatic:
        self.temporarily_disconnect_edge( b)
    self.temporarily_strip_bridge_edges()
    ring_clusters = map( list, [sub for sub in self.get_connected_components() if len( sub) > 1])
    self.reconnect_temporarily_disconnected_edges()
    # now proceed in localizing double bonds in each one
    for cluster in ring_clusters:
      els = [self._get_atoms_possible_aromatic_electrons( a, cluster) for a in cluster]
      if () in els:
        continue  # misuse of aromatic bonds (e.g. by smiles) or e.g. tetrahydronaphtalene
      current_bonds = self.vertex_subgraph_to_edge_subgraph( cluster)
      combinations = list( common.gen_combinations_of_series( els))
      # try most promissing combinations first (4n+2 with maximum number of 1's and maximum sum)
      combinations.sort( reverse=True, key=lambda x: sum(x)%4==2 and 1000000+1000*len([v for v in x if v==1])+sum(x) or sum(x)+1000*len([v for v in x if v==1]))
      for comb in combinations:
        if sum( comb) % 2 == 0:
          # assign the aromatic electrons to each vertex for this run
          for i,v in enumerate( cluster):
            v.properties_['arom_els'] = comb[i]
          # at first we process the bonds that are surely single (like C-S-C in thiophene)
          to_process = []
          for b in current_bonds:
            v1,v2 = b.vertices
            if v1.properties_['arom_els'] + v2.properties_['arom_els'] != 2 or \
               (b.order in (1,4) and (not v2.get_highest_possible_free_valency() or not v1.get_highest_possible_free_valency())):
              # these bonds must be single
              # either the number of electrons between the atoms is not 2 or
              # the atoms have no possible way of accomodating different than single bond
              b.order = 1
            else:
              to_process.append( b)
          to_process_atoms = self.edge_subgraph_to_vertex_subgraph( to_process)
          processed = [v for v in cluster if not v in to_process_atoms]
          # now we use the maximum_matching graph algorithm to localize the bonds
          work_graph = self.get_induced_copy_subgraph_from_vertices_and_edges( to_process_atoms, to_process, add_back_links=True)
          #return work_graph
          mate,nrex = work_graph.get_maximum_matching()
          ok = True
          if not nrex:
            # there are no exposed vertices
            for work_v,work_n in mate.items():
              if work_n != 0:
                work_e = work_v.get_edge_leading_to( work_n)
                work_e.properties_['original'].order = 2
            for e in to_process:
              if e.order == 4:
                e.order = 1
          else:
            ok = False
          if ok:
            # additional check
            for v in cluster:
              if v.properties_['arom_els'] == 1:
                for e in v.neighbor_edges:
                  if e.order == 2:
                    # this is ok
                    break
                else:
                  # this is not good
                  ok = False
                  break
          #if not (2*len( [b for b in to_process if b.order == 2])+sum([v.properties_['arom_els'] for v in cluster if v not in to_process_atoms]) == sum( comb)):
          #  ok = False
          if not ok:
            # there are some exposed vertices left
            for b in self.vertex_subgraph_to_edge_subgraph( cluster):
              b.order = 4
          else:
            break
      if [b for b in current_bonds if b.order == 4]:
        # we did not find a matching
        raise ValueError( "Localization of aromatic bonds failed")
    self.reconnect_temporarily_disconnected_edges()
    self.localize_fake_aromatic_bonds()
    



  def localize_aromatic_bonds_old( self):
    """localizes aromatic bonds (does not relocalize already localized ones),
    for those that are not aromatic but marked so
    (it is for instance possible to misuse 'cccc' in smiles to create butadiene)
    they will be properly localized but marked as non-aromatic"""
    erings = self.get_smallest_independent_cycles_e()
    # filter off rings without aromatic bonds
    erings = filter( lambda x: len( [b for b in x if b.aromatic]), erings)
    rings = map( self.edge_subgraph_to_vertex_subgraph, erings)
    # sort rings
    rings.sort( lambda x,y: len(y)%2 - len(x)%2) # odd size rings first
    last_rings = []
    while rings:
      # we have to continue with the neighbor rings of the last_ring (or the rings before)
      intersection = None
      if last_rings:
        aring = None
        while last_rings:
          last_ring = last_rings.pop( -1)
          for ring in rings:
            intersection = ring & last_ring
            if intersection:
              aring = ring
              last_rings.append( last_ring)
              break
          if aring:
            break
        if not aring:
          aring = rings.pop(0)
        else:
          rings.remove( aring)
      else:
        aring = rings.pop(0)
      last_rings.append( aring)
      # convert ring from set to list
      aring = self.sort_vertices_in_path( aring, start_from=intersection and intersection.pop() or None)
      # taken from mark_aromatic_bonds
      els = [self._get_atoms_possible_aromatic_electrons( a, aring) for a in aring]
      if () in els:
        continue  # misuse of aromatic bonds (e.g. by smiles) or e.g. tetrahydronaphtalene
      for comb in common.gen_combinations_of_series( els):
        if sum( comb) % 4 == 2:
          aring.append( aring[0])
          comb.append( comb[0])
          # at first we process the bonds that are surely single (like C-S-C in thiophene)
          already_set = None
          for i in range( len( aring) -1):
            if comb[i] + comb[i+1] != 2:
              # these bonds must be single
              b = aring[i].get_edge_leading_to( aring[i+1])
              b.order = 1
              already_set = i
          if already_set != None:
            # we reorder the comb and aring to start from the already set bonds
            j = already_set + 1
            aring = aring[j:len(aring)] + aring[1:j]
            aring.insert( 0, aring[-1])
            comb = comb[j:len(comb)] + comb[1:j]
            comb.insert( 0, comb[-1])
          i = 0
          while i+1 < len( aring):
            if comb[i] + comb[i+1] == 2:
              b = aring[i].get_edge_leading_to( aring[i+1])
              assert b != None # should be
              # to assure alternating bonds
              bs1 = [bo for bo in aring[i].get_neighbor_edges() if bo.order == 2 and bo.aromatic and bo!=b]
              bs2 = [bo for bo in aring[i+1].get_neighbor_edges() if bo.order == 2 and bo.aromatic and bo!=b]
              if len( bs1) == 0 and len( bs2) == 0:
                b.order = 2
              else:
                b.order = 1
            i += 1
          break
##           # ensure chemical meaning
##           if sum( comb) == 2*sum( [b.order-1 for b in self.vertex_subgraph_to_edge_subgraph( aring[:-1])]):
##             break
##           else:
##             print "AAAA"

    self.localize_fake_aromatic_bonds()




  def localize_fake_aromatic_bonds( self):
    to_go = [b for b in self.bonds if b.order == 4]

    processed = []
    for b in to_go:
      if not min( [a.free_valency for a in b.vertices]):
        b.order = 1
        processed.append( b)
    to_go = misc.difference( to_go, processed)

    while to_go:
      # find the right bond
      b = None
      for bo in to_go:
        bs1, bs2 = bo.get_neighbor_edges2()
        if not bs1 or len( [e for e in bs1 if e.order == 1]) > 0 and len( [e for e in bs1 if e.order == 2]) == 0 \
           or not bs2 or len( [e for e in bs2 if e.order == 1]) > 0 and len( [e for e in bs2 if e.order == 2]) == 0:
          b = bo
          break
        # new start for iteration
      if not b:
        for bo in self.edges:
          if not [e for e in bo.get_neighbor_edges() if e.order != 4]:
            b = bo
            break
      if not b:
        b = to_go.pop(0)
      # the code itself
      b.order = 2
      b.aromatic = 0
      for bo in b.get_neighbor_edges():
        if bo.order == 4:
          bo.order = 1
          bo.aromatic = 0
      # next turn
      to_go = [b for b in self.bonds if b.order == 4]


  def remove_unimportant_hydrogens( self):
    """removes all H atoms that don't have any special properties
    or are not part of a hydrogen only molecule."""
    only_hydrogens = True
    for v in self.vertices:
      if v.symbol != "H":
        only_hydrogens = False
        break
    if not only_hydrogens:
      for v in copy.copy( self.vertices):
        if v.symbol == 'H' and v.degree <= 1 and not v.isotope:
          self.remove_vertex( v)



  def remove_all_hydrogens( self):
    """removes all H atoms"""
    for v in copy.copy( self.vertices):
      if v.symbol == 'H' and v.degree <= 1:
        self.remove_vertex( v)



  def _get_atom_distance_matrix( self, a):
    self.mark_vertices_with_distance_from( a)
    big_out = []
    i = 0
    while 1:
      out = [a.symbol_number for a in self.vertices if a.properties_['d'] == i]
      if i > 0:
        out2 = []
        for a in self.vertices:
          if a.properties_['d'] == i-1:
            out2 += a.get_hydrogen_count() * [1]
        out += out2
      if out:
        out.sort()
        # out.reverse()
        big_out.append( tuple( out))
      else:
        return big_out
      i += 1




  def get_symmetry_unique_atoms( self):
    out = []
    vs = copy.copy( self.vertices)
    for v in vs:
      out.append( [i for i in self._get_atom_distance_matrix( v)])
    while out:
      k = vs.pop( 0)
      v = out.pop( 0)
      o = [k]
      proc_vs = []
      proc_out = []
      while out:
        k2 = vs.pop( 0)
        v2 = out.pop( 0)
        if v2 == v:
          o.append( k2)
        else:
          proc_vs.append( k2)
          proc_out.append( v2)
      yield o
      out = proc_out
      vs = proc_vs




  def number_atoms_uniquely( self):
    out = {}
    for v in self.vertices:
      out[v] = [i for i in self._get_atom_distance_matrix( v)]
      v.properties_['distance_matrix'] = out[v]
    ms = out.values()
    #self.remove_all_hydrogens()
    ms.sort()
    i = 1
    ret = []
    for m in ms:
      ret.append( out.keys()[ out.values().index( m)])
      #if ret.symbol != 'H':
      #  print "%3d: " % i, ret, ret.free_valency
      i += 1
    return ret




  def _read_file( self, name="/home/beda/oasa/oasa/mol.graph"):
    self.vertices = []
    self.edges = set()
    f = file( name, 'r')
    vs = f.readline()
    for i in vs.split(' '):
      if i != '\n':
        v = self.create_vertex()
        v.symbol = i
        self.add_vertex( v)
    for l in f.readlines():
      o, a, b = l.split(' ')
      e = self.create_edge()
      e.order = int( o)
      self.add_edge( self.vertices[int(a)], self.vertices[int(b)], e=e)
    f.close()




  def find_longest_mostly_carbon_chain( self):
    if len( self.vertices) < 2:
      return copy.copy( self.vertices)
    ends = [v for v in self.vertices if v.get_degree() == 1]
    paths = []
    for e1 in ends:
      for e2 in ends:
        if e1 != e2:
          path = self.find_path_between( e1, e2)
          paths.append( path)  # could be optimized by not storing all the paths (once we know how to weight the chains)
    # for now we don't take care of chain composition
    l = 0
    p = None
    for path in paths:
      if len( path) > l:
        l = len( path)
        p = path
    return p


  def remove_zero_order_bonds( self):
    """zero order bonds are created for example when SMILES is read"""
    zero_bonds = [e for e in self.edges if e.order == 0]
    for b in zero_bonds:
      self.disconnect_edge( b)



  # --- the fragment matching routines ---

  def select_matching_substructures( self, other, implicit_freesites=False, auto_cleanup=True):
    """select fragments that match the complete molecule 'other' and yield them
    as lists of atoms in the order of other.vertices; however when other has
    explicit hydrogens that match implicit hydrogens on self the length of the
    returned fragment might be shorter of the matched implicit hydrogens;
    
    IF YOU DON'T LET THE GENERATOR RUN TO THE END OR SET THE auto_cleanup TO FALSE, YOU HAVE TO
    RUN clean_after_search MANUALLY NOT TO CLUTTER THE STRUCTURE"""

    # at first decide if we need to add implicit hydrogens to self
    add_implicit = False
    for v in other.vertices:
      if (isinstance( v, atom) and v.symbol == 'H') or \
         (isinstance( v, atom) and v.explicit_hydrogens > 0) or \
         (isinstance( v, query_atom) and ('H' in v.symbols or 'R' in v.symbols)):
        add_implicit = True
        break
    if add_implicit:
      added_hs = self.add_missing_hydrogens()
      for h in added_hs:
        h.properties_['implicit_hydrogen'] = True
    # here we add explicit hydrogens to other, we call them implicit_hydrogen to make it
    # work the same way as implicit_hydrogen
    for v in other.vertices:
      for i in range( v.explicit_hydrogens):
        h = other.create_vertex()
        other.add_vertex( h)
        h.symbol = 'H'
        e = other.create_edge()
        e.order = 1
        other.add_edge( v, h, e=e)
        h.properties_['implicit_hydrogen'] = True
      v.explicit_hydrogens = 0  # make it 0 and after search put it back
    # for implicit_freesites we add the free_sites here
    if implicit_freesites:
      for v in other.vertices:
        v.properties_['old_free_sites'] = v.free_sites
        v.free_sites = v.free_valency

    # then we create the dicts for storing threads in each of the atoms
    i = 0 
    for a in other.vertices:
      a.properties_['subsearch'] = {}
    for e in other.edges | self.edges:
      e.properties_['subsearch'] = {}
    # here we select the vertex to start from
    vs = [v for v in other.vertices if isinstance( v, atom)]
    sym = common.least_common_item( [v.symbol for v in vs])
    v = [v for v in vs if v.symbol == sym][0]
    for a in self.vertices:
      a.properties_['subsearch'] = {}
      if v.matches( a):
        i += 1
        a.properties_['subsearch'][i] = v
        v.properties_['subsearch'][i] = a

    # now we can proceed with the search
    yielded = set()
    for thread in self._mark_matching_threads( v, other):
      vs = [v.properties_['subsearch'][thread] for v in other.vertices]
      # for symetrical fragments we have to get rid of copies (O1=N=O2 and O2=N=O1)
      vsset = frozenset( vs)
      if vsset not in yielded:
        if self._freesites_match( other, thread):
          yield [v for v in vs if not 'implicit_hydrogen' in v.properties_.keys()]
      yielded.add( vsset)

    if auto_cleanup:
      self.clean_after_search( other)



  def clean_after_search( self, other):
    # for implicit_freesites we restore original free_sites here
    for v in other.vertices:
      if "old_free_sites" in v.properties_:
        v.free_sites = v.properties_['old_free_sites']
        del v.properties_['old_free_sites']
    # finally we remove the added implicit hydrogens
    hs = [v for v in self.vertices if 'implicit_hydrogen' in v.properties_.keys()]
    for v in hs:
      del v.properties_['implicit_hydrogen']
      self.remove_vertex( v)
    # we also remove implicit_hydrogen (that represent atom.explicit_hydrogens) on other
    # and also set the explicit_hydrogens count back to the original atoms
    for v in other.vertices:
      if 'implicit_hydrogen' not in v.properties_:
        v.explicit_hydrogens = len( [h for h in v.neighbors if 'implicit_hydrogen' in h.properties_])
    hs = [v for v in other.vertices if 'implicit_hydrogen' in v.properties_.keys()]
    for v in hs:
      del v.properties_['implicit_hydrogen']
      other.remove_vertex( v)
    for v in self.vertices + other.vertices:
      del v.properties_['subsearch']
    for e in self.edges | other.edges:
      del e.properties_['subsearch']
    
      
  def _mark_matching_threads( self, v, other):
    """v is other vertex, other is the other molecule"""
    thread = 0
    threads = v.properties_['subsearch'].keys()
    while threads:
      thread = min( threads)
      threads.remove( thread)

      mirror = v.properties_['subsearch'][thread]
      for e, n in v.get_neighbor_edge_pairs():
        if thread not in n.properties_['subsearch']:
          candidates = set()
          for me, mn in mirror.get_neighbor_edge_pairs():
            if thread not in mn.properties_['subsearch'] and mn.matches( n) and me.matches( e) and thread not in e.properties_['subsearch']:
              candidates.add( (mn, me, e))

          if candidates:
            if len( candidates) > 1:
              new_threads = self._spawn_thread( other, thread, len( candidates)-1)
              ths = [thread] + new_threads
            else:
              ths = [thread]
            for c, me, e in candidates:
              th = ths.pop()
              n.properties_['subsearch'][th] = c
              c.properties_['subsearch'][th] = n
              e.properties_['subsearch'][th] = me
              me.properties_['subsearch'][th] = e
            [x for x in self._mark_matching_threads( n, other)] # just make the generator run
            if thread not in v.properties_['subsearch']:
              # the thread already died
              break
            else:
              pass
          else:
            self._delete_thread( other, thread)
            break
        # for proper handling of rings we have to check also the ones that are in this thread already
        elif thread not in e.properties_['subsearch']:
          me = self.get_edge_between( mirror, n.properties_['subsearch'][thread]) 
          if me and e.matches( me):
            e.properties_['subsearch'][thread] = me
            me.properties_['subsearch'][thread] = e
          else:
            self._delete_thread( other, thread)
            break
        else:
          pass
         
      threads = [i for i in v.properties_['subsearch'].keys() if i >= thread]
      if thread in threads:
        threads.remove( thread)
        yield thread
          


  def _spawn_thread( self, other, thread, number):
    my_vs = [v for v in self.vertices if thread in v.properties_['subsearch']]
    my_es = [e for e in self.edges if thread in e.properties_['subsearch']]
    other_vs = [v for v in other.vertices if thread in v.properties_['subsearch']]
    other_es = [e for e in other.edges if thread in e.properties_['subsearch']]
    max_thread = max( [max( v.properties_['subsearch'].keys()) for v in other.vertices if v.properties_['subsearch']])
    for i in range( max_thread +1, max_thread +number +1, 1):
      for v in my_vs + my_es + other_vs + other_es:
        v.properties_['subsearch'][ i] = v.properties_['subsearch'][thread]
    return range( max_thread +1, max_thread +number +1, 1)



  def _delete_thread( self, other, thread):
    for v in self.vertices + other.vertices:
      try:
        del v.properties_['subsearch'][thread]
      except KeyError:
        pass


  def _freesites_match( self, other, thread):
    for v in other.vertices:
      mirror = v.properties_['subsearch'][thread]
      unmatched_ns = [n for n in mirror.neighbors if thread not in n.properties_['subsearch']] # and not n.symbol == 'H']
      # if there is more unmatched neighbors then free-site it does not match (Hs don't count)
      if not len( unmatched_ns)+mirror.explicit_hydrogens <= v.free_sites:
        return False
    return True


  def contains_substructure( self, other, implicit_freesites=True):
    found = False
    for i in self.select_matching_substructures( other, implicit_freesites=implicit_freesites, auto_cleanup=False):
      found = True
      break
    self.clean_after_search( other)
    return found


  # // --- end of the fragment matching routines ---

  def detect_stereochemistry_from_coords( self, omit_rings=True):
    from . import stereochemistry, geometry
    def add_neighbor_double_bonds( bond, path):
      for _e in bond.get_neighbor_edges():
        if _e.order == 2 and _e not in path:
          path.append( _e)
          add_neighbor_double_bonds( _e, path)
    # double bonds
    # detect clusters of double bonds
    double_paths = []
    processed = set()
    for e in self.edges:
      if e.order == 2 and e not in processed:
        if omit_rings and not self.is_edge_a_bridge( e):
          continue
        path = [e]
        add_neighbor_double_bonds( e, path)
        if len( path) % 2:
          double_paths.append( path)
          processed |= set( path)
    # detect config on these paths
    for path in double_paths:
        vertices = []
        for bond in path:
          vertices.extend( bond.vertices)
        ends = [v for v in vertices if vertices.count(v) == 1]
        if len( ends) != 2: # two ends is the only thing we are prepared to handle
          continue
        end1, end2 = ends
        # set stereochemistry for all neighbors of both ends
        for e1,n1 in end1.get_neighbor_edge_pairs():
          plane1 = geometry.plane_normal_from_3_points( (n1.x,n1.y,n1.z),(end1.x,end1.y,end1.z),(end2.x,end2.y,end2.z))
          if plane1 == None:
            continue # some coords were missing
          if not e1 in path:
            for e2,n2 in end2.get_neighbor_edge_pairs():
              if not e2 in path:
                plane2 = geometry.plane_normal_from_3_points( (end1.x,end1.y,end1.z),(end2.x,end2.y,end2.z),(n2.x,n2.y,n2.z))
                #cos_angle = geometry.same_or_oposite_side( plane1, plane2)
                cos_angle = geometry.angle_between_planes( plane1, plane2)
                if cos_angle < 0:
                  value = stereochemistry.cis_trans_stereochemistry.OPPOSITE_SIDE
                else:
                  value = stereochemistry.cis_trans_stereochemistry.SAME_SIDE
                if len( path) == 1:
                  center = path[0]
                else:
                  center = None
                refs = [n1,end1,end2,n2]
                st = stereochemistry.cis_trans_stereochemistry( center=center, value=value, references=refs)
                to_remove = None
                to_add = None
                for st1 in self.stereochemistry:
                  if set( st1.references) == set( st.references):
                    if st.value == st1.value:
                      break
                    else:
                      to_remove = st1
                      break
                else:
                  self.add_stereochemistry( st)
                if to_remove:
                  self.remove_stereochemistry( to_remove)


  def mark_morgan( self):
    old_morgs = 0
    for v in self.vertices:
      v.properties_['new_morgan'] = v.degree #sum( [n.degree for n in v.neighbors])
    new_morgs = len( set( [v.properties_['new_morgan'] for v in self.vertices]))

    while new_morgs > old_morgs:
      for v in self.vertices:
        v.properties_['morgan'] = v.properties_['new_morgan'] 

      for v in self.vertices:
        v.properties_['new_morgan'] = v.properties_['morgan'] + sum( [n.properties_['morgan'] for n in v.neighbors])
      old_morgs = new_morgs
      new_morgs = len( set( [v.properties_['new_morgan'] for v in self.vertices]))
      
    print("morgan", old_morgs, [v.properties_['morgan'] for v in self.vertices])


  ## some geometry related things

  def normalize_bond_length( self, bond_length=30):
    """make the average bond-length be bond_length by scaling the structure up"""
    if not self.edges or len( self.vertices) < 2:
      return False

    minx = None
    maxx = None
    miny = None
    maxy = None
    # atoms
    for v in self.vertices:
      if not maxx or v.x > maxx:
        maxx = v.x
      if not minx or v.x < minx:
        minx = v.x
      if not miny or v.y < miny:
        miny = v.y
      if not maxy or v.y > maxy:
        maxy = v.y

    scale = bond_length / self.get_mean_bond_length()
    movex = (maxx+minx)/2
    movey = (maxy+miny)/2
    trans = transform3d.transform3d()
    trans.set_move( -movex, -movey, 0)
    trans.set_scaling( scale)
    trans.set_move( movex, movey, 0)
    for v in self.atoms:
      v.x, v.y, v.z = trans.transform_xyz( v.x, v.y, v.z)
    return True


  def get_mean_bond_length( self):
    """returns the mean bond length of bonds in the molecule"""
    if len( self.edges) == 0:
      return None
    bond_length = 0.0
    for e in self.edges:
      v1, v2 = e.vertices
      bond_length += math.sqrt( (v1.x-v2.x)**2 + (v1.y-v2.y)**2)
    bl = bond_length / len( self.edges)
    return bl

  def get_structure_hash( self):
    vs = self.number_atoms_uniquely()
    ret = []
    for v in vs:
      ret.append( str(v.properties_['distance_matrix']))
    res = "*".join( ret)
    import sha
    ss = sha.new()
    ss.update( res)
    return ss.hexdigest()

  def create_CIP_digraph( self, center):
    """creates a digraph according to rules described in CIP paper."""
    assert center in self.vertices
    from graph.digraph import digraph
    dg = digraph()
    dg.add_vertex( center.copy())
    return dg


def the_right_sorting_function( t1, t2):
  for i,l in enumerate( t1):
    k = t2[i]
    for j in range( len( l)):
      if j > len( k)-1:
        return 1
      if l[j] > k[j]:
        return 1
      elif l[j] < k[j]:
        return -1
    if len( l) < len( k):
      return -1
  return 0




def equals( mol1, mol2, level=0):
  """don't forget to put all hydrogens and bond orders to the bonds
     level 1 - number of atoms and bonds,
     level 2 - the number of atoms with same symbols is the same,
     level 3 - the whole connecivity (no stereo)
     """
  # level 1
  if not level or level >= 1:
    if not len( mol1.vertices) == len( mol2.vertices):
      return False
    if not len( mol1.edges) == len( mol2.edges):
      return False
  # level 2
  if not level or level >= 2:
    symbols1 = [v.symbol for v in mol1.vertices]
    symbols1.sort()
    symbols2 = [v.symbol for v in mol2.vertices]
    symbols2.sort()
    if symbols2 != symbols1:
      return False
  # level 3
  if not level or level >= 3:
    vs1 = mol1.number_atoms_uniquely()
    vs2 = mol2.number_atoms_uniquely()
    for i, v1 in enumerate( vs1):
      v2 = vs2[ i]
      if v1.properties_['distance_matrix'] != v2.properties_['distance_matrix']:
        return False
  return True

#import psyco
#psyco.profile()

##################################################
# DEMO

if __name__ == '__main__':

  def main():
    g = molecule()
    g._read_file()

    for b in g.edges:
      print(g.is_edge_a_bridge( b), end="")

##     g.add_missing_hydrogens()

##     for v in g.vertices:
##   #    print v.symbol, v.is_chiral()
##       if v.symbol:
##         ## for a in v.get_neighbors_CIP_sorted():
##         ##   print a, [na.symbol for na in a.get_neighbors()]
##         if v.is_chiral():
##           print v.symbol, [na.symbol for na in v.get_neighbors_CIP_sorted()]



  import time
  t = time.time()
  main()
  print(round( 1000*(time.time()-t), 2), 'ms')

# DEMO END
##################################################



##################################################
# TODO

# azulene is wrongly localized
# 0,2 must alternate in aromatic compound
# radical has 1 el in _get_atoms_possible_aromatic_electrons
# test the equals function

##################################################


