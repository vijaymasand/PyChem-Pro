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

from .plugin import plugin
from .molecule import molecule, equals
from . import periodic_table as PT
from . import oasa_exceptions
from . import reaction
from . import stereochemistry


from .config import Config

import re

import operator
import time
from functools import reduce

class smiles( plugin):

  name = "smiles"
  read = 1
  write = 1

  smiles_to_oasa_bond_recode = {'-': 1, '=': 2, '#': 3, ':': 4, ".": 0, "\\": 1, "/": 1}
  oasa_to_smiles_bond_recode = {1: '', 2: '=', 3: '#', 4:''}

  organic_subset = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']

  def __init__( self, structure=None):
    self.structure = structure

  def recode_oasa_to_smiles_bond( self, b):
    if b.aromatic:
      return ''
    elif b in self._stereo_bonds_to_others:
      others = [(e,st) for e,st in self._stereo_bonds_to_others[b] if e in self._stereo_bonds_to_code]
      if not others:
        code = "\\"
      else:
        # other bonds enforce encoding of this one, we select the first one,
        # bacause there is nothing we can do if there are clashing constrains anyway
        other, st = others[0]
        end1,inside1,inside2,end2 = st.references
        if set( other.vertices) == set( [end1,inside1]):
          v1 = inside1
          v2 = end1
        else:
          v1 = inside2
          v2 = end2
        last_order = self._processed_atoms.index( v2) - self._processed_atoms.index( v1)
        last_code = self._stereo_bonds_to_code[ other] == "\\" and 1 or -1
        relation = st.value == st.OPPOSITE_SIDE and -1 or 1
        if relation*last_code*last_order < 0:
          code = "/"
        else:
          code = "\\"
      self._stereo_bonds_to_code[ b] = code
      return code
    else:
      if b.order == 1:
        a1, a2 = b.vertices
        if 'aromatic' in a1.properties_ and 'aromatic' in a2.properties_:
          # non-aromatic bond connecting two aromatic rings, we need to return -
          return '-'
      return self.oasa_to_smiles_bond_recode[ b.order]

  def set_structure( self, structure):
    self.structure = structure

  def get_structure( self):
    return self.structure

  def read_smiles( self, text, explicit_hydrogens_to_real_atoms=False):
    self.explicit_hydrogens_to_real_atoms = explicit_hydrogens_to_real_atoms
    mol = Config.create_molecule()
    text = "".join( text.split())
    is_text = re.compile("^[A-Z][a-z]?$")
    # internally revert \/ bonds before numbers, this makes further processing much easier
    text = re.sub( r"([\\/])([0-9])", lambda m: (m.group(1)=="/" and "\\" or "/")+m.group(2), text)
    # // end
    chunks = re.split( r"(\[.*?\]|[A-Z][a-z]?|%[0-9]{1,2}|[^A-Z]|[a-z])", text)
    chunks = self._check_the_chunks( chunks)
    last_atom = None
    last_bond = None
    numbers = {}
    bracket_openings = []
    for c in chunks:
      # atom
      if is_text.match( c) or c.islower() or c[0] == "[":
        a = mol.create_vertex()
        if c[0] == "[":
          # atom spec in square brackets
          self._parse_atom_spec( c, a)
        else:
          # just atom symbol
          if c.islower():
            symbol = c.capitalize()
            a.properties_['aromatic'] = 1
          else:
            symbol = c
          a.symbol = symbol

        mol.add_vertex( a)
        if last_bond: # and not (not 'aromatic' in a.properties_ and last_bond.aromatic):
          # make last bond aromatic if it was stereo and atoms are aromatic
          if 'stereo' in last_bond.properties_ and \
             'aromatic' in last_atom.properties_ and \
             'aromatic' in a.properties_ and \
             not last_bond.aromatic:
            last_bond.aromatic = True
          mol.add_edge( last_atom, a, e=last_bond)
          last_bond = None
        elif last_atom:
          b = mol.add_edge( last_atom, a)
          if 'aromatic' in a.properties_:
            # aromatic bond
            b.order = 4
            b.type = 'n'
        last_atom = a
        last_bond = None
      # bond
      elif c in r'-=#:.\/':
        order = self.smiles_to_oasa_bond_recode[ c]
        last_bond = mol.create_edge()
        last_bond.order = order
        last_bond.type = 'n'
        if c in r'\/':
          last_bond.properties_['stereo'] = c
      # ring closure
      elif c.isdigit():
        if c in numbers:
          if last_bond:
            b = last_bond
          else:
            b = mol.create_edge()
            if "aromatic" in numbers[c].properties_:
              b.order = 4
          mol.add_edge( last_atom, numbers[c], e=b)
          last_bond = None
          del numbers[ c]
        else:
          numbers[c] = last_atom
          last_bond = None
      elif c == '(':
        bracket_openings.append( last_atom)
      elif c == ')':
        last_atom = bracket_openings.pop(-1)

    ## FINISH
    # deal with explicit valency, etc.
    for a in mol.vertices:
      if not "explicit_valency" in a.properties_:
        a.raise_valency_to_senseful_value()
      else:
        # detect radicals (but not biradicals - problem of triplet vs. singlet)
        del a.properties_['explicit_valency']
        if a.valency - a.occupied_valency == 1:
          a.multiplicity += 1
        else:
          a.valency = a.occupied_valency
      try:
        del a.properties_['aromatic']
      except:
        pass

    # stereochemistry
    self._process_stereochemistry( mol)

    if len(mol.vertices) == 0:
      mol = None
    self.structure = mol 


  def _parse_atom_spec( self, c, a):
    """c is the text spec,
    a is an empty prepared vertex (atom) instance"""
    bracketed_atom = re.compile(r"^\[(\d*)([A-z][a-z]?)(.*?)\]")
    m = bracketed_atom.match( c)
    if m:
      isotope, symbol, rest = m.groups()
    else:
      raise ValueError( "unparsable square bracket content '%s'" % c)
    if symbol.islower():
      symbol = symbol.capitalize()
      a.properties_['aromatic'] = 1
    a.symbol = symbol
    if isotope:
      a.isotope = int( isotope)
    # hydrogens
    _hydrogens = re.search( r"H(\d*)", rest)
    h_count = 0
    if _hydrogens:
      if _hydrogens.group(1):
        h_count = int( _hydrogens.group(1))
      else:
        h_count = 1
    a.explicit_hydrogens = h_count
    # charge
    charge = 0
    # one possible spec of charge
    _charge = re.search( "[-+]{2,10}", rest)
    if _charge:
      charge = len( _charge.group(0))
      if _charge.group(0)[0] == "-":
        charge *= -1
    # second one, only if the first one failed
    else:
      _charge = re.search( r"([-+])(\d?)", rest)
      if _charge:
        if _charge.group(2):
          charge = int( _charge.group(2))
        else:
          charge = 1
        if _charge.group(1) == "-":
          charge *= -1
    a.charge = charge
    # stereo
    _stereo = re.search( "@+", rest)
    if _stereo:
      stereo = _stereo.group(0)
      a.properties_['stereo'] = stereo
    # using [] means valency is explicit
    a.properties_['explicit_valency'] = True

  def _check_the_chunks( self, chunks):
    is_text = re.compile("^[A-Z][a-z]?$")
    is_long_num = re.compile( "^%[0-9]{1,2}$")
    ret = []
    for c in chunks:
      if not c:
        continue
      if is_text.match( c):
        if (not c in PT.periodic_table and len(c)==2) or c == "Sc": # Sc is S-c not scandium
          a,b = c
          ret.append( a)
          ret.append( b)
        else:
          ret.append( c)
      elif is_long_num.match( c):
        ret.append( str( int( c[1:])))
      else:
        ret.append( c)
    return ret



  def _process_stereochemistry( self, mol):
    ## process stereochemistry
    ## double bonds
    def get_stereobond_direction( end_atom, inside_atom, bond, init):
      position = mol.vertices.index( end_atom) - mol.vertices.index( inside_atom)
      char = bond.properties_['stereo'] == "\\" and 1 or -1
      direction = (position * char * init) < 0 and "up" or "down"
      return direction
    def get_end_and_inside_vertex_from_edge_path( edge, path):
      a1,a2 = edge.vertices
      if len( [e for e in a1.neighbor_edges if e in path]) == 1:
        return a1, a2
      return a2, a1
    
    stereo_edges = [e for e in mol.edges if "stereo" in e.properties_]
    paths = []
    for i,e1 in enumerate( stereo_edges):
      for e2 in stereo_edges[i+1:]:
        path = mol.get_path_between_edges( e1, e2)
        path2 = path[1:-1]
        if len( path2)%2 and not [_e for _e in path2 if _e.order != 2]:
          # only odd number of double bonds, double bonds only
          for _e in path[1:-1]:
            if not mol.is_edge_a_bridge_fast_and_dangerous( _e):
              break
          else:
            # only stereo related to non-cyclic bonds
            paths.append( path)
    
    for path in paths:
      bond1 = path[0]
      end_atom1,inside_atom1 = get_end_and_inside_vertex_from_edge_path( bond1, path)
      bond2 = path[-1]
      end_atom2,inside_atom2 = get_end_and_inside_vertex_from_edge_path( bond2, path)
      d1 = get_stereobond_direction( end_atom1, inside_atom1, bond1, -1)
      d2 = get_stereobond_direction( end_atom2, inside_atom2, bond2, -1)
      if d1 == d2:
        value = stereochemistry.cis_trans_stereochemistry.SAME_SIDE
      else:
        value = stereochemistry.cis_trans_stereochemistry.OPPOSITE_SIDE
      if len( path) == 3:
        center = path[1]
      else:
        center = None
      refs = [end_atom1,inside_atom1,inside_atom2,end_atom2]
      st = stereochemistry.cis_trans_stereochemistry( center=center, value=value, references=refs)
      mol.add_stereochemistry( st)

    # tetrahedral stereochemistry
    for v in mol.vertices:
      refs = None
      if 'stereo' in v.properties_:
        idx = [mol.vertices.index( n) for n in v.neighbors]
        idx.sort()
        if len( idx) < 3:
          pass # no stereochemistry with less then 3 neighbors
        elif len( idx) == 3:
          if v.explicit_hydrogens == 0:
            pass # no stereochemistry without adding hydrogen here
          else:
            if self.explicit_hydrogens_to_real_atoms:
              hs = mol.explicit_hydrogens_to_real_atoms( v)
              h = hs.pop()
            else:
              h = stereochemistry.explicit_hydrogen()
            v_idx = mol.vertices.index( v)
            idx1 = [i for i in idx if i < v_idx]
            idx2 = [i for i in idx if i > v_idx]
            refs = [mol.vertices[i] for i in idx1] + [h] + [mol.vertices[i] for i in idx2]
        elif len( idx) == 4:
          refs = [mol.vertices[i] for i in idx]
        else:
          pass # unhandled stereochemistry
      if refs:
        if v.properties_["stereo"] == "@":
          direction = stereochemistry.tetrahedral_stereochemistry.ANTICLOCKWISE
        elif v.properties_['stereo'] == "@@":
          direction = stereochemistry.tetrahedral_stereochemistry.CLOCKWISE
        else:
          continue # no meaning
        st = stereochemistry.tetrahedral_stereochemistry( center=v, value=direction, references=refs)
        mol.add_stereochemistry( st)

    # delete the data after processing
    for e in mol.edges:
      if 'stereo' in e.properties_:
        del e.properties_['stereo']
    for v in mol.vertices:
      if 'stereo' in v.properties_:
        del v.properties_['stereo']


  def get_smiles( self, mol):
    if not mol.is_connected():
      raise oasa_exceptions.oasa_not_implemented_error( "SMILES", "Cannot encode disconnected compounds, such as salts etc. HINT - use molecule.get_disconnected_subgraphs() to divide the molecule to individual parts.")
    #mol = molec.copy()
    self.molecule = mol
    self.ring_joins = []
    self._processed_atoms = []
    self.branches = {}
    self._stereo_bonds_to_code = {} # for bond it will contain character it uses
    self._stereo_bonds_to_others = {} # for bond it will contain the other bonds
    self._stereo_centers = {}
    # at first we mark all the atoms with aromatic bonds
    # it is much simple to do it now when all the edges are present
    # we can make use of the properties attribute of the vertex
    for e in mol.edges:
      if e.aromatic:
        for v in e.vertices:
          v.properties_[ 'aromatic'] = 1
    # stereochemistry information preparation
    for st in mol.stereochemistry:
      if isinstance( st, stereochemistry.cis_trans_stereochemistry):
        end1, inside1, inside2, end2 = st.references
        e1 = end1.get_edge_leading_to( inside1)
        e2 = end2.get_edge_leading_to( inside2)
        self._stereo_bonds_to_others[ e1] = self._stereo_bonds_to_others.get( e1, []) + [(e2, st)]
        self._stereo_bonds_to_others[ e2] = self._stereo_bonds_to_others.get( e2, []) + [(e1, st)]
      elif isinstance( st, stereochemistry.tetrahedral_stereochemistry):
        self._stereo_centers[st.center] = st
      else:
        pass # we cannot handle this
        
    ret = ''.join( [i for i in self._get_smiles( mol)])
    mol.reconnect_temporarily_disconnected_edges()
    # this is needed because the way temporarily_disconnected edges are handled is not compatible with the way smiles
    # generation works - it splits the molecule while reusing the same atoms and bonds and thus disconnected bonds accounting fails
    for e in mol.edges:
      e.disconnected = False
    # here tetrahedral stereochemistry is added
    for v, st in self._stereo_centers.items():
      processed_neighbors = []
      for n in self._processed_atoms:
        if n in v.neighbors:
          processed_neighbors.append( n)
        elif v.explicit_hydrogens and n is v:
          processed_neighbors.append( stereochemistry.explicit_hydrogen())
      count = match_atom_lists( st.references, processed_neighbors)
      clockwise = st.value == st.CLOCKWISE
      if count % 2 == 1:
        clockwise = not clockwise
      ch_symbol = clockwise and "@@" or "@"
      ret = ret.replace( "{{stereo%d}}" % mol.vertices.index(v), ch_symbol)
    return ret



  def _get_smiles( self, mol, start_from=None):
    # single atoms
    if len( mol.vertices) == 1:
      v = mol.vertices[0]
      yield self._create_atom_smiles( v)
      for e in self.ring_joins:
        if v in e.get_vertices():
          yield self.recode_oasa_to_smiles_bond( e)
          yield self._create_ring_join_smiles( self.ring_joins.index( e))
      return
    while not (is_line( mol) and (not start_from or start_from.get_degree() <= 1)):
      if is_pure_ring( mol):
        if start_from:
          self.ring_joins.append( mol.temporarily_disconnect_edge( start_from.neighbor_edges[0]))
        else:
          self.ring_joins.append( mol.temporarily_disconnect_edge( list( mol.edges)[0]))
      else:
        e, mol, branch_vertex, branch = self.disconnect_something( mol, start_from=start_from)
        if branch_vertex:
          if branch_vertex in self.branches:
            self.branches[ branch_vertex].append((e, branch))
          else:
            self.branches[ branch_vertex] = [(e, branch)]
        else:
          self.ring_joins.append( e)
    try:
      start, end = filter( lambda x: x.get_degree() == 1, mol.vertices)
    except:
      #print filter( lambda x: x.get_degree() == 1, mol.vertices)
      raise Exception("shit")
    if start_from == end:
      start, end = end, start
    v = start
    last = None
    was_end = False
    while True:
      yield self._create_atom_smiles( v)
      # the atom
      for e in self.ring_joins:
        if v in e.get_vertices():
          _b = self.recode_oasa_to_smiles_bond( e)
          if _b not in "/\\":
            yield _b
          yield self._create_ring_join_smiles( self.ring_joins.index( e))
      # branches
      if v in self.branches:
        for edg, branch in self.branches[ v]:
          yield '('
          yield self.recode_oasa_to_smiles_bond( edg)
          v1, v2 = edg.vertices
          vv = (v1 != v) and v1 or v2
          for i in self._get_smiles( branch, start_from=vv):
            yield i
          yield ')'
      # bond leading to the neighbor
      for e, neighbor in v.get_neighbor_edge_pairs():
        if neighbor != last:
          yield self.recode_oasa_to_smiles_bond( e)
          last = v
          v = neighbor
          break
      if was_end:
        break
      if v == end:
        was_end = True


  def _create_atom_smiles( self, v):
    self._processed_atoms.append( v)
    if 'aromatic' in v.properties_.keys():
      symbol = v.symbol.lower()
    else:
      symbol = v.symbol

    stereo = self._stereo_centers.get( v, None)

    if (v.symbol not in self.organic_subset) or (v.isotope) or (v.charge != 0) or (v.valency != PT.periodic_table[ v.symbol]['valency'][0]) or (stereo) or (v.multiplicity != 1):
      # we must use square bracket
      isotope = v.isotope and str( v.isotope) or ""
      # charge
      if v.charge:
        sym = v.charge < 0 and "-" or "+"
        charge = sym + (abs( v.charge) > 1 and str( abs( v.charge)) or "")
      else:
        charge = ""
      # explicit hydrogens
      num_h = v.valency - v.occupied_valency + v.explicit_hydrogens
      h_spec = (num_h and "H" or "") + (num_h > 1 and str( num_h) or "")
      # stereo
      if stereo:
        stereo = "{{stereo%d}}" % self.molecule.vertices.index( v)
      else:
        stereo = ""
      return "[%s%s%s%s%s]" % (isotope, symbol, stereo, h_spec, charge)
    else:
      # no need to use square brackets
      return symbol


  def disconnect_something( self, mol, start_from=None):
    """returns (broken edge, resulting mol, atom where mol was disconnected, disconnected branch)"""
    # we cannot do much about this part
    if start_from and start_from.get_degree() != 1:
      for e,n in start_from.get_neighbor_edge_pairs():
        if n.get_degree() > 2:
          mol.temporarily_disconnect_edge( e)
          return e, mol, None, None
      mol.temporarily_disconnect_edge( e)
      return e, mol, None, None
    # at first try to find edges for which degree of neighbors is bigger
    # than [2,2] and at best they are not bridges
    # when no non-bridges are present use the other ones
    #
    # the edges with crowded atoms
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if d1 > 2 and d2 > 2 and not mol.is_edge_a_bridge_fast_and_dangerous( e):
        mol.temporarily_disconnect_edge( e)
        return e, mol, None, None
    # the other valuable non-bridge edges
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if (d1 > 2 or d2 > 2) and not mol.is_edge_a_bridge_fast_and_dangerous( e):
        mol.temporarily_disconnect_edge( e)
        return e, mol, None, None
    # there are no non-bridges
    # we want to split off the smallest possible chunks
    min_size = None
    the_right_edge = None
    the_right_mol = None
    the_right_branch = None
    the_right_branch_atom = None
    ring_joints_in_branch = 1000
    ring_join_vertices = set( reduce( operator.add, [e.vertices for e in self.ring_joins], []))
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if d1 > 2 or d2 > 2:
        ps = mol.get_pieces_after_edge_removal( e)
        if len( ps) == 1:
          print("impossible")
          continue
        lenghts = map( len, ps)
        ms = min( lenghts)
        p1, p2 = ps
        the_mol = (len( p1) < len( p2)) and p2 or p1
        the_branch = (p1 == the_mol) and p2 or p1
        ring_joints = len( [i for i in the_branch if i in ring_join_vertices])
        if not min_size or ms < min_size or ring_joints_in_branch > ring_joints:
          min_size = ms
          the_right_edge = e
          the_right_mol = the_mol
          the_right_branch = the_branch
          ring_joints_in_branch = ring_joints
    if the_right_edge:
      # what is possible to make here instead in the loop is made here
      # it saves time
      v1, v2 = the_right_edge.vertices
      mol.temporarily_disconnect_edge( the_right_edge)
      the_right_branch_atom = (v1 in the_right_mol) and v1 or v2
      the_right_mol = mol.get_induced_subgraph_from_vertices( the_right_mol)
      the_right_branch = mol.get_induced_subgraph_from_vertices( the_right_branch)
      return (the_right_edge,
              the_right_mol,
              the_right_branch_atom,
              the_right_branch)
    # Handle symmetric molecules like benzene - fallback to simple disconnection
    # If no suitable edge found, try the first edge as fallback
    if mol.edges:
      e = list(mol.edges)[0]  # Get the first edge
      mol.temporarily_disconnect_edge(e)
      return e, mol, None, None
    else:
      raise Exception("Unable to generate SMILES: molecule has no edges")




  def disconnect_something_simple( self, mol, start_from=None):
    """returns (broken edge, resulting mol, atom where mol was disconnected, disconnected branch)"""
    # we cannot do much about this part
    if not mol.is_connected():
      print("unconnected ", mol)
    if start_from and start_from.get_degree() > 1:
      e = start_from._neighbors.keys()[0]
      mol.disconnect_edge( e)
      ps = [i for i in mol.get_connected_components()]
      if len( ps) == 1:
        return e, mol, None, None
      else:
        p1 = mol.get_induced_subgraph_from_vertices( ps[0])
        p2 = mol.get_induced_subgraph_from_vertices( ps[1])
        v1, v2 = e.get_vertices()
        p = (start_from in p1.vertices) and p1 or p2
        px = (p1 == p) and p2 or p1
        v = (v1 in p.vertices) and v1 or v2
        return e, p, v, px

    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if (d1 > 2 or d2 > 2):
        mol.disconnect_edge( e)
        ps = [i for i in mol.get_connected_components()]
        if len( ps) == 1:
          return e, mol, None, None
        else:
          p1 = mol.get_induced_subgraph_from_vertices( ps[0])
          p2 = mol.get_induced_subgraph_from_vertices( ps[1])
          v1, v2 = e.get_vertices()
          if start_from:
            p = (start_from in p1.vertices) and p1 or p2
            px = (p1 == p) and p2 or p1
            v = (v1 in p.vertices) and v1 or v2
            return e, p, v, px

          v = (v1 in p1.vertices) and v1 or v2
          return e, p1, v, p2
    print(mol, mol.is_connected(), ',', map( len, mol.get_connected_components()), ',', start_from)
    raise Exception("fuck, how comes!?")

  @staticmethod
  def _create_ring_join_smiles( index):
    i = index +1
    if i > 9:
      return "%%%d" % i
    else:
      return str( i)



def is_line( mol):
  """all degrees are 2 except of two with degree 1"""
  if len( mol.vertices) == 1:
    return True
  ones = 0
  for v in mol.vertices:
    d = v.get_degree() 
    if d == 1:
      if ones == 2:
        return False
      ones += 1
    elif d != 2:
      return False
  if ones == 2:
    return True
  return False

def is_pure_ring( mol):
  return filter( lambda x: x.get_degree() != 2, mol.vertices) == []

def match_atom_lists( l1, l2):
  """sort of bubble sort with counter"""
  count = 0
  for i1, v1 in enumerate( l1):
    for i2 in range( i1, len( l2)):
      v2 = l2[i2]
      if v2 == v1:
        if i1 != i2:
          l2[i1],l2[i2] = l2[i2],l2[i1]
          count += 1
        break
  assert l1 == l2
  return count


##################################################
## MODULE INTERFACE - newstyle

from .converter_base import converter_base

class smiles_converter( converter_base):

  # standard converter attrs
  reads_text = True
  writes_text = True
  reads_files = True
  writes_files = True

  default_configuration = {"R_GENERATE_COORDS": True,
                           "R_BOND_LENGTH": 1,
                           "R_LOCALIZE_AROMATIC_BONDS": True,
                           "R_EXPLICIT_HYDROGENS_TO_REAL_ATOMS": False,
                           
                           "W_AROMATIC_BOND_AUTODETECT": True,
                           "W_INDIVIDUAL_MOLECULE_SEPARATOR": ".",
                           "W_DETECT_STEREO_FROM_COORDS": True,
                           }

  def __init__( self):
    converter_base.__init__( self)

  def mols_to_text( self, structures):
    converter_base.mols_to_text( self, structures)
    sm = smiles()
    ret = []
    for mol in structures:
      if self.configuration["W_AROMATIC_BOND_AUTODETECT"]:
        mol.mark_aromatic_bonds()
      if self.configuration["W_DETECT_STEREO_FROM_COORDS"]:
        mol.detect_stereochemistry_from_coords()
      ret.append( sm.get_smiles( mol))
    self.last_status = self.STATUS_OK
    return self.configuration["W_INDIVIDUAL_MOLECULE_SEPARATOR"].join( ret)

  def read_text( self, text):
    # reaction support
    if ">" in text:
      if text.count(">") != 2:
        raise oasa_exceptions.oasa_smiles_error( "Only two '>' characters are supported in SMILES reactions. %d present" % text.count(">"))
      parts = text.split( ">")
      res = []
      for part in parts:
        res.append( self._read_string( part))
        if not self.last_status == self.STATUS_OK:
          res = []
          break
      if not res:
        self.result = []
        return []
      else:
        react = reaction.reaction()
        for p in res[0]:
          react.reactants.append( reaction.reaction_component( p))
        for p in res[2]:
          react.products.append( reaction.reaction_component( p))
        for p in res[1]:
          react.reagents.append( reaction.reaction_component( p))
        self.result = [react]
        return self.result
    else:
      return self._read_string( text)

  def _read_string( self, text):
    converter_base.read_text( self, text)
    sm = smiles()
    sm.read_smiles( text, explicit_hydrogens_to_real_atoms=self.configuration['R_EXPLICIT_HYDROGENS_TO_REAL_ATOMS'])
    mol = sm.structure
    if mol == None:
      return []
    mol.remove_zero_order_bonds()
    mols = mol.get_disconnected_subgraphs()
    for mol in mols:
      if self.configuration["R_LOCALIZE_AROMATIC_BONDS"]:
        mol.localize_aromatic_bonds()
        for b in mol.bonds:
          b.aromatic = 0
      if self.configuration["R_GENERATE_COORDS"]:
        coords_generator.calculate_coords( mol, bond_length=self.configuration['R_BOND_LENGTH'])
    self.result = mols
    self.last_status = self.STATUS_OK
    return mols

  def mols_to_file( self, structures, f):
    converter_base.mols_to_file( self, structures, f)
    f.write( self.mols_to_text( structures))

  def read_file( self, f):
    converter_base.read_file( self, f)
    mols = []
    for line in f:
      mol = self.text_to_mols( line)
      mols.extend( mol)
    return mols
    
converter = smiles_converter

# END OF MODULE INTERFACE
##################################################

##################################################
## MODULE INTERFACE - oldstyle

from . import coords_generator

reads_text = True
writes_text = True
reads_files = True
writes_files = True

def mol_to_text( structure):
  sm = smiles()
  structure.mark_aromatic_bonds()
  # Remove unimportant hydrogens for hydrogen-free SMILES
  if hasattr(structure, 'remove_unimportant_hydrogens'):
    structure.remove_unimportant_hydrogens()
  return sm.get_smiles( structure)

def text_to_mol( text, calc_coords=1, localize_aromatic_bonds=True):
  sm = smiles()
  sm.read_smiles( text)
  mol = sm.structure
  if localize_aromatic_bonds:
    mol.localize_aromatic_bonds()
    for b in mol.bonds:
      b.aromatic = 0
  if calc_coords:
    coords_generator.calculate_coords( mol, bond_length=calc_coords)
  return mol

def mol_to_file( mol, f):
  f.write( mol_to_text( mol))
  
def file_to_mol( f):
  return text_to_mol( f.read())


# END OF MODULE INTERFACE
##################################################



##################################################
# DEMO

if __name__ == '__main__':

  import sys

  def main( text, cycles):
    t = time.time()
    conv = converter()
    conv.configuration['R_EXPLICIT_HYDROGENS_TO_REAL_ATOMS'] = True
    for j in range( cycles):
      mols = conv.read_text( text)
      for mol in mols:
        #mol.remove_unimportant_hydrogens()
        print("  summary formula:   ", mol.get_formula_dict())
      text = conv.mols_to_text( mols)
      print("  generated SMILES:   %s" % text)
      print("  --")
    t = time.time()-t
    print('time per cycle', round( 1000*t/cycles, 2), 'ms')

  repeat = 3

  if not len( sys.argv) > 1:
    text = "COc5ccc4c2sc(cc2nc4c5)-c(cc1nc3c6)sc1c3ccc6OC"  #"ccc4ccc2cc1cc3ccccc3cc1cc2c4"
    text = "[C@@H](C)(CC)O"
    text = "O=C1c2c3c4c(cc2)c2ccc5c6c2c(ccc6C(=O)c2c5cccc2)c4ccc3c2ccccc12"
  else:
    text = sys.argv[1]

  print("oasa::SMILES DEMO")
  print("converting following smiles to smiles (%d times)" % repeat)
  print("  starting with:      %s" % text)
  print("  --------------------")
  main( text, repeat)

# DEMO END
##################################################






##################################################
# TODO

# last branch does not need to be branch
# optimize for either speed or human-readability
# at first get rid of the edges common to two or more rings, not critical
# handling of the start_from in disconnect_something as in disconnect_something_simple
# could in disconnect_something happen that a start_from will be returned in a branch?
# if I do "print m.get_smiles( m.structure)" the ring counting does not work after that
## the transformation can be destructive!!!  - check

## THIS IS A PROBLEM : C=1ccC=2C=1C=CC=CC=2  (should be azulene)

# E/Z stereo is ignored in rings
