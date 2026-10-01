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
from .molecule import molecule
from . import periodic_table as pt
from . import molfile
from . import misc

import re

import operator
import time
import xml.dom.minidom as dom
from . import dom_extensions
import string
import os
from . import coords_generator
from .oasa_exceptions import oasa_not_implemented_error, oasa_inchi_error, oasa_unsupported_inchi_version_error
import select
import sys
from .stereochemistry import cis_trans_stereochemistry

class inchi( plugin):

  name = "inchi"
  read = 1
  write = 0

  proton_donors = ['O','F','Cl','Br','I','N','S']
  proton_acceptors = ['P','N','S','O']
  electron_donors = ['O','N','S','P']
  electron_acceptors = ['B']

  stereo_remap = {"+": cis_trans_stereochemistry.OPPOSITE_SIDE,
                  "-": cis_trans_stereochemistry.SAME_SIDE,
                  "?": cis_trans_stereochemistry.UNDEFINED,
                  "u": cis_trans_stereochemistry.UNDEFINED,
                  }


  def __init__( self, structure=None):
    self.structure = structure
    self._no_possibility_to_improve = False
    self.cleanup()

  def cleanup( self):
    self.charge = 0
    self.forced_charge = 0
    self._protonation_dealt_with_already = 0
    self._added_hs = set()

  def set_structure( self, structure):
    self.structure = structure

  def get_structure( self):
    return self.structure



  def split_layers( self, text):
    # try if we have a xml document in hand
    try:
      doc = dom.parseString( text)
    except:
      # it seems not to be the case
      if text.startswith( "InChI="):
        text = text[6:]
      layers = re.split( "/", text.strip())
      # fixed hydrogen layer is marked stupidly by f/h :(, we have to deal with it
      i = None
      for l in layers:
        if l.startswith("f"):
          i = layers.index( l)
          break
      if i != None:
        del layers[i]
        layers[i] = 'f'+layers[i][1:]
      return layers
    structs = doc.getElementsByTagName( 'structure')
    if not structs:
      raise Exception("no structures found in xml string %s" % text)
    struct = structs[0]
    layers = []
    for name in ('version','formula','connections','H'):
      try:
        layer = struct.getElementsByTagName( name)[0]
      except IndexError:
        raise Exception("no '%s' tag found in xml string" % name)
      layers.append( dom_extensions.getTextFromElement( layer))
    return layers



  def get_layer( self, prefix):
    for l in self.layers:
      if l.startswith( prefix):
        return l[1:]

  def get_atom_with_inchi_number( self, number):
    assert type( number) == int
    for v in self.structure.vertices:
      if v.properties_['inchi_number'] == number:
        return v

  def read_inchi( self, text):
    try:
      self._read_inchi( text)
    except AssertionError:
      raise oasa_inchi_error( "Localization of bonds, charges or movable hydrogens failed")
    except:
      raise
    

  def _read_inchi( self, text):
    if not text:
      raise oasa_inchi_error( "No inchi was given")
    self.structure = molecule()
    self.layers = self.split_layers( text)
    # version support (very crude)
    self.version = self._check_version( self.layers[0])
    if not self.version:
      raise oasa_unsupported_inchi_version_error( self.layers[0])
    elif str( self.version[0]) != '1' or str( self.version[1]) != '0':
      raise oasa_unsupported_inchi_version_error( self.layers[0])
    
    self.hs_in_hydrogen_layer = self.get_number_of_hydrogens_in_hydrogen_layer()
    self.read_sum_layer()
    self.read_connectivity_layer()
    self._charge_mover = self._move_charge_somewhere_else()
    repeat = True
    run = 0
    # we have to repeat this step in order to find the right positioning of movable hydrogens
    while repeat and not self._no_possibility_to_improve:
      # cleanup
      for h in self._added_hs:
        self.structure.remove_vertex( h)
      for b in self.structure.edges:
        b.order = 1
      for v in self.structure.vertices:
        v.symbol = v.symbol
        v.charge = 0
      self.cleanup()

      # the code itself
      run += 1
      assert run < 50
      self._deal_with_notorious_groups()
      self.process_forced_charges()
      self.read_hydrogen_layer( run=run)
      self.read_charge_layer()
      self.read_p_layer()
      self.deal_with_da_bonds()
      #self._deal_with_valencies()
      self.compensate_for_forced_charges()
      self.structure.add_missing_bond_orders()
      #self.read_double_bond_stereo_layer()

      # here we check out if the molecule seems ok
      fvs = [v for v in self.structure.vertices if v.free_valency]
      if not fvs and not filter( None, [not v.order for v in self.structure.edges]):
        repeat = False
      else:
        if len( fvs) == 1:
          a = fvs[0]
          a.symbol = a.symbol # this resets the valency
          a.raise_valency_to_senseful_value()
          if not a.free_valency:
            repeat = False
          # deal with CO and other that require adding both + and -
          elif a.free_valency == 2 and [n for n in a.neighbors if n.free_valency==0 and n.symbol in self.electron_donors]:
            neighbor = [n for n in a.neighbors if n.free_valency==0 and n.symbol in self.electron_donors][0]
            a.charge -= 1
            neighbor.charge += 1
            a.get_edge_leading_to( neighbor).order = 1
            self.structure.add_missing_bond_orders()
            if not [v for v in self.structure.vertices if v.free_valency > 0]:
              # we fixed it
              repeat = False
            

      if repeat and self._no_possibility_to_improve and self.charge:
        try:
          self.next(_charge_mover)
        except StopIteration:
          pass
        else:
          self._no_possibility_to_improve = False
          run = 0

    if repeat and self._no_possibility_to_improve:
##       if len( filter( None, [v.free_valency for v in self.structure.vertices])) == 1:
##         print
##         print [(v.symbol, v.valency, v.free_valency)  for v in self.structure.vertices if v.free_valency], filter( None, [not v.order for v in self.structure.edges]), text
##       if sum( [v.charge for v in self.structure.vertices]) != self.charge:
##         print "Charge problem", sum( [v.charge for v in self.structure.vertices]), self.charge
      #pass
      raise oasa_inchi_error( "Localization of bonds, charges or movable hydrogens failed")
#    print >> sys.stderr, "run:", run



  def read_sum_layer( self):
    if "." in self.layers[1]:
      raise oasa_not_implemented_error( "INChI", "multiple compound systems are not supported by the library")

    form = pt.formula_dict( self.layers[1])
    processed_hs = 0 #for diborane and similar compounds we must process some Hs here
    j = 0
    for k in form.sorted_keys():
      for i in range( form[k]):
        if k == 'H':
          # we want to process only the Hs that are not in the h-layer
          if processed_hs >= form[k] - self.hs_in_hydrogen_layer:
            continue
          else:
            processed_hs += 1
        j += 1
        a = self.structure.create_vertex()
        a.symbol = k
        self.structure.add_vertex( a)
        a.properties_['inchi_number'] = j



  def read_connectivity_layer( self):
    layer = self.get_layer( "c")
    if not layer:
      return 
    chunks = re.split( "([0-9]*)", layer)
    chunks = filter( None, chunks)
    chunks = filter( lambda x: x!='-', chunks)
    last_atom = None
    bracket_openings = []
    for c in chunks:
      if c == '(':
        bracket_openings.append( last_atom)
      elif c == ')':
        last_atom = bracket_openings.pop(-1)
      elif c == ",":
        last_atom = bracket_openings.pop(-1)
        bracket_openings.append( last_atom)
      else:
        try:
          i = int( c)
        except:
          raise ValueError("unexpected character %s in the connectivity layer" % c)
        # atom
        if last_atom:
          self.structure.add_edge( last_atom-1, i-1)
        last_atom = i


  def read_hydrogen_layer( self, run=0):
    layer = self.get_layer( "h")
    # check presence of the layer
    if not layer:
      return

    fixed_layer = self.get_layer( "f")
    if fixed_layer:
      # we can skip movable hydrogens
      self._read_simple_hydrogen_layer( fixed_layer)
    else:
      re_for_brackets = r"\([H\d,\-]+?\)"
      brackets = re.findall( re_for_brackets, layer)
      for bracket in brackets:
        self._process_moving_hydrogen( bracket[1:-1], run=run)

    # we improve only in the _process_moving_hydrogen so if it was not called there is no possibility for improvement
    if not self._added_hs or fixed_layer:
      self._no_possibility_to_improve = True
      
    self._read_simple_hydrogen_layer( layer)



  def _read_simple_hydrogen_layer( self, layer):
    """just takes the layer and adds hydrogens according to what it seas,
    it does not care about moving hydrogens and stuff"""

    re_for_brackets = r"\([H\d,\-]+?\)"
    layer = re.sub( re_for_brackets, "", layer)  # clean the brackets out

    for vs, num in self._parse_h_layer( layer):
      for v in vs:
        for j in range( num):
          h = self.structure.create_vertex()
          h.symbol = 'H'
          self.structure.add_vertex( h)
          self.structure.add_edge( v-1, h)
          self._added_hs.add( h)



  def read_charge_layer( self):
    layer = self.get_layer( "q")
    if not layer:
      return

    charge = int( layer[1:]) - self.forced_charge
    self.charge = int( layer[1:])


    if charge > 0:
      # at first we try to put positive charge on atoms with negative free_valency
      change = True
      steps = 0
      while charge and change:
        steps += 1
        assert steps < 10
        change = False
        for v in self.structure.vertices:
          if v.free_valency == -1:
            if charge > 0:
              v.charge = -v.free_valency
            else:
              v.charge = v.free_valency
            charge -= v.charge
            change = True
            break

    if charge < 0:
      # we try to put negative charge on orphan atoms
      steps = 0
      change = True
      while charge and change:
        steps += 1
        assert steps < 10
        change = False
        for v in self.structure.vertices:
          if v.free_valency and sum( [n.free_valency for n in v.neighbors]) < v.free_valency:
            charge = self._valency_to_charge( v, charge)
            change = True
            break

    # this could be rewritten to put charges to places where a segment with odd number of free_valencies is
    # this would be more general than counting the free_valencies for the whole molecule
    free_valencies = sum( [a.free_valency for a in self.structure.vertices])
    if free_valencies % 2:
      # odd number of free_valencies means we have to put the charge somewhere to make a free_valency there
      change = True
      steps = 0

      while charge and change:
        change = False
        steps += 1
        assert steps < 10

        for v in self.structure.vertices:
          if v.symbol != 'C' and not v.free_valency and filter( None, [n.free_valency for n in v.neighbors]):
            v.charge = charge
            charge = 0
            change = True
            break

    # then we do the others, at first trying heteroatoms
    steps = 0
    change = True
    while charge and change:
      change = False
      steps += 1
      assert steps < 10

      for v in self.structure.vertices:
        if v.free_valency and v.symbol != 'C':
          charge = self._valency_to_charge( v, charge)
          change = True
          break


    # then we try the carbon atoms as well
    steps = 0
    change = True
    while charge and change:
      steps += 1
      assert steps < 10

      change = False
      for v in self.structure.vertices:
        if v.free_valency:
          charge = self._valency_to_charge( v, charge)
          change = True
          break

    if charge:
      raise oasa_exceptions.oasa_inchi_error( "The molecular charge could not be allocated to any atom (%d)." % charge)



  def read_p_layer( self):
    layer = self.get_layer( "p")
    if not layer:
      return


    p = int( layer[1:]) - self._protonation_dealt_with_already
    self.charge += p
    charges = self.forced_charge
    old_p = p

    # at first add protons to negative charges
    for v in self.structure.vertices:
      if p > 0:
        if v.charge < 0:
          h = self.structure.create_vertex()
          h.symbol = 'H'
          self.structure.add_vertex( h)
          self.structure.add_edge( v, h)
          p -= 1
          v.charge += 1
          self._protonation_dealt_with_already += 1

    if self._protonation_dealt_with_already:
      # p was solved above
      dp = self._protonation_dealt_with_already - sum( [v.charge for v in self.structure.vertices])
      if dp > 0:
        # there were no forced charges, so we are not in a zwittrion state, huh, who does not understand that?!
        for v in self.structure.vertices:
          if dp > 0 and v.symbol in self.proton_acceptors:
            v.charge += 1
            dp -= 1

    # then come the other additional protones
    while p:
      old_p = p
      for v in self.structure.vertices:
        if p > 0:
          if v.symbol in self.proton_acceptors:
            h = self.structure.create_vertex()
            h.symbol = 'H'
            self.structure.add_vertex( h)
            self.structure.add_edge( v, h)
            self._added_hs.add( h)
            p -= 1
            if charges:
              # the charge is already present as a forced charge
              charges -= 1
            else:
              # the charge can normaly be delt with
              v.charge += 1
            break
        elif p < 0:
          if v.symbol in self.proton_donors:
            hs = [n for n in v.neighbors if n.symbol=='H']
            if hs:
              h = hs[0]
              self.structure.remove_vertex( h)
              p += 1
              if charges:
                # the charge is already present as a forced charge
                charges += 1
              else:
                # the charge can normaly be delt with
                v.charge -= 1
              break
      assert p < old_p


  def read_double_bond_stereo_layer( self):
    def get_lowest_numbered_neighbor( atom, verbotten_atom):
      neighs = [(n.properties_.get('inchi_number',1000000), n)
                for n in atom.neighbors
                if n is not verbotten_atom and "inchi_number" in n.properties_]
      if neighs:
        neighs.sort( reverse=True)
        return neighs[0][1]
      else:
        # we take the only one that is there (if any)
        neighs = [n for n in atom.neighbors if n is not verbotten_atom]
        if neighs:
          return neighs[0]
        else:
          raise oasa_inchi_error( "No neigbors on atom with stereo information!")
        

    
    layer = self.get_layer( "b")
    if not layer:
      return
    for a1,a2,sign in re.findall( r"(\d+)-(\d+)([-+?u])", layer):
      atom1 = self.get_atom_with_inchi_number( int( a1))
      atom2 = self.get_atom_with_inchi_number( int( a2))
      bond = self.structure.get_edge_between( atom1, atom2)
      # we need neighbors with lowest inchi_number
      neigh1 = get_lowest_numbered_neighbor( atom1, atom2)
      neigh2 = get_lowest_numbered_neighbor( atom2, atom1)
      st = cis_trans_stereochemistry( references=[neigh1,atom1,atom2,neigh2],
                                      value=self.stereo_remap[sign])
      self.structure.add_stereochemistry( st)



  def _valency_to_charge( self, v, charge):
    """sets charge of the atom so that it has minimal free_valency,
    returns the 'unused' charge in case it is not comsumed whole"""
    ch = misc.signum( charge) * min( abs( charge), v.free_valency)
    if charge < 0:
      charge += ch
      v.charge = -ch
    else:
      charge -= ch
      v.charge = ch
    return charge
        

  def get_number_of_hydrogens_in_hydrogen_layer( self):
    # version check
    layer = self.get_layer( "h")
    if not layer:
      return 0

    # check if we can handle it
    if "*" in layer or ";" in layer:
      raise oasa_not_implemented_error( "INChI", "multiple compound systems are not supported by the library")

    ret = 0

    re_for_brackets = r"\([H\d,\-]+?\)"
    brackets = re.findall( re_for_brackets, layer)
    for bracket in brackets:
      ret += self._get_hs_in_moving_hydrogen( bracket[1:-1])
    layer = re.sub( re_for_brackets, "", layer)  # clean the brackets out

    for vs, num in self._parse_h_layer( layer):
      ret += len( vs) * num

    return ret


  def _get_hs_in_moving_hydrogen( self, chunk):
    chks = chunk.split( ',')
    if len( chks[0]) > 1:
      if chks[0][1:].endswith( "-"):
        hs = int( chks[0][1:-1] or 1)  # number of atoms 
      else:
        hs = int( chks[0][1:])  # number of atoms
    else:
      hs = 1
    return hs
    


  def _check_version( self, ver):
    """returns a tuple of two version numbers"""
    if "Beta" in ver:
      return None
    v = ver.strip(string.ascii_lowercase+string.ascii_uppercase)
    if "." in v:
      return v.split('.')
    else:
      return v, 0


  def _process_moving_hydrogen( self, chunk, run=0):
    chks = chunk.split( ',')
    if len( chks[0]) > 1:
      if chks[0][1:].endswith("-"):
        hs = int( chks[0][1:-1] or 1) + 1  # number of atoms 
        self._protonation_dealt_with_already += 1
      else:
        hs = int( chks[0][1:])  # number of atoms
    else:
      hs = 1

    vertices = [self.structure.vertices[ int( i)-1] for i in chks[1:]]
    vs = vertices
    take = hs

    # here we generate still shorter combinations of vertices that receive hydrogens
    # because of applying them in circles in the next step, we simulate the case
    # where more than one hydrogen is put to one atom
    while (run > misc.x_over_y( len( vs), take)) and take:
      run -= misc.x_over_y( len( vs), take)
      take -= 1
    if not take:
      self._no_possibility_to_improve = True
      return []
    variations = misc.gen_variations( vs, take)
    for i in range( run):
      vs = next(variations)
        
    while hs:
      for v in vs:
        h = self.structure.create_vertex()
        h.symbol = 'H'
        self.structure.add_vertex( h)
        self.structure.add_edge( v, h)
        self._added_hs.add( h)
        hs -= 1
        if not hs:
          break



  def _split_h_layer( self, layer):
    was_h = False
    chunks = []
    chunk = ""
    for c in layer:
      if c == 'H':
        was_h = True
        chunk += c
      elif c == "," and was_h:
        was_h = False
        chunks.append( chunk)
        chunk = ""
      else:
        chunk += c
    if chunk:
      chunks.append( chunk)
    return chunks


  def _parse_h_layer( self, layer):
    chunks = self._split_h_layer( layer)
    for chunk in chunks:
      try:
        head, tail = chunk.split( 'H')
      except Exception as e:
        raise oasa_inchi_error( "error in hydrogen layer - missing H symbol")
      num_h = tail and int( tail) or 1
      vertices = []
      for p in head.split( ","):
        if "-" in p:
          try:
            a, b = map( int, p.split("-"))
          except Exception as e:
            raise oasa_inchi_error( "error in hydrogen layer - non-number character(s) present in atom range specification")
          vertices.extend( range( a, b+1))
        else:
          try:
            vertices.append( int( p))
          except Exception as e:
            raise oasa_inchi_error( "error in hydrogen layer - non-number character(s) present in atom specification")

      yield vertices, num_h


  def process_forced_charges( self):
    """this marks the charges that are forced by the connectivity and thus helps
    process zwitrions"""
    forced_charge = 0
    for v in self.structure.vertices:
      if v.symbol in ("N","S") and v.free_valency == -1:
        v.charge = 1
        forced_charge += 1
    self.forced_charge = forced_charge


  def compensate_for_forced_charges( self):
    """if there were foced charges and the molecule should not have any charge,
    we have to take care of it here"""
    charge = self.charge - self.forced_charge
    old_charge = charge
    while not self.charge and charge:
      for chunk in self.structure._gen_free_valency_connected_components():
        if len( chunk) % 2:
          if len( chunk) == 1:
            vs = chunk
          else:
            vs = [v for v in chunk if len( set( v.neighbors) & chunk) == 1]
          for v in vs:
            if charge < 0:
              if v.symbol in self.proton_donors and v.free_valency:
                v.charge -= 1
                charge += 1
            elif charge > 0:
              if v.symbol in self.proton_acceptors:
                v.charge += 1
                charge -= 1
            if not charge:
              return
##       if old_charge == charge:
##         print "AAAAA", self.layers
      assert old_charge != charge
                
        
  def deal_with_da_bonds( self):
    """deal with donor-acceptor bonds, this fixes mostly boron containing compounds"""
    for v in self.structure.vertices:
      if v.free_valency < 0 and v.symbol in self.electron_acceptors:
        ns = [n for n in v.neighbors if n.symbol in self.electron_donors]
        if ns:
          n = ns[0]
          n.charge += 1
          v.charge -= 1

  def _move_charge_somewhere_else( self):
    if self.charge < 0:
      pass
    elif self.charge == 1:
      for v in self.structure.vertices:
        if v.symbol in self.proton_acceptors:
          v.charge = 1
          yield None

          
  def _deal_with_notorious_groups( self):
    """some groups such as NO2, SO3H etc. need the valency of the central atom to be risen,
    this is done here"""
    for v in self.structure.vertices:
      if v.symbol == "N" and v.degree > 2:
        for n in v.neighbors:
          if n.symbol == "O" and n.degree == 1:
            v.raise_valency()
            break
      elif v.symbol == "S" and v.symbol > 2:
        for n in v.neighbors:
          if n.symbol == "O" and n.degree == 1:
            v.raise_valency()


  def _deal_with_valencies( self):
    go = True
    i = 0
    while go:
      i += 1
      if i> 100:
        print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        break
      go = False
      # set all valency_raise to zero
      for v in self.structure.vertices:
        v.properties_['valency_raise'] = 0
      # mark vertices with the demands to raise valency
      for comp in self.structure._gen_free_valency_connected_components():
#        print map( str, comp)
        if sum( [v.free_valency for v in comp]) % 2:
          go = True
          for v in comp:
            for ne in v.neighbors:
              if ne not in comp:
                ne.properties_['valency_raise'] += 1
      # raise the valency of the most asked for atom
      vs = [(v.properties_['valency_raise'], v) for v in self.structure.vertices]
      vs.sort()
      vs.reverse()
      for (x, v) in vs:
#        print x, v, v.valency, v.free_valency
        if x and v.raise_valency():
          break
##           if v.free_valency == 0 and v.symbol in self.proton_acceptors:
##             if v.charge == 0:
##               v.charge = 1
##             elif v.charge == 1:
##               if v.raise_valency():
##                 v.charge = 0
##                 break
##           elif v.raise_valency():
##             break
#          print x, v, v.symbol, v.degree, v.valency, v.free_valency
#          print [e.order for e in v.neighbor_edges]
#          break
#        elif x == 0:
#          break
        
import subprocess

def _run_command( command, input):
  p = subprocess.Popen( command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  o, er = p.communicate( str( input))
  return o


def generate_inchi_and_inchikey( m, program=None, fixed_hs=True, ignore_key_error=False):
  """ignore the case when InChIKey cannot be generated for some reason
  (no mhash library and old InChI program)"""
  if not program:
    from . import config
    program = config.Config.inchi_binary_path
  mf = molfile.mol_to_text( m)
  if os.name == 'nt':
    options = "/AUXNONE /STDIO /Key" + (fixed_hs and " /FixedH" or "")
  else:
    options = "-AUXNONE -STDIO -Key" + (fixed_hs and " -FixedH" or "")
  #print options
  command = [os.path.abspath( program)] + options.split()
  text = _run_command( command, mf)
  inchi = ""
  key = ""
  warnings = []
  for line in text.splitlines():
    if line.startswith( "Warning"):
      warnings.append(line.strip() + "\n")
    elif line.startswith( "End of file detected"):
      pass
    elif line.startswith( "InChIKey="):
      key = line.strip()[9:]
      break
    elif line.startswith( "InChI="):
      inchi = inchi + line.strip()
    elif line.startswith( "Error"):
      break
  if not inchi:
    raise oasa_inchi_error( "InChI program did not create any output InChI")
  if not key:
    # probably old version of the InChI software
    try:
      from . import inchi_key
    except ImportError:
      if ignore_key_error:
        key = None
      else:
        raise oasa_inchi_error( "InChIKey could not be generated - inchi_key module failed to load properly.")
    else:
      key = inchi_key.key_from_inchi( inchi)
  return inchi, key, warnings


def generate_inchi( m, program=None, fixed_hs=True):
  inchi, key, warnings = generate_inchi_and_inchikey( m, program=program, fixed_hs=fixed_hs, ignore_key_error=True)
  return inchi, warnings

def generate_inchi_key( m, program=None, fixed_hs=True):
  inchi, key, warnings = generate_inchi_and_inchikey( m, program=program, fixed_hs=fixed_hs)
  return key, warnings


  
##################################################
# MODULE INTERFACE

reads_text = 1
reads_files = 1
writes_text = 1
writes_files = 1

def text_to_mol( text, include_hydrogens=True, mark_aromatic_bonds=False, calc_coords=1):
  try:
    inc = inchi()
    inc.read_inchi( text)
    mol = inc.structure
    #mol.add_missing_bond_orders()
    if not include_hydrogens:
      mol.remove_unimportant_hydrogens()
    if calc_coords:
      coords_generator.calculate_coords( mol, bond_length=calc_coords)
    if mark_aromatic_bonds:
      mol.mark_aromatic_bonds()
    return mol
  except oasa_inchi_error:
    raise
  except oasa_not_implemented_error:
    raise
  except:
    raise

def file_to_mol( f):
  return text_to_mol( f.read())


def mol_to_text( mol, program=None, fixed_hs=True):
  return generate_inchi( mol, program=program, fixed_hs=fixed_hs)[0]


def mol_to_file( mol, f):
  f.write( mol_to_text( mol))



#
##################################################




##################################################
# DEMO

if __name__ == '__main__':

  from . import smiles

  def main( text, cycles):
    t1 = time.time()
    for jj in range( cycles):
      mol = text_to_mol( text, calc_coords=True, include_hydrogens=False)
      print(map( str, [b for b in mol.bonds if b.order == 0]))
      print("  smiles: ", smiles.mol_to_text( mol))
      print("  inchi:  ", generate_inchi( mol, fixed_hs=False, program="/home/beda/bin/stdinchi-1"))
      print("  charge: ", sum( [a.charge for a in mol.vertices]))
      print("  mw:     ", mol.weight)
    print(generate_inchi_and_inchikey( mol, fixed_hs=False, program="/home/beda/bin/stdinchi-1"))
    t1 = time.time() - t1
    print('time per cycle', round( 1000*t1/cycles, 2), 'ms')

  repeat = 3
  inch = "InChI=1S/C34H16O2/c35-33-25-7-3-1-5-17(25)19-9-11-21-24-14-16-28-32-20(18-6-2-4-8-26(18)34(28)36)10-12-22(30(24)32)23-13-15-27(33)31(19)29(21)23/h1-16H"

  #"InChI=1/C34H16O2/c35-33-25-7-3-1-5-17(25)19-9-11-21-24-14-16-28-32-20(18-6-2-4-8-26(18)34(28)36)10-12-22(30(24)32)23-13-15-27(33)31(19)29(21)23/h1-16H"
  #inch = "InChI=1S/CO/c1-2"
  print("oasa::INCHI DEMO")
  print("converting following inchi into smiles (%d times)" % repeat)
  print("  inchi:   %s" % inch)
  
  #import profile
  #profile.run( 'main( inch, repeat)')
  main( inch, repeat)


# DEMO END
##################################################



##################################################
## TODO

## rename layer to sublayer


##################################################

