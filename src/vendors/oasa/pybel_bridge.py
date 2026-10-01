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

"""This module is intended for integration of OpenBabel/Pybel into OASA.
It provides functionality for translation of OASA molecules into Pybel molecules."""

import pybel, openbabel
from .molecule import molecule
from .atom import atom
from .bond import bond
from .periodic_table import periodic_table as PT
num2symbol = dict( [(v['ord'],k) for k,v in PT.items()])


class PybelConverter( object):

  ## -------------------- Pybel to OASA --------------------

  @classmethod
  def pybel_to_oasa_atom( self, patom):
    oatom = atom()
    psymbol = num2symbol[patom.atomicnum]
    oatom.symbol = psymbol
    oatom.x, oatom.y, oatom.z = patom.coords
    oatom.charge = patom.formalcharge
    return oatom

  @classmethod
  def pybel_to_oasa_bond( self, pbond):
    obond = bond()
    obond.order = pbond.GetBondOrder()
    return obond

  @classmethod
  def pybel_to_oasa_molecule_with_atom_map( self, pmol):
    omol = molecule()
    patom_idx2oatom = {}
    for pa in pmol.atoms:
      oa = self.pybel_to_oasa_atom( pa)
      omol.add_vertex( oa)
      patom_idx2oatom[ pa.idx] = oa
    for pb in openbabel.OBMolBondIter( pmol.OBMol):
      ob = self.pybel_to_oasa_bond( pb)
      i1 = pb.GetBeginAtomIdx()
      i2 = pb.GetEndAtomIdx()
      oa1 = patom_idx2oatom[ i1]
      oa2 = patom_idx2oatom[ i2]
      omol.add_edge( oa1, oa2, ob)
    return omol, patom_idx2oatom

  @classmethod
  def pybel_to_oasa_molecule( self, *args):
    omol, patom_idx2oatom = self.pybel_to_oasa_molecule_with_atom_map( *args)
    return omol

  ## // -------------------- Pybel to OASA --------------------

  ## -------------------- OASA to Pybel --------------------

  @classmethod
  def oasa_to_ob_atom( self, oatom):
    num = PT[ oatom.symbol]['ord']
    obatom = openbabel.OBAtom()
    obatom.SetAtomicNum( num)
    obatom.SetVector( oatom.x, oatom.y, oatom.z)
    obatom.SetFormalCharge( oatom.charge)
    #patom = pybel.Atom()
    #patom.atomicnum = num
    return obatom

  @classmethod
  def oasa_to_ob_bond( self, obond):
    pbond = openbabel.OBBond()
    pbond.SetBondOrder( obond.order)
    return pbond

  @classmethod
  def oasa_to_pybel_molecule_with_atom_map( self, omol):
    obmol = openbabel.OBMol()
    oatom2obatom_idx = {}
    for i,oa in enumerate( omol.atoms):
      oba = self.oasa_to_ob_atom( oa)
      oba.SetIdx( i+1)
      obmol.AddAtom( oba)
      oatom2obatom_idx[ oa] = i+1
    for ob in omol.bonds:
      oa1, oa2 = ob.vertices
      idx1 = oatom2obatom_idx[ oa1]
      idx2 = oatom2obatom_idx[ oa2]
      obmol.AddBond( idx1, idx2, ob.order)
    pmol = pybel.Molecule( obmol)
    return pmol, oatom2obatom_idx

  @classmethod
  def oasa_to_pybel_molecule( self, *args):
    omol, patom_idx2oatom = self.oasa_to_pybel_molecule_with_atom_map( *args)
    return omol

  ## // -------------------- OASA to Pybel --------------------

  ## -------------------- conversion --------------------

  @classmethod
  def read_text( self, format, text):
    """returns a list of OASA molecules from a string"""
    obc = openbabel.OBConversion()
    if not obc.SetInFormat( format):
        raise ValueError("invalid format %s" % format)
    obmol = openbabel.OBMol()
    ok = obc.ReadString( obmol, text)
    mols = []
    while ok:
      mol = self.pybel_to_oasa_molecule( pybel.Molecule( obmol))
      if len( mol.vertices) > 0:
        # throw away empty molecules - CDX seems to cause these to appear
        mols.append( mol)
      obmol = openbabel.OBMol()
      ok = obc.Read( obmol)
    return mols

  ## // -------------------- conversion --------------------

  ## -------------------- conversion support --------------------
  ## not needed - pybel provides this functionality, I just overlooked it
  
  @classmethod
  def get_supported_input_formats( self):
    conv = openbabel.OBConversion()
    ret = []
    for form in conv.GetSupportedInputFormat():
      form_id, desc = [x.strip() for x in form.split("--")]
      ret.append( (form_id, desc))
    return ret

  @classmethod
  def get_supported_output_formats( self):
    conv = openbabel.OBConversion()
    ret = []
    for form in conv.GetSupportedOutputFormat():
      form_id, desc = [x.strip() for x in form.split("--")]
      ret.append( (form_id, desc))
    return ret

  ## // -------------------- conversion support --------------------

class ForceFieldOptimizer( object):
  """provides the forcefield interface for oasa molecules"""
  
  def __init__( self, mol, forcefield_name="Ghemical"):
    self.mol = mol
    self._amol, self._oatom2patom_idx = PybelConverter.oasa_to_pybel_molecule_with_atom_map( self.mol)
    self._patom_idx2oatom = dict( [(v,k) for k,v in self._oatom2patom_idx.items()])
    self.forcefield_name = forcefield_name
    self.ff = openbabel.OBForceField.FindForceField( self.forcefield_name)
    self.ff.Setup( self._amol.OBMol)
    
  def conjugate_gradients( self, step_size=10, max_steps=100):
    """partial result will be yielded after the 'step_size' number of iteration has
    been performed. When false is yielded, we did not reach convergence criteria yet.
    'max_steps' is the number of repetitions of this, so the real number of iterations
    is 'max_steps' * 'step_size'
    """
    self.ff.ConjugateGradientsInitialize( step_size * max_steps)
    for i in range( max_steps):
      do_continue = self.ff.ConjugateGradientsTakeNSteps( step_size)
      yield not do_continue
      if not do_continue:
        break

  def update_coords( self):
    self.ff.UpdateCoordinates( self._amol.OBMol)
    for patom in self._amol:
      oatom = self._patom_idx2oatom[ patom.idx]
      oatom.coords = patom.coords
    


if __name__ == "__main__":
  import time
  
  pmol = pybel.readstring("smi", "CC(=O)O")
  omol = PybelConverter.pybel_to_oasa_molecule( pmol)
  print(omol)
  from . import smiles
  c = smiles.converter()
  print(c.mols_to_text( [omol]))

  #print get_supported_input_formats()
  #print get_supported_output_formats()

  if False:
    from . import molfile
    f = file( "tbu-benzen.mol", "r")
    mol = molfile.file_to_mol( f)
    f.close()

    print(PybelConverter.oasa_to_pybel_molecule( mol))

    t = time.time()
    ff = ForceFieldOptimizer( mol)
    gen = ff.conjugate_gradients()
    for step in gen:
      print(step, ff.ff.Energy())
      ff.update_coords()
      print([a.coords for a in mol.atoms])
    print(time.time() - t)

  a = pybel.readstring( "smi", "CCCC\nCCC")

  for mol in PybelConverter.read_text( "smi", "CCCC\nCCC"):
    print(mol)

  print("AAAA")
  t = file( "/zaloha/temp/3-chlorbenzoic_acid.cdx", "r").read()
  pmol = pybel.readstring( "cdx", t)
  print("XXX", pmol)

  for mol in PybelConverter.read_text( "cdx", t):
    print(mol)

