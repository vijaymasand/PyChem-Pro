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

import unittest
from . import linear_formula
from . import smiles
from . import molecule

# helper functions

def create_test(i,name):
  def test(obj):
    getattr(obj, name)(i)
  return test


## linear formula testing

class TestLinearFormula(unittest.TestCase):

  formulas = [("CH3COOH","CC(=O)O",0,0),
              ("CH3C(CH3)3","CC(C)(C)C",0,0),
              ("CH2C(CH3)3","CC(C)(C)C",1,0),
              ("(CH2)7Cl","CCCCCCCCl",1,0),
              ]
    
  def _testformula(self, num):
    l = linear_formula.linear_formula()
    linear, smile, start_valency, end_valency = self.formulas[num]
    m1 = l.parse_text( linear, start_valency=start_valency, end_valency=end_valency) 
    m2 = smiles.text_to_mol( smile)
    self.assert_(molecule.equals(m1,m2,level=3))

# this creates individual test for formulas
for i in range( len( TestLinearFormula.formulas)):
  setattr( TestLinearFormula, "testformula"+str(i+1), create_test(i,"_testformula"))

## // linear formula testing



## substructure testing

class TestSubstructure(unittest.TestCase):

  formulas = [("CC(=O)C","CC(=O)C",True),
              ("CCCC","CCC",True),
              ("CC(=O)C","C=O",True),
              ("C=O","CC(=O)C",False),
              ("C(=O)H","C=O",True),
              ("C=O","C(=O)H",True),      # implicit hydrogens work
              ("C(=O)O","C(=O)OH",True),  # implicit hydrogens work
              ("C(=O)OH","C(=O)O",True),
              ("C(=O)OC","C(=O)OH",False),  # explicit hydrogens work
              ("C(=O)OC","C(=O)[OH]",False),  # explicit hydrogens in atom specs work
              ("C(=O)OC","C(=O)O",True),
              #
              ("CC(=O)H","C(=O)H",True),
              # charges
              ("C(=O)[O-]","C(=O)O",True),              
              ("C(=O)O","C(=O)[O-]",False),
              ("C(=O)[O-]","C(=O)[OH]",False),
              ("C(=O)[OH]","C(=O)[O-]",False),
              ("C(=O)[O-]","C(=O)OH",False),
              ("C(=O)OCC","C(=O)[O-]",False),
              ]
    
  def _testformula(self, num):
    smile1, smile2, result = self.formulas[num]
    m1 = smiles.text_to_mol( smile1)
    m2 = smiles.text_to_mol( smile2)
    self.assert_(m1.contains_substructure(m2)==result)

# this creates individual test for substructures
for i in range( len( TestSubstructure.formulas)):
  setattr( TestSubstructure, "testformula"+str(i+1), create_test(i,"_testformula"))

## // substructure testing


## SMILES equality testing

class TestEqualSMILES(unittest.TestCase):

  formulas = [("Sc1ccccc1","S-c1ccccc1",True),  # check Sc (in PT scandium) bug
              ("Oc1ccccc1","O-c1ccccc1",True),
              ("c1ccccc1","C:1:C:C:C:C:C:1", True),
              ("c1ccccc1","[CH]:1:[CH]:[CH]:[CH]:[CH]:[CH]:1", True),
              ("c1cscc1","C1=C-S-C=C1", True),
              ("c1ccccc1","C1=CC=CC=C1", True),
              ("c1ccccc1","C=1C=CC=CC=1", True),
              ("C=1CC=1","C1C=C1", True),
              ("C12CC%01CC2","C12CC1CC2", True),
              ("C%01%02CC%01CC2","C12CC1CC2", True),
              ("Oc%11ccccc%11","Oc1ccccc1",True),
              ("C%99CCCCC%992CCCCC2","C1CCCCC12CCCCC2",True),
              ("C%99CCCCC%88%99CCCCC%88","C1CCCCC12CCCCC2",True),
              ("H","[H]", False),
              ("C","[2H]C", False),
              ("[C]","[CH0]", True),
              ]
    
  def _testformula(self, num):
    smile1, smile2, result = self.formulas[num]
    m1 = smiles.text_to_mol( smile1)
    m2 = smiles.text_to_mol( smile2)
    self.assert_(molecule.equals(m1,m2,level=3)==result)

# this creates individual test for substructures
for i in range( len( TestEqualSMILES.formulas)):
  setattr( TestEqualSMILES, "testformula"+str(i+1), create_test(i,"_testformula"))

## // SMILES equality testing


## SMILES reading testing

class TestSMILESReading(unittest.TestCase):

  formulas = [("Sc1ccccc1",("C6H6S",)),
              ("Oc1ccccc1",("C6H6O",)),
              ("Oc%11ccccc%11",("C6H6O",)),
              ("[Na+].[Cl-]", ("Na","Cl")),
              ("[O-]c1ccccc1.[Na+]",("C6H5O","Na")),
              ("O=C[O-].[NH4+]",("CHO2","H4N")),
              ("c1ccccc1-c1ccccc1",("C12H10",)),
              ("c1cscc1",("C4H4S",)),
              ]
    
  def _testformula(self, num):
    smile1, sum_forms = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    for i,mol in enumerate( mols):
      self.assert_( i < len( sum_forms))
      self.assertEqual( str( mol.get_formula_dict()), sum_forms[i]) 

  def test_empty_smiles( self):
    conv = smiles.converter()
    for text in ("", "  "):
      mols = conv.read_text( "")
      self.assertEqual( mols, [])


# this creates individual test for substructures
for i in range( len( TestSMILESReading.formulas)):
  setattr( TestSMILESReading, "testformula"+str(i+1), create_test(i,"_testformula"))

## // SMILES reading testing


## SMILES Reaction support

class TestSMILESReactionSupport(unittest.TestCase):
  
  def test1(self):
    """tests handling of reactions by the SMILES reader on a preparation of methyl-formate"""
    c = smiles.converter()
    reacts = c.read_text( "O=CO.CO>[H+]>O=COC.O")
    self.assertEqual( reacts, c.result)
    self.assertEqual( c.last_status, c.STATUS_OK)
    self.assertEqual( len( reacts), 1)
    react = reacts[0]
    self.assertEqual( len( react.reactants), 2)
    self.assertEqual( len( react.reactants[0].molecule.atoms), 3)
    self.assertEqual( len( react.reactants[1].molecule.atoms), 2)
    self.assertEqual( len( react.reagents[0].molecule.atoms), 1)
    self.assertEqual( len( react.products[0].molecule.atoms), 4)
    self.assertEqual( len( react.products[1].molecule.atoms), 1)
    self.assertEqual( str( react.products[0].molecule.get_formula_dict()), "C2H4O2")

  def test2(self):
    """test reactions with some empty parts"""
    c = smiles.converter()
    reacts = c.read_text( "C=C.[H][H]>>CC")
    self.assertEqual( len( reacts), 1)
    react = reacts[0]
    self.assertEqual( len( react.reagents), 0)
    self.assertEqual( len( react.reactants), 2)
    self.assertEqual( len( react.products), 1)




## // SMILES Reaction support

## Reaction test

from . import reaction

class TestReactionComponent(unittest.TestCase):

  def test1(self):
    mol = smiles.text_to_mol( "CCCO")
    rc = reaction.reaction_component( mol, 2)
    self.assertEqual( rc.stoichiometry, 2)
    self.assertRaises( Exception, reaction.reaction_component, mol, "x")
    self.assertRaises( Exception, rc._set_stoichiometry, "x")
    self.assertRaises( Exception, reaction.reaction_component, 2, 2)
    self.assertRaises( Exception, rc._set_molecule, "x")    


## Explicit hydrogens, occupied_valency and free_valency testing

class TestValency(unittest.TestCase):

  # (formula, (N explicit_hydrogens, N occupied_valency, N free_valency))

  formulas = [("CN",(0,1,2)),
              ("C[NH2]",(2,3,0)),
              ("C[NH]C",(1,3,0)),
              ("CN(C)C",(0,3,0)),
              ]
    
  def _testformula(self, num):
    smile1, (explicit_hs, occupied_v, free_v) = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    self.assertEqual( len( mols), 1)
    mol = mols[0]
    ns = [v for v in mol.vertices if v.symbol == 'N']
    self.assertEqual( len( ns), 1)
    n = ns[0]
    self.assertEqual( n.explicit_hydrogens, explicit_hs)
    self.assertEqual( n.free_valency, free_v)
    self.assertEqual( n.occupied_valency, occupied_v)

# this creates individual test
for i in range( len( TestValency.formulas)):
  setattr( TestValency, "testformula"+str(i+1), create_test(i,"_testformula"))

## // Explicit hydrogens, occupied_valency and free_valency testing


## Charge computation testing

class TestCharge(unittest.TestCase):

  formulas = [("",()),
              ("c1ccccc1",(0,)),
              ("[Na+].[Cl-]", (1,-1)),
              ("[O-]c1ccccc1.[Na+]",(-1,1)),
              ("O=C[O-].[NH4+]",(-1,1)),
              ("[O-]CC[O-].[K+].[K+]", (-2,1,1)),
              ("[OH-].[OH-].[Ca+2]",(-1,-1,2)),
              ]
    
  def _testformula(self, num):
    smile1, charges = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    comp_charges = [mol.charge for mol in mols]
    comp_charges.sort()
    charges = list( charges)
    charges.sort()
    self.assertEqual( comp_charges, charges)

# this creates individual test for substructures
for i in range( len( TestCharge.formulas)):
  setattr( TestCharge, "testformula"+str(i+1), create_test(i,"_testformula"))

## // Charge computation testing


## Stereochemistry testing

class TestStereo(unittest.TestCase):

  formulas = [("c1ccccc1",()),
              (r"C\C=C/C", (1,)),
              (r"C\C=C/C=C/C=C\C", (-1,1,1)),
              (r"C\C(\O)=C/C=C/C=C\C", (-1,-1,1,1)),
              (r"O\C=C=C=C/N=C/Br", (-1,1)),
              (r"O\C(\N)=C/C=C\C=C\Cl", (-1,-1,1,1))
              ]
    
  def _testformula(self, num):
    smile1, directions = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    self.assertEqual( len( mols), 1)
    mol = mols[0]
    sts = [st.value==st.SAME_SIDE and 1 or -1 for st in mol.stereochemistry]
    sts.sort()
    sts = tuple( sts)
    self.assertEqual( sts, directions)

# this creates individual test for substructures
for i in range( len( TestStereo.formulas)):
  setattr( TestStereo, "testformula"+str(i+1), create_test(i,"_testformula"))


class TestStereo2(unittest.TestCase):
  """tests if stereochemistry of structure read from smiles
  remains the same when recoded back to smiles and read again."""

  formulas = [r"C\C=C/C",
              r"N\C=C/C=C/Cl",
              r"O\C(\N)=C/C=C\C=C\Cl",
              r"O\C=C=C=C/N=C/Br",
              r"C/C(Cl)=C(\O)C",
              ]
    
  def _testformula(self, num):
    def create_st_sum( st):
      symbols = [a.symbol for a in (st.references[0],st.references[-1])]
      symbols.sort()
      return tuple( symbols + [st.value==st.SAME_SIDE and 1 or -1])
      
    # round 1
    smile1 = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    self.assertEqual( len( mols), 1)
    mol = mols[0]
    sts1 = [create_st_sum( st) for st in mol.stereochemistry]
    # round 2
    smile2 = conv.mols_to_text( [mol])
    mols = conv.read_text( smile1)
    self.assertEqual( len( mols), 1)
    mol = mols[0]
    sts2 = [create_st_sum( st) for st in mol.stereochemistry]
    for stsum in sts2:
      self.assertEqual( (stsum in sts1), True)

# this creates individual test for substructures
for i in range( len( TestStereo2.formulas)):
  setattr( TestStereo2, "testformula"+str(i+1), create_test(i,"_testformula"))


class TestStereo3(unittest.TestCase):
  """tests if stereochemistry of structure read from smiles
  remains the same when coordinates are calculates, stereo information
  thrown away and recalculated from the coords."""

  formulas = [r"C\C=C/C",
              r"N\C=C/C=C/Cl",
              r"O\C(\N)=C/C=C\C=C\Cl",
              r"O\C=C=C=C/N=C/Br",
              r"C/C(Cl)=C(\O)N",
              ]
    
  def _testformula(self, num):
    def create_st_sum( st):
      symbols = [a.symbol for a in (st.references[0],st.references[-1])]
      symbols.sort()
      return tuple( symbols + [st.value==st.SAME_SIDE and 1 or -1])
      
    # round 1
    smile1 = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    self.assertEqual( len( mols), 1)
    mol = mols[0]
    sts1 = [create_st_sum( st) for st in mol.stereochemistry]
    # round 2
    mol.stereochemistry = []
    mol.detect_stereochemistry_from_coords()
    sts2 = [create_st_sum( st) for st in mol.stereochemistry]
    for stsum in sts2:
      if not stsum in sts1:
        # not necessarily an error, we must count changes
        changes = 0
        for (v1,v2) in zip( stsum, sts1[0]):
          if v1 != v2:
            changes += 1
        self.assertEqual( changes % 2, 0)
      else:
        self.assertEqual( (stsum in sts1), True)

# this creates individual test for substructures
for i in range( len( TestStereo3.formulas)):
  setattr( TestStereo3, "testformula"+str(i+1), create_test(i,"_testformula"))


## // Stereo testing


## Graph matching algorithm testing

class TestGraphMatching(unittest.TestCase):
  """tests if graph algorithm to find maximum matchings works"""

  formulas = [("C1CCCCC1",3,0),  # smiles, matching pairs, exposed vertices
              ("C1CCCCC1C",3,1),
              ("C1CCCCC1CC",4,0),
              ("C1CC1C(C)C",2,2),
              ]
    
  def _testformula(self, num):
    smile1,match_pair_num,exposed_vs_num = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    self.assertEqual( len( mols), 1)
    mol = mols[0]
    mate, nrex = mol.get_maximum_matching()
    self.assertEqual( nrex,exposed_vs_num)
    self.assertEqual( len( [v for v,m in mate.items() if m!=0]), 2*match_pair_num)

# this creates individual test for substructures
for i in range( len( TestGraphMatching.formulas)):
  setattr( TestGraphMatching, "testformula"+str(i+1), create_test(i,"_testformula"))


## // Graph matching algorithm testing




if __name__ == '__main__':
  import sys
  if len( sys.argv) > 1 and sys.argv[1] == "-v":
    tests = [globals()[k] for k in dir() if type(globals()[k])==type and issubclass( globals()[k], unittest.TestCase)]
    ss = []
    for test in tests:
      s1 = unittest.makeSuite( test)
      unittest.TextTestRunner(verbosity=2).run(s1)
  else:
    unittest.main()
