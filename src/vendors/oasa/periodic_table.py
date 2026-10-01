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
#
#
#
#--------------------------------------------------------------------------


#from generators import lazy_map, take
import re
import operator
import types
from functools import reduce


"""periodic table as a dictionary, plus functions for molecular
formula manipulation and computation"""

# data come from kalzium and http://bodr.svn.sourceforge.net/viewvc/*checkout*/bodr/trunk/bodr/elements/elements.xml?revision=34&content-type=text%2Fplain with some manual changes

periodic_table = {
"H"  : {"ord":   1, "weight":   1.0079, "exact_mass":   1.00782503, "en": 2.10, "els":  1, "valency": (1,)},
"He" : {"ord":   2, "weight":   4.0026, "exact_mass":   4.00260325, "en": 0.00, "els":  8, "valency": (0, 2)},
"Li" : {"ord":   3, "weight":   6.9410, "exact_mass":   7.01600455, "en": 0.98, "els":  1, "valency": (1,)},
"Be" : {"ord":   4, "weight":   9.0122, "exact_mass":   9.01218220, "en": 1.57, "els":  2, "valency": (2,)},
"B"  : {"ord":   5, "weight":  10.8110, "exact_mass":  11.00930540, "en": 2.04, "els":  3, "valency": (3,)},
"C"  : {"ord":   6, "weight":  12.0107, "exact_mass":  12.00000000, "en": 2.55, "els":  4, "valency": (4, 2)},
"N"  : {"ord":   7, "weight":  14.0067, "exact_mass":  14.00307400, "en": 3.04, "els":  5, "valency": (3, 5)},
"O"  : {"ord":   8, "weight":  15.9994, "exact_mass":  15.99491462, "en": 3.44, "els":  6, "valency": (2,)},
"F"  : {"ord":   9, "weight":  18.9984, "exact_mass":  18.99840322, "en": 3.98, "els":  7, "valency": (1,)},
"Ne" : {"ord":  10, "weight":  20.1797, "exact_mass":  19.99244018, "en": 0.00, "els":  8, "valency": (0, 2)},
"Na" : {"ord":  11, "weight":  22.9898, "exact_mass":  22.98976928, "en": 0.93, "els":  1, "valency": (1,)},
"Mg" : {"ord":  12, "weight":  24.3050, "exact_mass":  23.98504170, "en": 1.31, "els":  2, "valency": (2,)},
"Al" : {"ord":  13, "weight":  26.9815, "exact_mass":  26.98153863, "en": 1.61, "els":  3, "valency": (3,)},
"Si" : {"ord":  14, "weight":  28.0855, "exact_mass":  27.97692653, "en": 1.91, "els":  4, "valency": (4,)},
"P"  : {"ord":  15, "weight":  30.9738, "exact_mass":  30.97376163, "en": 2.19, "els":  5, "valency": (3, 5)},
"S"  : {"ord":  16, "weight":  32.0650, "exact_mass":  31.97207100, "en": 2.58, "els":  6, "valency": (2, 4, 6)},
"Cl" : {"ord":  17, "weight":  35.4530, "exact_mass":  34.96885268, "en": 3.16, "els":  7, "valency": (1, 3, 5, 7)},
"Ar" : {"ord":  18, "weight":  39.9480, "exact_mass":  39.96238312, "en": 0.00, "els":  8, "valency": (0, 2)},
"K"  : {"ord":  19, "weight":  39.0983, "exact_mass":  38.96370668, "en": 0.82, "els":  1, "valency": (1,)},
"Ca" : {"ord":  20, "weight":  40.0780, "exact_mass":  39.96259098, "en": 1.00, "els":  2, "valency": (2,)},
"Sc" : {"ord":  21, "weight":  44.9559, "exact_mass":  44.95591190, "en": 1.36, "els":  3, "valency": (3, 1)},
"Ti" : {"ord":  22, "weight":  47.8670, "exact_mass":  47.94794630, "en": 1.54, "els":  4, "valency": (4, 3)},
"V"  : {"ord":  23, "weight":  50.9415, "exact_mass":  50.94395950, "en": 1.63, "els":  5, "valency": (2, 4, 5)},
"Cr" : {"ord":  24, "weight":  51.9961, "exact_mass":  51.94050750, "en": 1.66, "els":  6, "valency": (2, 3, 6)},
"Mn" : {"ord":  25, "weight":  54.9380, "exact_mass":  54.93804510, "en": 1.55, "els":  7, "valency": (2, 3, 4, 6, 7)},
"Fe" : {"ord":  26, "weight":  55.8450, "exact_mass":  55.93493750, "en": 1.83, "els":  8, "valency": (0, 2, 3)},
"Co" : {"ord":  27, "weight":  58.9332, "exact_mass":  58.93319500, "en": 1.88, "els":  8, "valency": (2, 3)},
"Ni" : {"ord":  28, "weight":  58.6934, "exact_mass":  57.93534290, "en": 1.91, "els":  8, "valency": (2, 3)},
"Cu" : {"ord":  29, "weight":  63.5460, "exact_mass":  62.92959750, "en": 1.90, "els":  1, "valency": (0, 1, 2)},
"Zn" : {"ord":  30, "weight":  65.3900, "exact_mass":  63.92914220, "en": 1.65, "els":  2, "valency": (2,)},
"Ga" : {"ord":  31, "weight":  69.7230, "exact_mass":  68.92557360, "en": 1.81, "els":  3, "valency": (3,)},
"Ge" : {"ord":  32, "weight":  72.6400, "exact_mass":  73.92117780, "en": 2.01, "els":  4, "valency": (4,)},
"As" : {"ord":  33, "weight":  74.9216, "exact_mass":  74.92159650, "en": 2.18, "els":  5, "valency": (3, 5)},
"Se" : {"ord":  34, "weight":  78.9600, "exact_mass":  79.91652130, "en": 2.55, "els":  6, "valency": (2, 4, 6)},
"Br" : {"ord":  35, "weight":  79.9040, "exact_mass":  78.91833710, "en": 2.96, "els":  7, "valency": (1, 3, 5)},
"Kr" : {"ord":  36, "weight":  83.8000, "exact_mass":  83.91150700, "en": 0.00, "els":  8, "valency": (0, 2)},
"Rb" : {"ord":  37, "weight":  85.4678, "exact_mass":  84.91178974, "en": 0.82, "els":  1, "valency": (1,)},
"Sr" : {"ord":  38, "weight":  87.6200, "exact_mass":  87.90561210, "en": 0.95, "els":  2, "valency": (2,)},
"Y"  : {"ord":  39, "weight":  88.9059, "exact_mass":  88.90584830, "en": 1.22, "els":  3, "valency": (3,)},
"Zr" : {"ord":  40, "weight":  91.2240, "exact_mass":  89.90470440, "en": 1.33, "els":  4, "valency": (4,)},
"Nb" : {"ord":  41, "weight":  92.9064, "exact_mass":  92.90637810, "en": 1.60, "els":  5, "valency": (3, 5)},
"Mo" : {"ord":  42, "weight":  95.9400, "exact_mass":  97.90540820, "en": 2.16, "els":  6, "valency": (3, 5, 6)},
"Tc" : {"ord":  43, "weight":  98.9063, "exact_mass":  97.90721600, "en": 1.90, "els":  7, "valency": (5, 7)},
"Ru" : {"ord":  44, "weight": 101.0700, "exact_mass": 101.90434930, "en": 2.20, "els":  8, "valency": (3, 4, 6, 8)},
"Rh" : {"ord":  45, "weight": 102.9055, "exact_mass": 102.90550400, "en": 2.28, "els":  8, "valency": (3, 4)},
"Pd" : {"ord":  46, "weight": 106.4200, "exact_mass": 105.90348600, "en": 2.20, "els":  8, "valency": (2, 4)},
"Ag" : {"ord":  47, "weight": 107.8682, "exact_mass": 106.90509700, "en": 1.93, "els":  1, "valency": (1,)},
"Cd" : {"ord":  48, "weight": 112.4110, "exact_mass": 113.90335850, "en": 1.69, "els":  2, "valency": (2,)},
"In" : {"ord":  49, "weight": 114.8180, "exact_mass": 114.90387800, "en": 1.78, "els":  3, "valency": (3,)},
"Sn" : {"ord":  50, "weight": 118.7100, "exact_mass": 119.90219470, "en": 1.96, "els":  4, "valency": (2, 4)},
"Sb" : {"ord":  51, "weight": 121.7600, "exact_mass": 120.90381570, "en": 2.05, "els":  5, "valency": (3, 5)},
"Te" : {"ord":  52, "weight": 127.6000, "exact_mass": 129.90622440, "en": 2.10, "els":  6, "valency": (2, 4, 6)},
"I"  : {"ord":  53, "weight": 126.9045, "exact_mass": 126.90447300, "en": 2.66, "els":  7, "valency": (1, 3, 5, 7)},
"Xe" : {"ord":  54, "weight": 131.2930, "exact_mass": 131.90415350, "en": 2.60, "els":  8, "valency": (0, 2)},
"Cs" : {"ord":  55, "weight": 132.9055, "exact_mass": 132.90545190, "en": 0.79, "els":  1, "valency": (1,)},
"Ba" : {"ord":  56, "weight": 137.2370, "exact_mass": 137.90524720, "en": 0.89, "els":  2, "valency": (2,)},
"La" : {"ord":  57, "weight": 138.9055, "exact_mass": 138.90635330, "en": 1.10, "els":  3, "valency": (3,)},
"Ce" : {"ord":  58, "weight": 140.1160, "exact_mass": 139.90543870, "en": 1.12, "els":  4, "valency": (3, 4)},
# conflicts with propyl
"Pr" : {"ord":  59, "weight": 140.9076, "exact_mass": 140.90765280, "en": 1.13, "els":  5, "valency": (3, 4)},
"Nd" : {"ord":  60, "weight": 144.2420, "exact_mass": 141.90772330, "en": 1.14, "els":  6, "valency": (3,)},
"Pm" : {"ord":  61, "weight": 145.0000, "exact_mass": 144.91274900, "en": 1.13, "els":  7, "valency": (3,)},
"Sm" : {"ord":  62, "weight": 150.3600, "exact_mass": 151.91973240, "en": 1.17, "els":  8, "valency": (2, 3)},
"Eu" : {"ord":  63, "weight": 151.9640, "exact_mass": 152.92123030, "en": 1.20, "els":  8, "valency": (2, 3)},
"Gd" : {"ord":  64, "weight": 157.2500, "exact_mass": 157.92410390, "en": 1.20, "els":  8, "valency": (3,)},
"Tb" : {"ord":  65, "weight": 158.9254, "exact_mass": 158.92534680, "en": 1.10, "els":  1, "valency": (3, 4)},
"Dy" : {"ord":  66, "weight": 162.5000, "exact_mass": 163.92917480, "en": 1.22, "els":  2, "valency": (3,)},
"Ho" : {"ord":  67, "weight": 164.9303, "exact_mass": 164.93032210, "en": 1.23, "els":  3, "valency": (3,)},
"Er" : {"ord":  68, "weight": 167.2590, "exact_mass": 165.93029310, "en": 1.24, "els":  4, "valency": (3,)},
"Tm" : {"ord":  69, "weight": 168.9342, "exact_mass": 168.93421330, "en": 1.25, "els":  5, "valency": (2, 3)},
"Yb" : {"ord":  70, "weight": 173.0400, "exact_mass": 173.93886210, "en": 1.10, "els":  6, "valency": (2, 3)},
"Lu" : {"ord":  71, "weight": 174.9670, "exact_mass": 174.94077180, "en": 1.27, "els":  7, "valency": (3,)},
"Hf" : {"ord":  72, "weight": 178.4900, "exact_mass": 179.94655000, "en": 1.30, "els":  4, "valency": (4,)},
"Ta" : {"ord":  73, "weight": 180.9479, "exact_mass": 180.94799580, "en": 1.50, "els":  5, "valency": (5,)},
"W"  : {"ord":  74, "weight": 183.8400, "exact_mass": 183.95093120, "en": 2.36, "els":  6, "valency": (6,)},
"Re" : {"ord":  75, "weight": 186.2070, "exact_mass": 186.95575310, "en": 1.90, "els":  7, "valency": (7,)},
"Os" : {"ord":  76, "weight": 190.2300, "exact_mass": 191.96148070, "en": 2.20, "els":  8, "valency": (4, 6, 8)},
"Ir" : {"ord":  77, "weight": 192.2170, "exact_mass": 192.96292640, "en": 2.20, "els":  8, "valency": (3, 4, 6)},
"Pt" : {"ord":  78, "weight": 195.0780, "exact_mass": 194.96479110, "en": 2.28, "els":  8, "valency": (2, 4)},
"Au" : {"ord":  79, "weight": 196.9665, "exact_mass": 196.96656870, "en": 2.54, "els":  1, "valency": (1, 3)},
"Hg" : {"ord":  80, "weight": 200.5900, "exact_mass": 201.97064300, "en": 2.00, "els":  2, "valency": (1, 2)},
"Tl" : {"ord":  81, "weight": 204.3833, "exact_mass": 204.97442750, "en": 2.04, "els":  3, "valency": (1, 3)},
"Pb" : {"ord":  82, "weight": 207.2000, "exact_mass": 207.97665210, "en": 2.33, "els":  4, "valency": (2, 4)},
"Bi" : {"ord":  83, "weight": 208.9804, "exact_mass": 208.98039870, "en": 2.02, "els":  5, "valency": (3, 5)},
"Po" : {"ord":  84, "weight": 208.9824, "exact_mass": 208.98243040, "en": 2.00, "els":  6, "valency": (2, 4, 6)},
"At" : {"ord":  85, "weight": 209.9871, "exact_mass": 209.98714800, "en": 2.20, "els":  7, "valency": (1, 7)},
"Rn" : {"ord":  86, "weight": 222.0176, "exact_mass": 222.01757770, "en": 0.00, "els":  8, "valency": (0, 2)},
"Fr" : {"ord":  87, "weight": 223.0000, "exact_mass": 223.01973590, "en": 0.70, "els":  1, "valency": (1,)},
"Ra" : {"ord":  88, "weight": 226.0254, "exact_mass": 226.02540980, "en": 0.90, "els":  2, "valency": (2,)},
# conflicts with acetyl
"Ac" : {"ord":  89, "weight": 227.0000, "exact_mass": 227.02775210, "en": 1.10, "els":  3, "valency": (2, 3)},
"Th" : {"ord":  90, "weight": 232.0381, "exact_mass": 232.03805530, "en": 1.30, "els":  4, "valency": (3, 4)},
"Pa" : {"ord":  91, "weight": 231.0359, "exact_mass": 231.03588400, "en": 1.50, "els":  5, "valency": (4, 5)},
"U"  : {"ord":  92, "weight": 238.0289, "exact_mass": 238.05078820, "en": 1.38, "els":  6, "valency": (3, 4, 5, 6)},
"Np" : {"ord":  93, "weight": 237.0000, "exact_mass": 237.04817340, "en": 1.36, "els":  7, "valency": (4, 3, 5, 7)},
"Pu" : {"ord":  94, "weight": 244.0000, "exact_mass": 244.06420400, "en": 1.28, "els":  8, "valency": (3, 4, 5, 6, 7)},
"Am" : {"ord":  95, "weight": 243.0000, "exact_mass": 243.06138110, "en": 1.30, "els":  8, "valency": (2, 3, 4, 5, 6)},
"Cm" : {"ord":  96, "weight": 247.0000, "exact_mass": 247.07035400, "en": 1.30, "els":  8, "valency": (3, 4)},
"Bk" : {"ord":  97, "weight": 247.0000, "exact_mass": 247.07030700, "en": 1.30, "els":  1, "valency": (3, 4)},
"Cf" : {"ord":  98, "weight": 251.0000, "exact_mass": 251.07958700, "en": 1.30, "els":  2, "valency": (3, 4)},
"Es" : {"ord":  99, "weight": 252.0000, "exact_mass": 252.08298000, "en": 1.30, "els":  3, "valency": (2, 3)},
"Fm" : {"ord": 100, "weight": 257.0000, "exact_mass": 257.09510500, "en": 1.30, "els":  4, "valency": (2, 3)},
"Md" : {"ord": 101, "weight": 258.0000, "exact_mass": 258.09843100, "en": 1.30, "els":  5, "valency": (2, 3)},
"No" : {"ord": 102, "weight": 259.0000, "exact_mass": 259.10103000, "en": 1.30, "els":  6, "valency": (2, 3)},
"Lr" : {"ord": 103, "weight": 262.0000, "exact_mass": 262.10963000, "en": 1.30, "els":  7, "valency": (3,)},
"Rf" : {"ord": 104, "weight": 261.0000, "exact_mass": 261.10877000, "en": 0.00, "els":  4, "valency": (3,)},
"Db" : {"ord": 105, "weight": 262.0000, "exact_mass": 262.11408000, "en": 0.00, "els":  5, "valency": (0,)},
"Sg" : {"ord": 106, "weight": 266.0000, "exact_mass": 263.11832000, "en": 0.00, "els":  6, "valency": (0,)},
"Bh" : {"ord": 107, "weight": 264.0000, "exact_mass": 264.12460000, "en": 0.00, "els":  7, "valency": (0,)},
"Hs" : {"ord": 108, "weight": 277.0000, "exact_mass": 265.13009000, "en": 0.00, "els":  8, "valency": (0,)},
"Mt" : {"ord": 109, "weight": 268.0000, "exact_mass": 268.13873000, "en": 0.00, "els":  8, "valency": (0,)},
"Ds" : {"ord": 110, "weight": 281.0000, "exact_mass": 271.14606000, "en": 0.00, "els":  8, "valency": (0,)},
"X": {'query': True, 'ord': 300, 'key': '"X"', 'weight': 0, 'valency': (1,)},          # halogen
"Q": {'query': True, 'ord': 301, 'key': '"Q"', 'weight': 0, 'valency': (1, 2, 3, 4)},  # anything not H or C
"A": {'query': True, 'ord': 302, 'key': '"A"', 'weight': 0, 'valency': (1, 2, 3, 4)},  # anything not H
"R": {'query': True, 'ord': 303, 'key': '"R"', 'weight': 0, 'valency': (1, 2, 3, 4)},  # anything
}

# elements that accept cations and raise their valency; for each element a valency is specified
# because this property is valency (oxidation state) specific 
accept_cation = {'N': 3, 'P': 3, 'O': 2, 'S': 2, 'Se': 2}

# elements that accept anions and raise their valency; for each element a valency is specified
# because this property is valency (oxidation state) specific 
accept_anion = {'B': 3, 'Al': 3, 'P': 5}


class composition_dict( dict):
  """special dict that automatically converts itself to human readable composition on str()"""
  def __str__( self):
    ret = ''
    for n in ('C','H'):
      if n in self:
        if ret:
          ret += ', '
        ret += "%s: %2.3f%%" % (n, self[n])
    k = self.keys()
    k.sort()
    for n in self:
      if n not in ('C','H'):
        if ret:
          ret += ', '
        ret += "%s: %2.3f%%" % (n, self[n])
    return ret

class formula_dict( dict):
  """special dict that automatically converts itself to human readable
  formula on str(). Implements += for convenient formula concatenation"""

  def __init__( self, form=None):
    dict.__init__( self)
    ## incomplete means that there were some problems to fully convert a formula to this dict
    self.incomplete = 0
    if type( form) in (str, str):
      self.read_formula_string( form)
    elif type( form) == dict:
      for key, val in form.items():
        if key in periodic_table and type( val) == int:
          self[ key] = val
        else:
          raise ValueError("some of the dictionary entries are not valid for formula_dict (%s => %s)" % (str( key), str( val)))
  
  def __str__( self, reverse=0):
    sum = ''
    k = self.sorted_keys()
    if reverse:
      k.reverse()
    for s in k:
      if self[s] == 1:
        num = ''
      else:
        num = str( self[s])
      sum += s+num
    return sum

  def __iadd__( self, other):
    for s in other:
      if s in self:
        self[s] += other[s]
      else:
        self[s] = other[s]
    return self

  def __add__( self, other):
    ret = formula_dict()
    for form in (self, other):
      for s in form:
        if s in ret:
          ret[s] += form[s]
        else:
          ret[s] = form[s]
    return ret

  def __mul__( self, other):
    if not type( other) == int:
      raise TypeError("formula_dict can be only multiplied by an integer")
    res = formula_dict()
    for key in self.keys():
      res[key] = other * self[key]
    return res


  def get_element_fraction( self, element):
    if element in self:
      return self[element]*periodic_table[element]['weight']/self.get_molecular_weight()
    return 0

  def get_molecular_weight( self):
    tot = 0
    for i in self:
      tot += self[i]* periodic_table[i]['weight']
    return tot

  def get_exact_molecular_mass( self):
    tot = 0
    for i in self:
      tot += self[i]* periodic_table[i]['exact_mass']
    return tot

  
  def keys_in_order( self):
    return self.sorted_keys()

  def sorted_keys( self):
    k = self.keys()
    ret = []
    if 'C' in k:
      for a in ('C','H'):
        if a in k:
          ret.append( a)
          k.remove( a)
      k.sort()
      return ret+k
    else:
      k.sort()
      return k
    

  def read_formula_string( self, form):
    is_formula = re.compile("^([A-Z][a-z]?[0-9]*)*$")
    #form = "".join( form.split("."))
    form = form.replace( ".", "")
    if not is_formula.match( form):
      return None
    chunks = re.split( "([A-Z][a-z]*)", form)
    del chunks[0]
    for i in range( 0, len( chunks), 2):
      if chunks[i] in self:
        if chunks[i+1] == '':
          j = 1
        else:
          j = int( chunks[i+1])
        self[ chunks[i]] += j
      elif chunks[i] in periodic_table:
        if chunks[i+1] == '':
          j = 1
        else:
          j = int( chunks[i+1])
        self[ chunks[i]] = j
      else:
        self.incomplete = 1

  def get_html_repr_as_string( self, outer_element=None, reverse=0):
    sum = ''
    k = self.sorted_keys()
    if reverse:
      k.reverse()
    for s in k:
      if self[s] == 1:
        num = ''
      else:
        num = '<sub>%d</sub>' % self[s]
      sum += s+num
    if outer_element:
      return '<%s>%s</%s>' % (outer_element, sum, outer_element)
    return sum
    
  def is_saturated_alkyl_chain( self):
    if (self.sorted_keys() == ['C','H']) and (self['H'] == 2*self['C']+1):
      return 1
    else:
      return 0

  def to_tuple( self):
    return tuple( reduce( operator.add, [[k,self[k]] for k in self.sorted_keys()], []))


    
def dict_to_composition( form):
  w = form.get_molecular_weight()
  ret = composition_dict()
  for s in form:
    ret[ s] = form.get_element_fraction(s) * 100
  return ret

def formula_to_weight( formula):
  return formula_dict( formula).get_molecular_weight()

def formula_to_formula( formula):
  return str( formula_dict( formula))

def formula_to_composition( formula):
  return dict_to_composition( formula_to_dict( formula))


## other support functions

def text_to_hydrogenated_atom( text):
  a = re.match( r'^([a-z]{1,2})(h)(\d*)$', text.lower())
  if a:
    atom = a.group( 1)
    hydrogens = a.group( 3)
  else:
    a = re.match( r'^(h)(\d*)([a-z]{1,2})$', text.lower())
    if a:
      atom = a.group( 3)
      hydrogens = a.group( 2)
    else:
      return None

  if atom.capitalize() in periodic_table:
    ret = formula_dict()
    ret[ atom.capitalize()] = 1
    if hydrogens:
      ret[ 'H'] = int( hydrogens)
    else:
      ret[ 'H'] = 1
    return ret
  else:
    return None



def gen_bit_masks( length):
  ret = length * [0]
  yield ret
  for i in range( 2 ** length):
    ret[0] += 1
    for j in range( length):
      if ret[j] == 2:
        ret[j] = 0
        if j == length-1:
          raise StopIteration
        else:
          ret[j+1] += 1
      else:
        break
    yield ret
  

