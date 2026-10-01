#--------------------------------------------------------------------------
#     This file is part of BKChem - a chemical drawing program
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

"""Module containing miscelanous functions used in BKChem that don't
fit anywhere else. Does not contain any objects."""

import math
import time
import string
from warnings import warn
import re
import types
import base64


def intersection( a, b):
  "returns intersection of 2 lists"
  ret = []
  for i in a:
    if i in b:
      ret.append( i)
  return ret


def difference( a,b):
  "returns difference of 2 lists ( a-b)"
  ret = list( a)  # needed for type conversion of tuple for instance
  for i in b:
    if i in ret:
      ret.remove( i)
  return ret


def signum( a):
  if a == 0:
    return 0
  elif a < 0:
    return -1
  else:
    return 1


def filter_unique( items):
  #t1 = time.time()
  ret = []
  for item in items:
    if item not in ret:
      ret.append( item)
  #print time.time() - t1
  return ret


#def sub_all_numbers( text):
#  return re.sub( "\d+", '<sub>\g<0><\\sub>', text)

def reverse_molecule_formula( formula):
  pass
  

def normalize_coords( coords):
  x1, y1, x2, y2 = coords
  if x2 < x1:
    x2, x1 = x1, x2
  if y2 < y1:
    y2, y1 = y1, y2
  return (x1, y1, x2, y2)


def list_difference( list):
  """return a list of differences between list members,
  the list is by 1 shorter than the original"""
  ret = []
  for i in range( len( list)-1):
    ret.append( list[i]-list[i+1])
  return ret


def split_number_and_unit( txt):
  try:
    v = float( txt)
    return v, ''
  except:
    pass
  cutter = re.compile( r"([+-]?\d*\.?\d*)\s*([a-zA-Z]*)")
  if txt:
    a = cutter.match( txt)
    if a and a.group(1):
      return float( a.group(1)), a.group(2)
  return None, None


def lazy_apply( function, arguments):
  """similar to apply but returns a callable (lambda) that performs the apply when called."""
  return lambda: apply( function, arguments)




def extend_bbox( bbox, pixels=1):
  minx = min( (bbox[0], bbox[2]))
  maxx = max( (bbox[0], bbox[2]))
  miny = min( (bbox[1], bbox[3]))
  maxy = max( (bbox[1], bbox[3]))
  return minx-pixels, miny-pixels, maxx+pixels, maxy+pixels


def smallest_common_bbox( bboxes):
  _x0, _y0, _x1, _y1 = None, None, None, None
  for (x0, y0, x1, y1) in bboxes:
    minx = min( x0, x1)
    maxx = max( x0, x1)
    miny = min( y0, y1)
    maxy = max( y0, y1)
    if not _x0 or minx < _x0:
      _x0 = minx
    if not _x1 or maxx > _x1:
      _x1 = maxx
    if not _y0 or miny < _y0:
      _y0 = miny
    if not _y1 or maxy > _y1:
      _y1 = maxy
  return _x1, _y1, _x0, _y0

      
    



def has_one_value_only( iterable):
  if not iterable:
    return 0
  a = iterable[0]
  for i in iterable:
    if a != i:
      return 0
  return 1



def plural_or_singular( iterable):
  """useful for string construction such as 'you have %d apple%s' % (len(apples), plural_or_singular( apples)"""
  if len( iterable) == 1:
    return ''
  else:
    return 's'


def len_and_ending( iterable):
  return (len( iterable), plural_or_singular( iterable))


def reverse( iterable):
  for i in range( len( iterable)-1, -1, -1):
    yield iterable[i]



# some helper, higher order functions

map_functions = lambda funcs, value: zip( apply, funcs, len(funcs)*[value])

something_true = lambda vals: len( filter( None, vals))

some_apply = lambda func, vals: something_true( map( func, vals))



# some results

isinstance_of_one = lambda obj, parents: some_apply( lambda x: isinstance( obj, x), parents)


  

def gen_variations(items, n):
  if n==0:
    yield []
  else:
    for i in range( len(items)-n+1):
      for v in gen_variations(items[i+1:],n-1):
        yield [items[i]]+v


def gen_combinations( items, n):
  if n==0:
    yield []
  else:
    for i in range( len( items)):
      for v in gen_combinations( [x for x in items if x!=items[i]], n-1):
        yield [items[i]]+v


def x_over_y( x, y):
  return factorial( x) / factorial( y) / factorial( x-y)


def factorial( x):
  ret = 1
  for i in range( 2, x+1):
    ret *= i
  return ret


def gen_variations_and_one( items, length):
  for i in items:
    vars = gen_variations( [x for x in items if x!=i], length-1)
    for vs in vars:
      yield vs + [i]


#print x_over_y( 11, 7)
#print len( list( gen_variations( [1,2,3,4,5,6,7,8,9,10,11,12], 8)))
