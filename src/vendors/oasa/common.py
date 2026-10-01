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




def gen_combinations_of_series( series):
  """series is a list of lists (tuples), the generator yields
  lists by combining each element of each list with each other"""
  counter = len( series) * [0]
  end = [i-1 for i in map( len, series)]
  counter[0] = -1 # dirty trick
  while counter != end:
    for i, e in enumerate( end):
      if counter[i] < e:
        counter[i] += 1
        for j in range( i):
          counter[j] = 0
        break
    yield [s[ counter[ j]] for j,s in enumerate( series)]


def is_uniquely_sorted( series, sorting_function=None):
  """you put a *sorted* series inside and get the information whether all the items
  are unique (there are no two same items)"""
  sf = sorting_function or cmp
  for i in range( len( series)-1):
    if sf( series[i], series[i+1]) == 0:
      return False
  return True



def least_common_item( series):
  d = {}
  for i in series:
    d[i] = d.get( i, 0) + 1
  return d.keys()[ d.values().index( min( d.values()))]


    
