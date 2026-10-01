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

import xml.dom.minidom as dom
from . import dom_extensions
from . import transform
from . import geometry
import math
from . import misc
import operator
import copy
from functools import reduce


class svg_out:

  show_hydrogens_on_hetero = False
  margin = 15
  line_width = 2
  bond_width = 6
  # should individual parts of an edge be grouped together
  group_items = True


  def __init__( self):
    pass

  def mol_to_svg( self, mol, before=None, after=None):
    """before and after should be methods or functions that will take one
    argument - svg_out instance and do whatever it wants with it - usually
    adding something to the resulting DOM tree"""
    self.document = dom.Document()
    top = dom_extensions.elementUnder( self.document,
                                       "svg",
                                       attributes=(("xmlns", "http://www.w3.org/2000/svg"),
                                                   ("version", "1.0")))
    self.top = dom_extensions.elementUnder( top, "g",
                                            attributes=(("stroke", "#000"),
                                                        ("stroke-width", "1.0")))

    x1, y1, x2, y2 = None, None, None, None
    for v in mol.vertices:
      if x1 == None or x1 > v.x:
        x1 = v.x
      if x2 == None or x2 < v.x:
        x2 = v.x
      if y1 == None or y1 > v.y:
        y1 = v.y
      if y2 == None or y2 < v.y:
        y2 = v.y

    w = int( x2 - x1 + 2*self.margin)
    h = int( y2 - y1 + 2*self.margin)
    top.setAttribute( "width", str( w))
    top.setAttribute( "height", str( h))

    self.transformer = transform.transform()
    self.transformer.set_move( -x1+self.margin, -y1+self.margin)

    self.molecule = mol

    if before:
      before( self)

    for e in copy.copy( mol.edges):
      self._draw_edge( e)
    for v in mol.vertices:
      self._draw_vertex( v)

    if after:
      after( self)

    return self.document


  def paper_to_canvas_coord( self, x):
    return x


  def prepare_dumb_transformer( self):
    tr = transform.transform()
    tr.set_scaling( self.paper_to_canvas_coord( 1))
    return tr


  def _create_parent( self, item, top):
    if self.group_items:
      parent = dom_extensions.elementUnder( top, "g")
      if 'svg_id' in item.properties_:
        parent.setAttribute( "id", item.properties_['svg_id'])
    else:
      parent = self.top
    return parent
      

  def _draw_edge( self, e):
    v1, v2 = e.vertices
    start = self.transformer.transform_xy( v1.x, v1.y)
    end = self.transformer.transform_xy( v2.x, v2.y)
    parent = self._create_parent( e, self.top)

    if e.order == 1:
      self._draw_line( parent, start, end, line_width=self.line_width)

    if e.order == 2:
      side = 0
      # find how to center the bonds
      # rings have higher priority in setting the positioning
      for ring in self.molecule.get_smallest_independent_cycles():
        if v1 in ring and v2 in ring:
          side += reduce( operator.add, [geometry.on_which_side_is_point( start+end, (self.transformer.transform_xy( a.x,a.y))) for a in ring if a!=v1 and a!=v2])
      # if rings did not decide, use the other neigbors
      if not side:
        for v in v1.neighbors + v2.neighbors:
          if v != v1 and v!= v2:
            side += geometry.on_which_side_is_point( start+end, (self.transformer.transform_xy( v.x, v.y)))
      if side:
        self._draw_line( parent, start, end, line_width=self.line_width)
        x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], self.bond_width*misc.signum( side))
        self._draw_line( parent, (x1, y1), (x2, y2), line_width=self.line_width)
      else:
        for i in (1,-1):
          x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], i*self.bond_width*0.5)
          self._draw_line( parent, (x1, y1), (x2, y2), line_width=self.line_width)
        
        
    elif e.order == 3:
      self._draw_line( parent, start, end, line_width=self.line_width)
      for i in (1,-1):
        x1, y1, x2, y2 = geometry.find_parallel( start[0], start[1], end[0], end[1], i*self.bond_width*0.7)
        self._draw_line( parent, (x1, y1), (x2, y2), line_width=self.line_width)
    

  def _draw_vertex( self, v):
    parent = self._create_parent( v, self.top)

    if v.symbol != "C" or v.charge != 0 or v.multiplicity != 1:
      x = v.x - 5
      y = v.y + 6
      x1 = x
      x2 = x + 12
      y1 = y - 12
      y2 = y + 2
      text = v.symbol
      if self.show_hydrogens_on_hetero:
        if v.free_valency == 1:
          text += "H"
        elif v.free_valency > 1:
          text += "H%d" % v.free_valency
      # charge
      if v.charge == 1:
        text += "+"
      elif v.charge == -1:
        text += "-"
      elif v.charge > 1:
        text += str( v.charge) + "+"
      elif v.charge < -1:
        text += str( v.charge)
      # radicals
      if v.multiplicity in (2,3):
        self._draw_circle( parent, self.transformer.transform_xy((x2+x1)/2,y-17), fill_color="#000", opacity=1, radius=3)
        if v.multiplicity == 3:
          self._draw_circle( parent, self.transformer.transform_xy((x2+x1)/2,y+5), fill_color="#000", opacity=1, radius=3)  
      self._draw_rectangle( parent, self.transformer.transform_4( (x1, y1, x2, y2)), fill_color="#fff")
      self._draw_text( parent, self.transformer.transform_xy(x,y), text)


  def _draw_line( self, parent, start, end, line_width=1, capstyle=""):
    x1, y1 = start
    x2, y2 = end
    line = dom_extensions.elementUnder( parent, 'line',
                                        (( 'x1', str( x1)),
                                         ( 'y1', str( y1)),
                                         ( 'x2', str( x2)),
                                         ( 'y2', str( y2))))


  def _draw_text( self, parent, xy, text, font_name="Arial", font_size=16):
    x, y = xy
    dom_extensions.textOnlyElementUnder( parent, "text", text,
                                         (( "x", str( x)),
                                          ( "y", str( y)),
                                          ( "font-family", font_name),
                                          ( "font-size", str( font_size)),
                                          ( 'fill', "#000")))


  def _draw_rectangle( self, parent, coords, fill_color="#fff", stroke_color="#fff"):
    x, y, x2, y2 = coords
    dom_extensions.elementUnder( parent, 'rect',
                                 (( 'x', str( x)),
                                  ( 'y', str( y)),
                                  ( 'width', str( x2-x)),
                                  ( 'height', str( y2-y)),
                                  ( 'fill', fill_color),
                                  ( 'stroke', stroke_color)))

  def _draw_circle( self, parent, xy, radius=5, fill_color="#fff", stroke_color="#fff", id="", opacity=0):
    x, y = xy
    el = dom_extensions.elementUnder( parent, 'ellipse',
                                      (( 'cx', str( x)),
                                       ( 'cy', str( y)),
                                       ( 'rx', str( radius)),
                                       ( 'ry', str( radius)),
                                       ( 'stroke-width', "1"),
                                       ( 'fill', fill_color),
                                       ( 'stroke', stroke_color),
                                       ( 'fill-opacity', str(opacity)),
                                       ( 'stroke-opacity', str(opacity)),
                                      ))
    if id:
      el.setAttribute( "id", id)


def mol_to_svg( mol, filename):
  c = svg_out()
  tree = c.mol_to_svg( mol)
  f = file( filename, "w")
  f.write( tree.toxml())
  f.close()


if __name__ == "__main__":

  #import inchi
  #mol = inchi.text_to_mol( "1/C7H6O2/c8-7(9)6-4-2-1-3-5-6/h1-5H,(H,8,9)", include_hydrogens=False, calc_coords=30)
  #import smiles
  #mol = smiles.text_to_mol( "CC[CH]", calc_coords=40)
  from . import molfile
  mol = molfile.file_to_mol( file( "/home/beda/bkchem/bkchem/untitled0.mol", "r"))
  mol_to_svg( mol, "output.svg")


