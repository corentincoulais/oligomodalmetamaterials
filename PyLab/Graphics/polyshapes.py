# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 18:11:14 2014

@author: coulais
"""

"""

Kirby Urner

4D Solutions

First published: Apr 29 2007

 

Suitable for spatial geometry and/or synergetics students.

 

Update May 10:  I'd forgotten 8 of the 20 Icosahedron triangles!  Added.

Update May 13:  Added Octahedron, Mite, Coupler

 

"""

 

from stickworks import Vector, Edge

from visual import color

from math import sqrt

 

phi = (sqrt(5) + 1)/2.0

        

class Polyhedron (object):

 

    # defaults may be overridden

    showfaces      = True

    showedges      = True

    showvertices   = True

    # default POV-Ray textures

    face_texture   = 'T_Stone14'   # from stones.inc

    vertex_texture = 'T_Silver_1A' # from metals.inc

    edge_texture   = 'T_Copper_4A' # from metals.inc

    

    def scale(self, scalefactor):

        newverts = {}

        for v in self.vertices:

            newverts[v] = self.vertices[v] * scalefactor

        return self.__class__(newverts)

 

    __mul__ = __rmul__ = scale

 

    def translate(self, vector):

        newverts = {}

        for v in self.vertices:

            newverts[v] = self.vertices[v] + vector

        return self.__class__(newverts)

 

    __add__ = __radd__ = translate

    

    def _distill(self):

 

        edges = []

        unique = set()

        

        for f in self.faces:

            for pair in zip(f , f[1:] + (f[0],)):

                unique.add( tuple(sorted(pair)) )

 

        for edge in unique:

            edges.append( Edge(self.vertices[edge[0]],self.vertices[edge[1]]) )

 

        return edges            

 

    def draw(self):

        # VPython wireframe view, native to stickworks.py

        for e in self.edges:

            e.draw()

            

 

class Amodule (Polyhedron) :

    pass

 

class Bmodule (Polyhedron) :

    pass

 

class Mite (Polyhedron) :

 

    def __init__(self,

                 verts  = dict(j = Vector(( 0,  1, 0)),

                               o = Vector(( 0,  0, 0)),

                               r = Vector(( 1,  0, 1)),

                               s = Vector(( 1,  0,-1)))):

                 

        # 4 vertices

        self.vertices = verts

 

        # 4 faces

        self.faces = (('j','o','r'),('j','r','s'),('j','s','o'),('o','r','s'))

 

        self.edges = self._distill()

        

class Smite (Polyhedron) :

    pass

 

class Coupler (Polyhedron) :

    

    def __init__(self,

                 verts  = dict(j = Vector(( 0,  1, 0)),

                               l = Vector(( 0, -1, 0)),

                               q = Vector((-1,  0, 1)),

                               r = Vector(( 1,  0, 1)),

                               s = Vector(( 1,  0,-1)),

                               t = Vector((-1,  0,-1)))):

 

        # 6 vertices

        self.vertices = verts

 

        # 8 faces

        self.faces = (('j','q','r'),('j','r','s'),('j','s','t'),('j','t','q'),

                      ('l','q','r'),('l','r','s'),('l','s','t'),('l','t','q'))

 

        self.edges = self._distill()

 

class Tetrahedron (Polyhedron) :

 

    def __init__(self,

                 verts  = dict(a = Vector((-1, -1, 1)),

                               b = Vector((-1,  1, -1)),

                               c = Vector((1, 1, 1)),

                               d = Vector((1, -1, -1)))):

        """

        Imagine a cube centered at the origin and with

        a positive octant vertex at (1,1,1).  Inscribe

        a regular tetrahedron as six face diagonals therein.

        """

        # 4 vertices

        self.vertices = verts

 

        # 4 faces

        self.faces = (('a','b','c'),('a','c','d'),

                      ('a','d','b'),('b','d','c'))

 

        self.edges = self._distill()

 

class Cube (Polyhedron):

 

    def __init__(self, verts = dict( a = Vector((-1, -1, 1)),

                                     b = Vector((-1,  1, -1)),

                                     c = Vector((1, 1, 1)),

                                     d = Vector((1, -1, -1)),

                                     e = Vector((1,  1, -1)),

                                     f = Vector((1, -1,  1)),

                                     g = Vector((-1, -1, -1)),

                                     h = Vector((-1, 1, 1)))):

 

        # 8 vertices

        self.vertices = verts

 

        # 6 faces

        self.faces = (('a','f','c','h'),('h','c','e','b'),

                      ('b','e','d','g'),('g','d','f','a'),

                      ('c','f','d','e'),('a','h','b','g'))

 

        self.edges = self._distill()

 

class Octahedron (Polyhedron):

 

    def __init__(self, verts = dict( i = Vector(( 0, 0, 1)),

                                     j = Vector(( 0, 1, 0)),

                                     k = Vector(( 0, 0,-1)),

                                     l = Vector(( 0,-1, 0)),

                                     m = Vector(( 1, 0, 0)),

                                     n = Vector((-1, 0, 0)))):

 

        # 6 vertices

        self.vertices = verts

 

        # 8 faces

        self.faces = (('i','l','m'),('i','m','j'),('i','j','n'),('i','n','l'),                      

                      ('k','l','m'),('k','m','j'),('k','j','n'),('k','n','l'))

 

        self.edges = self._distill()

                                     

class Dodecahedron (Polyhedron):

    pass

 

class Icosahedron (Polyhedron):

 

    def __init__(self, verts = dict(

            # 12 vertices at the corners of 3 mutually

            # orthogonal golden rectangles

            xya=Vector(( phi/2, 0.5, 0.0)), # phi rectangle in xy

            xyb=Vector(( phi/2,-0.5, 0.0)),

            xyc=Vector((-phi/2,-0.5, 0.0)),

            xyd=Vector((-phi/2, 0.5, 0.0)),

            #-----------------------------

            xza=Vector((-0.5, 0.0, phi/2)), # Phi rectangle in xz 

            xzb=Vector(( 0.5, 0.0, phi/2)),

            xzc=Vector(( 0.5, 0.0,-phi/2)),

            xzd=Vector((-0.5, 0.0,-phi/2)),

            #-----------------------------

            yza=Vector(( 0.0, phi/2, 0.5)), # Phi rectangle in yz 

            yzb=Vector(( 0.0, phi/2,-0.5)),

            yzc=Vector(( 0.0,-phi/2,-0.5)),

            yzd=Vector(( 0.0,-phi/2, 0.5)),

            )):

 

        # 12 vertices

        self.vertices = verts

 

        # 20 equiangular triangles

        self.faces = (

            ('xza','xzb','yzd'),

            ('yzd','xzb','xyb'),

            ('xyb','xzb','xya'),

            ('xya','yza','xzb'),

            ('xzb','yza','xza'),

 

            ('xzd','xzc','yzb'),

            ('yzb','xzd','xyd'),

            ('xyd','xzd','xyc'),

            ('xyc','xzd','yzc'),

            ('yzc','xzd','xzc'),

 

            ('xyd','yzb','yza'),

            ('yza','yzb','xya'),

            ('xya','yzb','xzc'),

            ('xzc','xya','xyb'),

            ('xyb','xzc','yzc'),

            ('yzc','xyb','yzd'),

            ('yzd','yzc','xyc'),

            ('xyc','yzd','xza'),

            ('xza','xyc','xyd'),

            ('xyd','xza','yza')

            )

           

 

        self.edges = self._distill()

        

        self.rectangles = (

            ('xya','xyb','xyc','xyd'),

            ('xza','xzb','xzc','xzd'),

            ('yza','yzb','yzc','yzd'))    

 

    def goldrects(self):

        Edge.color = green

        for r in self.rectangles:

            c0,c1,c2,c3 = [self.vertices[i] for i in r]

            Edge(c0,c1).draw()

            Edge(c1,c2).draw()

            Edge(c2,c3).draw()

            Edge(c3,c0).draw()

 

class Cuboctahedron (Polyhedron):

    pass

 

def test():

    """

    The Concentric Hierarchy by R. Buckminster Fuller

    """

    Edge.color = color.orange

    tetra = Tetrahedron() * 0.5

    tetra.draw()

    

    Edge.color = color.green

    cube = Cube() * 0.5

    cube.draw()

 

    Edge.color = color.red

    cube = Octahedron()

    cube.draw()

 

    Edge.color = color.cyan

    ico = Icosahedron() * sqrt(2)

    ico.draw()

 

def test2():

    """

    Coupler in a Cube (canonical volumes 1 and 3 respectively)

    """

 

    Edge.color = color.orange

    tetra = Tetrahedron()

    tetra.draw()

    

    Edge.color = color.blue

    coupler = Mite()

    coupler.draw()

 

    #Edge.color = color.blue

    #coupler = Coupler()

    #coupler.draw()

    

    Edge.color = color.green

    cube = Cube()

    cube.draw()    

 

if __name__ == '__main__':

    test()

    # test2()