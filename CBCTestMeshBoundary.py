# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 15:22:46 2012

@author: yye00
"""


from cbc.pdesys import *
set_log_active(True)
set_log_level(5)


mesh = UnitSquare(20, 20)
h = CellSize(mesh)

# setup the BCs
def left(x, on_boundary):
    print "Left:", x
    if on_boundary :
        print "on_boundary"
    if x[0] < DOLFIN_EPS:
        print "boundary"
    else:
        print "inside"
    return on_boundary and x[0] < 1.0 / 20 - DOLFIN_EPS
def right(x):  return x[0] > 1.0-DOLFIN_EPS
def top(x):    return x[1] > 1.0-DOLFIN_EPS
def bottom(x): return x[1] < DOLFIN_EPS
def allbounds(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0-DOLFIN_EPS \
        or x[1] < DOLFIN_EPS or x[1] > 1.0-DOLFIN_EPS


class Structure(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > DOLFIN_EPS and x[0] < 1.0 - DOLFIN_EPS \
           and x[1] > DOLFIN_EPS and x[1] < 1.0 - DOLFIN_EPS

sub_domains = MeshFunction("uint", mesh, mesh.topology().dim())
sub_domains.set_all(0)

# Mark structure domain as 1
structure = Structure()
structure.mark(sub_domains, 1)

# Extract sub meshes
boundary_mesh = SubMesh(mesh, sub_domains, 0)
structure_mesh = SubMesh(mesh, sub_domains, 1)

# Plot meshes
plot(boundary_mesh, title="Boundary")
plot(structure_mesh, title="Structure... or the inside")
interactive()
