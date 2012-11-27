# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 15:22:46 2012

@author: yye00
"""

from cbc.pdesys import *

set_log_active(True)
set_log_level(5)

mesh = UnitSquare(20, 20,"crossed")
h = CellSize(mesh)

# setup the BCs
def left(x):   return x[0] < DOLFIN_EPS
def right(x):  return x[0] > 1.0-DOLFIN_EPS
def top(x):    return x[1] > 1.0-DOLFIN_EPS
def bottom(x): return x[1] < DOLFIN_EPS
def allbounds(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0-DOLFIN_EPS \
        or x[1] < DOLFIN_EPS or x[1] > 1.0-DOLFIN_EPS

# add some optimization options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

# problem_parameters are defined in Problem.py
#problem_parameters['dt'] = 0.001
problem_parameters['time_integration'] = "Steady" # default='Steady'
problem = Problem(mesh, problem_parameters)

# Set up first PDESystem
solver_parameters['space']['u'] = VectorFunctionSpace #
default=FunctionSpace
solver_parameters['degree']['u'] = 2 # default=1
solver_parameters['degree']['p'] = 1 # default=1

solver_parameters['familyname'] = 'Scalar'
solver_parameters['iteration_type'] = 'Newton'
GloalFormulation = PDESystem([['u', 'p']], problem, solver_parameters)

class DarcyGlobal(PDESubSystem):
    def form(self, u, u_, v_u, p, v_p, **kwargs):
        F =  inner(u,v_u)*dx - 0.1*inner(p, div(v_u))*dx + inner(div(u),v_p)*dx
        return F

noslip = Constant((0.0, 0.0))
inlet = Constant((1.0, 0.0))
Sw_inlet = Constant(1.0 )
bc = [DirichletBC(GloalFormulation.V['up'].sub(1), Constant(0.0), right),
      DirichletBC(GloalFormulation.V['up'].sub(0), noslip, top),
      DirichletBC(GloalFormulation.V['up'].sub(0), noslip, bottom),
      DirichletBC(GloalFormulation.V['up'].sub(0), inlet, left)]

GloalFormulation.add_pdesubsystem(DarcyGlobal, ['u', 'p'], bcs=bc)

problem.solve()

plot(problem.pdesystems['Scalar'].u_, interactive=True, title="u", scale=0.1)
plot(problem.pdesystems['Scalar'].p_, interactive=True, title="p")

ufile  = File("results/velocity.pvd")
pfile  = File("results/pressure.pvd")

