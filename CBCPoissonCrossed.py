# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 18:50:34 2012

@author: yye00
"""

from cbc.pdesys import *

mesh = UnitSquare(10, 10,"crossed")
# Change desired items in the problem_parameters dict from cbc.pdesys
problem_parameters['time_integration'] = "Steady" # default='Steady'
problem = Problem(mesh, problem_parameters)

solver_parameters['familyname'] = 'Scalar'
solver_parameters['degree']['u']=1
solver_parameters['family']['u']='CG'
solver_parameters['iteration_type'] = 'Newton'
poisson = PDESystem([['u']], problem, solver_parameters) # Creates FunctionSpace, Functions etc.
poisson.f = Expression('-exp(-(pow(x[0]-0.5,2) + pow(x[1]-0.5,2))/0.2)')

class Poisson(PDESubSystem):
    def form(self, u, v_u, f, **kwargs):    # v_u is the TestFunction
        return inner(grad(u), grad(v_u))*dx + f*v_u*dx

def left(x):   return x[0] < DOLFIN_EPS
def right(x):  return x[0] > 1.0-DOLFIN_EPS
def top(x):    return x[1] > 1.0-DOLFIN_EPS
def bottom(x): return x[1] < DOLFIN_EPS
def allbounds(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0-DOLFIN_EPS \
        or x[1] < DOLFIN_EPS or x[1] > 1.0-DOLFIN_EPS


bcs = [DirichletBC(poisson.V['u'], (0.0), left),
       DirichletBC(poisson.V['u'], (0.0), right),
       DirichletBC(poisson.V['u'], (0.0), top),
       DirichletBC(poisson.V['u'], (0.0), bottom)]

poisson.pdesubsystems['u'] = Poisson(vars(poisson), ['u'], bcs=bcs)
problem.solve()

plot(problem.pdesystems['Scalar'].u_, interactive=True, title="u")

plot(project(mesh, div(problem.pdesystems['Scalar'].u_),
             FunctionSpace(mesh,"CG",1) ), interactive=True, title="div of u")





