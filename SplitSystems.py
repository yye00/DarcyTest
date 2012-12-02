# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 19:38:53 2012

@author: yye00
"""

from cbc.pdesys import *
#set_log_active(True)
#set_log_level(5)

# Create a mesh
mesh = UnitSquare(16, 16,"crossed")
# and plot it
# plot(mesh, interactive=True)

# setup the BCs
def left(x):   return x[0] < DOLFIN_EPS
def left2(x):   return x[0] < 4*DOLFIN_EPS
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
problem_parameters['time_integration'] = "Transient" # default='Steady'
problem = Problem(mesh, problem_parameters)

# Set up first PDESystem
solver_parameters['iteration_type'] = 'Newton'
solver_parameters['familyname'] = 'Scalar'
solver_parameters['space']['u'] = VectorFunctionSpace #
default=FunctionSpace
solver_parameters['degree']['u'] = 2 # default=1
solver_parameters['degree']['p'] = 1 # default=1
NStokes = PDESystem([['u', 'p']], problem, solver_parameters)

class NavierStokes(PDESubSystem):
    def form(self, u, v_u, u_, u_1, p, v_p, dt, **kwargs):
        return inner(u,v_u)*dx - inner(p, div(v_u))*dx   \
           + inner(div(u),v_p)*dx

noslip = Constant((0.0, 0.0))
inlet = Expression(('6.0*x[1]*(1.0-x[1])','0.0'))

bc = [DirichletBC(NStokes.V['up'].sub(0), inlet,  left),
      DirichletBC(NStokes.V['up'].sub(1), Constant(0.0), right),
      DirichletBC(NStokes.V['up'].sub(0), noslip, top),
      DirichletBC(NStokes.V['up'].sub(0), noslip, bottom)]

NStokes.pdesubsystems['up'] = NavierStokes(vars(NStokes), ['u', 'p'], bcs=bc, reassemble_lhs=False)

def update(self):
    plot(self.pdesystems['Scalar'].u_, title="u", scale=0.5)
    plot(self.pdesystems['Scalar'].p_, title="p", scale=1.0)

Problem.update = update

# Integrate the solution from t=0 to t=0.5
problem.prm['T'] = 0.1
problem.solve()

solver_parameters['familyname'] = 'Scalar2'
solver_parameters['degree']['Sw'] = 1 # default=1
solver_parameters['family']['Sw'] = "DG"
solver_parameters['iteration_type'] = 'Picard'

scalar = PDESystem([['Sw']], problem, solver_parameters)
                
class Scalar(PDESubSystem):
    def form(self, U_, dt, Sw, Sw_, Sw_1, v_Sw, **kwargs):
        Krw = Sw_
        Kro = 1.0 - Sw_
        fw = Krw/(Krw+Kro)
        F =  1.0/dt*inner((Sw-Sw_1),v_Sw)*dx - inner(div(fw*U_),v_Sw)*dx
        return F

bcSw = [DirichletBC(scalar.V['Sw'], Constant(1.0), left2)]

scalar.U_ = 0.5*(NStokes.u_+NStokes.u_1)

csub1 = Scalar(vars(scalar), ['Sw'], bcs=bcSw, max_inner_iter=5) # Iterate on c_
scalar.pdesubsystems['Sw'] = csub1

# Integrate both PDESystems from t=0.5 to t=1.0 using Picard
# iterations on each time step
def update2(self):
    #import pdb
    #pdb.set_trace()
    plot(self.pdesystems['Scalar2'].Sw_)
Problem.update = update2

# initialization has to be done this way according to CBC
problem.q0 = {'Sw': Expression(('0.5'), element=scalar.V['Sw'].ufl_element())}
problem.initialize(scalar)

# Integrate both PDESystems from T=1.0 to T=1.5 using Newton
# iterations on each time step for the scalar
problem.prm['T'] = 1.5
problem.solve()