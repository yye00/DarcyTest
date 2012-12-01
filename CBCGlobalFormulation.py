# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 16:59:45 2012

@author: yye00
"""


from cbc.pdesys import *
mesh = UnitSquare(20, 20)

# setup the BCs
def left(x):   return x[0] < DOLFIN_EPS
def right(x):  return x[0] > 1.0-DOLFIN_EPS
def top(x):    return x[1] > 1.0-DOLFIN_EPS
def bottom(x): return x[1] < DOLFIN_EPS
def allbounds(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0-DOLFIN_EPS \
        or x[1] < DOLFIN_EPS or x[1] > 1.0-DOLFIN_EPS

# point source simulation function
qw=Expression("+1.0*exp(-(pow(x[0]-0.25,2) + pow(x[1]-0.25,2))/0.02)")
qo=Expression("-1.0*exp(-(pow(x[0]-0.75,2) + pow(x[1]-0.75,2))/0.02)")

num_refinements = 3
radius = 0.05
for i in range(num_refinements):
    # Mark cells for refinement
    markers = MeshFunction("bool", mesh, mesh.topology().dim())
    markers.set_all(False)
    for cell in cells(mesh):
        if cell.midpoint().distance(Point(0.25,0.25)) < 2*radius or \
           cell.midpoint().distance(Point(0.75,0.75)) < 2*radius :
            markers[cell.index()] = True
    # Refine mesh
    mesh = refine(mesh, markers)



# problem_parameters are defined in Problem.py
problem_parameters['time_integration'] = "Transient" # default='Steady'
problem = Problem(mesh, problem_parameters)

# Set up first PDESystem
solver_parameters['space']['u'] = VectorFunctionSpace #
default=FunctionSpace
solver_parameters['degree']['u'] = 2 # default=1
solver_parameters['degree']['p'] = 1 # default=1
solver_parameters['degree']['Sw'] = 1 # default=1

solver_parameters['familyname'] = 'Scalar'
solver_parameters['iteration_type'] = 'Newton'
GloalFormulation = PDESystem([['u', 'p', 'Sw']], problem, solver_parameters)

class DarcyGlobal(PDESubSystem):
    def form(self, u, v_u, p, v_p, dt, Sw, Sw_, Sw_1, v_Sw, **kwargs):
        Krw = Sw_
        Kro = 1.0-Sw_
        Dlambda = Krw+Kro
        fw = Krw/(Krw+Kro)
        F =  inner(u,v_u)*dx - Dlambda*inner(p, div(v_u))*dx   \
           + inner(div(u),v_p)*dx - (qw+qo)*v_p*dx           \
           + inner(Sw-Sw_1, v_Sw)/dt*dx - inner(div(fw*u),v_Sw)*dx - qw*v_Sw*dx
        return F



problem.prm['T'] = 3.0

ufile  = File("results/velocity.pvd")
pfile  = File("results/pressure.pvd")
Swfile = File("results/Sw.pvd")

def update(self):
    plot(self.pdesystems['Scalar'].Sw_, title="Sw", scale=1.0)
    plot(self.pdesystems['Scalar'].u_, title="u", scale=2.0)
    plot(self.pdesystems['Scalar'].p_, title="p", scale=1.0, interactive=True)
    # Save to file
    ufile  << self.pdesystems['Scalar'].u_
    pfile  << self.pdesystems['Scalar'].p_
    Swfile << self.pdesystems['Scalar'].Sw_

Problem.update = update

# initialization has to be done this way according to CBC
problem.q0 = {'upSw': Expression(('0', '0', '0', '0.5'), element=GloalFormulation.V['upSw'].ufl_element())}

problem.initialize(GloalFormulation)


noslip = Constant((0.0, 0.0))
#inlet = Expression(('6.0*x[1]*(1.0-x[1])','0.0'))
#inlet = Expression(('0.5*sin(x[1]*3.141592654)','0.0'))
#inlet = Constant((1.0, 0.0))
Sw_inlet = Constant(1.0)
outlet = Constant((1.0,0.0))
bc = [DirichletBC(GloalFormulation.V['upSw'].sub(0), noslip,  allbounds),
      DirichletBC(GloalFormulation.V['upSw'].sub(1), Constant(0.0), allbounds)]


GloalFormulation.add_pdesubsystem(DarcyGlobal, ['u', 'p', 'Sw'])


problem.solve()
print "Done"
