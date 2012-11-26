# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 15:22:46 2012

@author: yye00
"""


Soc = 0.1
Swirr = 0.2

from cbc.pdesys import *
set_log_active(True)
set_log_level(5)

#set_log_active(True)
#set_log_level(5)

mesh = UnitSquare(20, 20)
h = CellSize(mesh)

print(DOLFIN_EPS)
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
structure.mark()

# Extract sub meshes
fluid_mesh = SubMesh(mesh, sub_domains, 0)
structure_mesh = SubMesh(mesh, sub_domains, 1)

# Plot meshes
plot(fluid_mesh, title="Fluid")
plot(structure_mesh, title="Structure")
interactive()



# add some optimization options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

# problem_parameters are defined in Problem.py
#problem_parameters['dt'] = 0.001
problem_parameters['time_integration'] = "Transient" # default='Steady'
problem = Problem(mesh, problem_parameters)

# Set up first PDESystem
solver_parameters['space']['u'] = VectorFunctionSpace #
default=FunctionSpace
solver_parameters['degree']['u'] = 2 # default=1
solver_parameters['degree']['p'] = 1 # default=1
solver_parameters['degree']['Sw'] = 1 # default=1
#solver_parameters['family']['Sw'] = 'DG'

solver_parameters['familyname'] = 'Scalar'
solver_parameters['iteration_type'] = 'Newton'
GloalFormulation = PDESystem([['u', 'p', 'Sw']], problem, solver_parameters)

class DarcyGlobal(PDESubSystem):
    def form(self, u, u_, u_1, v_u, p, v_p, dt, Sw, Sw_, Sw_1, v_Sw, **kwargs):
        # SUPG stablization
        #u_norm = sqrt(dot(u,u))
        #r =
        Krw = Sw_
        Kro = 1.0 - Sw_
        Dlambda = Krw+Kro
        fw = Krw/(Krw+Kro)
        F =  inner(u,v_u)*dx - Dlambda*inner(p, div(v_u))*dx   \
           + inner(div(u),v_p)*dx      \
           + 1.0/dt*inner((Sw-Sw_1),v_Sw)*dx - inner(div(fw*u),v_Sw)*dx

        return F


ufile  = File("results/velocity.pvd")
pfile  = File("results/pressure.pvd")
Swfile = File("results/Sw.pvd")

def update(self):
    plot(self.pdesystems['Scalar'].Sw_, title="Sw", scale=1.0)
    plot(self.pdesystems['Scalar'].u_, title="u", scale=0.05)
    plot(self.pdesystems['Scalar'].p_, title="p", scale=1.0)
    # Save to file
    ufile  << self.pdesystems['Scalar'].u_
    pfile  << self.pdesystems['Scalar'].p_
    Swfile << self.pdesystems['Scalar'].Sw_


Problem.update = update

problem.q0 = {'upSw': Expression(('0', '0', '0', '0.0'), element=GloalFormulation.V['upSw'].ufl_element())}

problem.initialize(GloalFormulation)


noslip = Constant((0.0, 0.0))
#inlet = Expression(('0.5*sin(x[1]*3.141592654)','0.0'))
inlet = Constant((1.0, 0.0))
Sw_inlet = Constant(1.0)
outlet = Constant((1.0,0.0))
bc = [DirichletBC(GloalFormulation.V['upSw'].sub(0), inlet,  left),
      DirichletBC(GloalFormulation.V['upSw'].sub(1), Constant(0.0), right),
      DirichletBC(GloalFormulation.V['upSw'].sub(2), Constant(1.0), left),
      DirichletBC(GloalFormulation.V['upSw'].sub(0), noslip, top),
      DirichletBC(GloalFormulation.V['upSw'].sub(0), noslip, bottom)]

#DirichletBC(GloalFormulation.V['upSw'].sub(2), Constant(0.0), right),
#      DirichletBC(GloalFormulation.V['upSw'].sub(2), 0.0, right)]

GloalFormulation.add_pdesubsystem(DarcyGlobal, ['u', 'p', 'Sw'], bcs=bc)

problem.prm['T'] = 0.5
problem.solve()
print "Done"
