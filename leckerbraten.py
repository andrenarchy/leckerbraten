from dolfin import *
import numpy

mesh = Mesh('coils3d.xml')
subdomains = MeshFunction('uint', mesh, 'coils3d_physical_region.xml')

# Manually specify this. In the future,
# Dolfin will be able to do this automatically.
subdomain_classification = {1: 'coil',
                            2: 'coil',
                            3: 'brick',
                            4: 'air'
                            }

# Setting heat conductivity.
V0 = FunctionSpace(mesh, 'DG', 0)
lambda_values = {'coil': 1.0,
                 'brick': 0.1,
                 'air': 0.8
                 }
# Define a custom Expression().
class Conductivity(Expression):
    def eval_cell(self, values, x, cell):
        subdomain_id = subdomains.array()[cell.index]
        clss = subdomain_classification[subdomain_id]
        values[0] = lambda_values[clss]
        return
lmbda = Function(V0)
lmbda.interpolate(Conductivity())
File('cont.pvd') << lmbda

# Define a custom Expression() for the heater.
heater_values = {'coil': 1.0,
                 'brick': 0.0,
                 'air': 0.0
                 }
class Heater(Expression):
    def eval_cell(self, values, x, cell):
        subdomain_id = subdomains.array()[cell.index]
        clss = subdomain_classification[subdomain_id]
        values[0] = heater_values[clss]
        return
heater = Function(V0)
heater.interpolate(Heater())

V = FunctionSpace(mesh, 'CG', 1)

File('cont.pvd') << lmbda
u0 = Constant(0.0)
def u0_boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, u0, u0_boundary)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = heater

# do the actual time stepping
T = 0.5
dt = 1.0e-2
y = project(Constant(0.0), V)
yfile = File('solution.pvd')
yfile << y
t = dt

# Previous time step.
y_1 = Function(V)
y_1.assign(y)

# backward Euler
a = u*v*dx + dt*inner(lmbda * nabla_grad(u), nabla_grad(v))*dx
L = y_1*v*dx + dt*f*v*dx

while t <= T:
    print t
    # Homogenous boundary conditions don't require changing the system.
    prm = parameters['krylov_solver'] # short form
    prm['absolute_tolerance'] = 1E-10
    prm['relative_tolerance'] = 1E-6
    prm['maximum_iterations'] = 100
    prm['monitor_convergence'] = False
    #solve(M, y.vector(), b,
    solve(a == L, y,
          solver_parameters={'linear_solver': 'cg',
                             'preconditioner': 'amg'})

    # Write the data out to files.
    #plot(y)
    yfile << y
    t += dt
    y_1.assign(y)
