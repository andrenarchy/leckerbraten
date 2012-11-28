from dolfin import *
import numpy

#parameters['linear_algebra_backend'] = 'uBLAS'

mesh = Mesh('msh_pan.xml')
subdomains = MeshFunction('uint', mesh, 'msh_pan_physical_region.xml')

# Manually specify this. In the future,
# Dolfin will be able to do this automatically.
subdomain_classification = {1: 'pan',
                            2: 'handle',
                            3: 'steak'
                            }


# Setting heat conductivity.
lambda_values = {'pan': 1.0, #1.172e-5,
                 'handle': 1.0, # 8.2e-8,
                 # steak: cf. P.S. Sheridana, N.C. Shilton, http://www.sciencedirect.com/science/article/pii/S0260877401000838
                 'steak': 1.0 #1.5e-7   
                 }
# Define a custom Expression().
class Diffusivity(Expression):
    def eval_cell(self, values, x, cell):
        subdomain_id = subdomains.array()[cell.index]
        clss = subdomain_classification[subdomain_id]
        values[0] = lambda_values[clss]
        return

# only needed for plotting
plot = True
if plot:
    V0 = FunctionSpace(mesh, 'DG', 0) 
    lmbda = Function(V0)
    lmbda.interpolate(Diffusivity())
    File('cont.pvd') << lmbda

mesh = UnitCube(10,10,10)
V = FunctionSpace(mesh, 'CG', 1)

u0 = Constant(350.0)
def u0_boundary(x, on_boundary):
    # TODO: use physical region
    return on_boundary and x[2]<1e-12
bc = DirichletBC(V, u0, u0_boundary)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

# do the actual time stepping
T = 1.0
dt = 1.0e-5
w = Function(V)
w.interpolate(Constant(293.0))
#w.interpolate(Expression('350 - 100*x[2]'))

# this works:
#w.interpolate(Expression('293+57*exp(-1*x[2])'))
# this doesn't work:
#w.interpolate(Expression('293+57*exp(-100*x[2])'))

wfile = File('solution.pvd')
wfile << w
t = dt

# Previous time step.
w_1 = Function(V)
w_1.assign(w)

# boundary condition:
# n \dot grad(u) = alpha * (u - Tenv)
alpha = 1
Tenv = 350

# backward Euler
lambd = Diffusivity()
a = u*v*dx + dt*inner( grad(u), grad(v))*dx + dt*alpha*u*v*ds(0)
L = w_1*v*dx  + dt*alpha*Tenv*v*ds # rhs == 0

# reactionary (forward) Euler
#a = u*v*dx 
#L = w_1*v*dx - dt*inner( grad(w_1), grad(v))*dx

while t <= T:
    print t
    # Homogenous boundary conditions don't require changing the system.
    prm = parameters['krylov_solver'] # short form
    prm['absolute_tolerance'] = 1E-10
    prm['relative_tolerance'] = 1E-14
    prm['maximum_iterations'] = 100
    prm['monitor_convergence'] = True
    #solve(M, y.vector(), b,
    solve(a == L, w,
          #bcs = bc,
          )#solver_parameters={'linear_solver': 'cg',
           #                  'preconditioner': 'amg'})

    # Write the data out to files.
    #plot(y)
    wfile << w
    t += dt
    w_1.assign(w)
