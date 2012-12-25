from dolfin import *
import numpy

def main():
    #parameters['linear_algebra_backend'] = 'uBLAS'

    mesh = Mesh('msh_pan.xml')
    subvolumes = MeshFunction('uint', mesh, 'msh_pan_physical_region.xml')
    #subfacets = MeshFunction('uint', mesh, 'msh_pan_facet_region.xml', )
    boundary_parts = MeshFunction("uint", mesh, mesh.topology().dim()-1)

    class LowerRobinBoundary(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14   # tolerance for coordinate comparisons
            return on_boundary and abs(x[2]) < tol
    Gamma_R = LowerRobinBoundary()
    Gamma_R.mark(boundary_parts, 0)

    class OtherRobinBoundary(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14   # tolerance for coordinate comparisons
            return on_boundary and abs(x[2]) >= tol
    Gamma_O = OtherRobinBoundary()
    Gamma_O.mark(boundary_parts, 1)

    # Manually specify this. In the future,
    # Dolfin will be able to do this automatically.
    subvolumes_classification = {1: 'pan',
                                2: 'handle',
                                3: 'steak'
                                }
    subfacets_classification = {
            0: 'air',
            4: 'heater'
            }

    # Setting heat conductivity.
    lambda_values = {'pan': 1.172e-5,
                     'handle': 8.2e-8,
                     # steak: cf. P.S. Sheridana, N.C. Shilton, http://www.sciencedirect.com/science/article/pii/S0260877401000838
                     'steak': 1.5e-7   
                     }

    tenv_values = { 
            'heater': (lambda t: min(350, 293 + 10*t)),
            'air': (lambda t: 293)
            }

    # Define a custom Expression().
    class Diffusivity(Expression):
        def eval_cell(self, values, x, cell):
            subvolume_id = subvolumes.array()[cell.index]
            clss = subvolumes_classification[subvolume_id]
            values[0] = lambda_values[clss]
            return

    # only needed for plotting
    plot = True
    if plot:
        V0 = FunctionSpace(mesh, 'DG', 0) 
        lmbda = Function(V0)
        lmbda.interpolate(Diffusivity())
        File('cont.pvd') << lmbda

    # mesh = UnitCube(10,10,10)
    V = FunctionSpace(mesh, 'CG', 1)

    bc = DirichletBC(V, Constant(350), boundary_parts, 0)

    #u0 = Constant(350.0)
    #def u0_boundary(x, on_boundary):
    #    # TODO: use physical region
    #    return on_boundary and x[2]<1e-12
    #bc = DirichletBC(V, u0, u0_boundary)
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)

    # do the actual time stepping
    T = 1e2
    dt = 1.0e-1
    w = Function(V)
    w.interpolate(Constant(293.0))
    #w.interpolate(Expression('350 - 100*x[2]'))

    # this works:
    w.interpolate(Expression('293+57*exp(-1*x[2])'))
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
    Tenv = 293.0

    # backward Euler
    lambd = 1e-7 #Diffusivity()
    a = u*v*dx + dt*inner( lambd*grad(u), grad(v))*dx #+ dt*lambd*alpha*u*v*ds

    # reactionary (forward) Euler
    #a = u*v*dx 
    #L = w_1*v*dx - dt*inner( grad(w_1), grad(v))*dx

    class Tenv(Expression):
        def eval_cell(self, values, x, cell):
            subfacet_id = subfacets.array()[cell.index]
            clss = subfacets_classification[subfacet_id]
            values[0] = tenv_values[clss](t)
            return
    tenv = Tenv()

    while t <= T:
        print t
        tenv = min(350, 293 + t*10)
        L = w_1*v*dx # + dt*lambd*alpha*tenv*v*ds(0) + dt*lambd*alpha*293*v*ds(1) # rhs == 0
        # Homogenous boundary conditions don't require changing the system.
        prm = parameters['krylov_solver'] # short form
        prm['absolute_tolerance'] = 1E-10
        prm['relative_tolerance'] = 1E-14
        prm['maximum_iterations'] = 100
        prm['monitor_convergence'] = True
        #solve(M, y.vector(), b,
        solve(a == L, w,
              bcs = bc,
              )#solver_parameters={'linear_solver': 'cg',
               #                  'preconditioner': 'amg'})

        # Write the data out to files.
        #plot(y)
        wfile << w
        t += dt
        w_1.assign(w)

if __name__ == '__main__':
    main()
