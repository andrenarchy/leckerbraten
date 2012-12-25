from dolfin import *
import numpy

# robin boundary -dot(grad(u),n) = alpha*(u - q)
class RobinBC:
    def __init__(self, alpha, q):
        self.alpha = alpha
        self.q = q

def solve_heat(mesh, u_init=Constant(0.0), lambd=Constant(1.0), f=Constant(0.0), boundaries=None, dbcs={}, rbcs={}, expressions=[], scale_dt=1.0, u_ex=None):
    '''Solve the heat equation.

    Solve the heat equation
        u_t - div(lambd * grad(u)) = f          in omega
    with Dirichlet boundary conditions
                                 u = dbcs[tag]  on boundary tagged as 'tag'
    and Robin boundary conditions
                   -dot(grad(u),n) = rbcs[tag].alpha*(u - rbcs[tag].q) on boundary tagged as 'tag'
    Arguments:
        mesh:   the discretized domain.
        u_init: initial temperature distribution.
        lambd:  diffusivity.
        f:      right hand side.
        boundaries: MeshFunction with boundary tags.
                If boundaries is None, the default tag is 0.
        dbcs:   dictionary mapping boundary tags to the corresponding 
                Dirichlet boundary function.
        rbcs:   dictionary mapping boundary tags to the corresponding 
                RobinBC objects.
        expressions: list of expressions that need the time to be set in 
                the time loop.
        scale_dt: scale of timestep, dt = scale_dt*mesh.hmax().
        u_ex:   exact solution (if known).
    '''

    if boundaries is None:
        boundaries = MeshFunction("uint", mesh, mesh.topology().dim()-1)
        boundaries.set_all(0)
    if len(dbcs)==0 and len(rbcs)==0:
        raise ValueError('Dirichlet or Robin boundary conditions have to be specified (both possible)!')

    # keep a list of all expressions that need the current time
    expressions = list(expressions)
    expressions += [lambd, f]
    for _, bc in dbcs.items():
        expressions += [bc]
    for _, bc in rbcs.items():
        expressions += [bc.alpha]
        expressions += [bc.q]
    if u_ex is not None:
        expressions += [u_ex]
    # set time in ALL THE expressions
    def set_time(expressions, t):
        for expr in expressions:
            expr.t = t
    
    # make sure you use high quality meshes, 
    # i.e. hmin/hmax should be close to 1. 
    # the following command can be used for gmsh:
    #     gmsh -clmax 0.1 -3 -optimize msh_pan.geo
    hmax = mesh.hmax()
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    
    ds = Measure("ds")[boundaries]
    dbc = [DirichletBC(V, bc, boundaries, tag) for tag, bc in dbcs.items()]

    t = 0.
    set_time(expressions, t)
    dt = scale_dt*hmax
    tend = 1.

    print("Solve with hmax=%e (%d unknowns) and dt=%e (hmin/hmax=%e)." % (hmax, V.dim(), dt, mesh.hmin()/hmax))

    Avar = u*v*dx + dt*inner(lambd*grad(u),grad(v))*dx
    # incorporate Robin boundary conditions into bilinear form
    for tag, bc in rbcs.items():
        Avar += dt*lambd*bc.alpha*u*v*ds(tag)

    L2errors = []

    u_old = Function(V)
    u_old.interpolate(u_init)
    u_new = Function(V)
    u_new = u_old.copy()

    while t<tend:
        t += dt
        set_time(expressions, t)

        bvar = (u_old + dt*f)*v*dx
        for tag, bc in rbcs.items():
            bvar += dt*lambd*bc.alpha*bc.q*v*ds(tag)

        A, b = assemble_system(Avar, bvar, dbc, exterior_facet_domains=boundaries)

        solve(A, u_new.vector(), b, "cg", "amg")

        u_old = u_new.copy()
        if u_ex is not None:
            L2errors += [errornorm(u_ex,u_new)]

    ret = {
            "u": u_new,
            "L2errors": L2errors
            }
    return ret

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
