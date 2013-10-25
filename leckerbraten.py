from dolfin import *
import numpy

# =============================================================================
# robin boundary -dot(grad(u),n) = alpha*(u - q)
class RobinBC:
    def __init__(self, alpha, q):
        self.alpha = alpha
        self.q = q
# =============================================================================
def solve_heat(mesh,
               u_init = Constant(0.),
               lambd = Constant(1.),
               f = Constant(0.),
               boundaries = None,
               dbcs = {},
               rbcs = {},
               expressions = [],
               t0 = 0.0,
               tend = 1.0,
               scale_dt = 1.0,
               wfile = None,
               u_ex = None
               ):
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
        boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
        boundaries.set_all(0)
    if not dbcs and not rbcs:
        raise ValueError('Dirichlet or Robin boundary conditions have to be specified (both possible)!')

    # keep a list of all expressions that need the current time
    expressions = list(expressions)
    expressions.append(lambd)
    expressions.append(f)
    expressions.extend(dbcs.values())
    for _, bc in rbcs.items():
        expressions += [bc.alpha]
        expressions += [bc.q]
    if u_ex is not None:
        expressions.append(u_ex)
    # Set time in ALL THE expressions.
    # TODO this breaks if any of the expressions has no .t available. Fix this.
    def set_time(expressions, t):
        for expr in expressions:
            # TODO Well, well. If any of the expressions has a Sum(), in which
            #      one component has the parameter 't', that parameter is not
            #      going to get set. With hasattr(), this fails silently,
            #      otherwise an exception would be raised. This is something
            #      that needs to be fixed in any case.
            if hasattr(expr, 't'):
                expr.t = t

    # Make sure you use high quality meshes,
    # i.e. hmin/hmax should be close to 1.
    # The following command can be used for gmsh:
    #     gmsh -clmax 0.1 -3 -optimize msh_pan.geo
    hmax = mesh.hmax()
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    ds = Measure('ds')[boundaries]
    dbc = [DirichletBC(V, bc, boundaries, tag) for tag, bc in dbcs.items()]

    t = t0
    set_time(expressions, t)
    dt = scale_dt * hmax

    print('Solve with hmax=%e (%d unknowns) and dt=%e (hmin/hmax=%e).'
          % (hmax, V.dim(), dt, mesh.hmin()/hmax)
         )

    Avar = u*v*dx + dt * inner(lambd*grad(u),grad(v))*dx
    # incorporate Robin boundary conditions into bilinear form
    for tag, bc in rbcs.items():
        Avar += dt * lambd * bc.alpha * u * v * ds(tag)

    L2errors = []

    u_old = Function(V)
    u_old.interpolate(u_init)

    u_new = Function(V, name='temperature')
    u_new.assign(u_old)
    if wfile is not None:
        wfile << (u_new, t)

    while t < tend:
        print('Time step %e -> %e   (dt = %e)...' % (t, t+dt, dt))
        t += dt
        set_time(expressions, t)

        bvar = (u_old + dt*f) * v * dx
        for tag, bc in rbcs.items():
            bvar += dt * lambd * bc.alpha * bc.q * v * ds(tag)

        A, b = assemble_system(Avar, bvar, dbc,
                               exterior_facet_domains = boundaries
                               )

        solve(A, u_new.vector(), b, 'cg', 'amg')

        u_old.assign(u_new)
        if u_ex is not None:
            L2errors.append(errornorm(u_ex, u_new))
        if wfile is not None:
            wfile << (u_new, t)

    ret = {'u': u_new,
           'L2errors': L2errors
           }
    return ret
# =============================================================================
def main():
    mesh = Mesh('msh_pan.xml')
    subvolumes = MeshFunction('size_t', mesh, 'msh_pan_physical_region.xml')
    #subfacets = MeshFunction('size_t', mesh, 'msh_pan_facet_region.xml', )
    boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    boundaries.set_all(0)
    rbcs = {
            0: RobinBC(Constant(1.), Constant(293.))
            }

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[2]<1e-12
    bottom = Bottom()
    bottom.mark(boundaries, 1)

    #class Tbottom(Expression):
    #    def eval(self, value, x):
    #        value[0] = 293 + min(self.t*10, 10)
    rbcs[1] = RobinBC(Constant(1.), Expression('293.0 + 10.0*t', t=0))

    # Manually specify this. In the future,
    # Dolfin will be able to do this automatically.
    subvolumes_classification = {1: 'pan',
                                 2: 'handle',
                                 3: 'steak'
                                 }

    # Setting heat conductivity.
    lambd_values = {
            'pan': 1.172e-5,
            'handle': 8.2e-8,
            # steak: cf. P.S. Sheridana, N.C. Shilton,
            # http://www.sciencedirect.com/science/article/pii/S0260877401000838
            'steak': 1.5e-7
            }

    # Define a custom Expression().
    class Lambd(Expression):
        def eval_cell(self, values, x, cell):
            subvolume_id = subvolumes.array()[cell.index]
            clss = subvolumes_classification[subvolume_id]
            values[0] = lambd_values[clss]
            #values[0] = 1.0
            return

    wfile = XDMFFile('solution.xdmf')
    wfile.parameters['flush_output'] = True;
    wfile.parameters['rewrite_function_mesh'] = False;
    sol = solve_heat(mesh,
            u_init=Constant(293.),
            lambd=Lambd(),
            boundaries=boundaries,
            rbcs=rbcs,
            tend=30.,
            scale_dt=0.5,
            wfile=wfile
            )
    return
# =============================================================================
if __name__ == '__main__':
    main()
# =============================================================================
