'''
Omega: [0,1]^3
solution: u(x,t) = sin(pi*t)*exp(-norm(x)^2)
=> rhs:   f(x,t) = (pi*cos(pi*t) - sin(pi*t)*(4*norm(x)^2 - 6))*exp(-norm(x)^2)
'''

from dolfin import *
import numpy as np
from matplotlib import pyplot as pp

normstr = '(x[0]*x[0]+ x[1]*x[1] + x[2]*x[2])'
u_ex = Expression('sin(pi*t)*exp(-' + normstr + ')', t=0)
f = Expression('(pi*cos(pi*t) - sin(pi*t)*(4*' + normstr + ' - 6))*exp(-' + normstr + ')', t=0)
u_ex_grad = Expression( (
        '-2*sin(pi*t)*x[0]*exp(-' + normstr + ')',
        '-2*sin(pi*t)*x[1]*exp(-' + normstr + ')',
        '-2*sin(pi*t)*x[2]*exp(-' + normstr + ')'
        ), t=0)
# robin boundary -dot(grad(u),n) = alpha*(u - q)
alpha = Constant(1.0)

class RobinBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0]<0.5

mesh = UnitCube(2,2,2)
#mesh = Mesh('msh_pan.xml')

def solve_heat(mesh, p, boundaries=None, scale_dt=1.0):
    hmax = mesh.hmax()
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    if boundaries is None:
        boundaries = MeshFunction("uint", mesh, mesh.topology().dim()-1)
        boundaries.set_all(1)
    # tag 0: Dirichlet
    # tag 1: Robin
    ds = Measure("ds")[boundaries]
    bc = DirichletBC(V, u_ex, boundaries, 0)

    n = FacetNormal(mesh)

    t = 0.
    u_ex.t = t
    f.t = t
    dt = scale_dt/(2**p)
    tend = 1.

    print("Solve with hmax=%e (%d unknowns) and dt=%e." % (hmax, V.dim(), dt))

    Avar = u*v*dx + dt*inner(grad(u),grad(v))*dx + dt*alpha*u*v*ds(1)
    # Robin boundary function
    q = u_ex + (1/alpha)*dot(u_ex_grad, n)

    L2errors = []

    u_old = Function(V)
    u_old.interpolate(u_ex)
    u_new = Function(V)

    while t<tend:
        t += dt
        u_ex.t = t
        f.t = t
        u_ex_grad.t = t

        bvar = (u_old + dt*f)*v*dx + dt*alpha*q*v*ds(1)

        A, b = assemble_system(Avar, bvar, bc, exterior_facet_domains=boundaries)

        solve(A, u_new.vector(), b, "cg", "amg")

        u_old = u_new.copy()
        L2errors += [errornorm(u_ex,u_new)]

    return hmax, max(L2errors)

hmaxs = []
errors = []
P = range(0,3)
for p in P:
    # mark Robin boundary condition
    boundaries = MeshFunction("uint", mesh, mesh.topology().dim()-1)
    boundaries.set_all(0)
    robin = RobinBoundary()
    robin.mark(boundaries, 1)

    hmax, err = solve_heat(mesh, p, boundaries=boundaries, scale_dt=0.2)
    hmaxs += [hmax]
    errors += [err]
    mesh = refine(mesh)

hmaxs = np.array(hmaxs)
errors = np.array(errors)
eoc = np.log(errors[1:]/errors[:-1])/np.log(hmaxs[1:]/hmaxs[:-1])
print(errors)
print(eoc)

pp.loglog( hmaxs, errors , '-o')
pp.loglog( [hmaxs[0],hmaxs[-1]], [errors[0], errors[0]*(hmaxs[-1]/hmaxs[0])])
pp.show()
