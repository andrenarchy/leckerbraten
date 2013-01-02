from dolfin import *
import numpy as np
from matplotlib import pyplot as pp
from leckerbraten import solve_heat, RobinBC

'''
Omega: [0,1]^3
solution: u(x,t) = sin(pi*t)*exp(-norm(x)^2)
=> rhs:   f(x,t) = (pi*cos(pi*t) - sin(pi*t)*( (1+norm(x)^2)*(4*norm(x)^2 - 6) - 4*norm(x)^2 )*exp(-norm(x)^2)
'''
normstr = '(x[0]*x[0]+ x[1]*x[1] + x[2]*x[2])'
# heat diffusivity
lambd = Expression('1 + ' + normstr)
u_ex = Expression('sin(pi*t)*exp(-' + normstr + ')', t=0)
f = Expression('(pi*cos(pi*t) - sin(pi*t)*( (1 + '+normstr+')*(4*' + normstr + ' - 6) - 4*'+normstr+' ) )*exp(-' + normstr + ')', t=0)
u_ex_grad = Expression( (
        '-2*sin(pi*t)*x[0]*exp(-' + normstr + ')',
        '-2*sin(pi*t)*x[1]*exp(-' + normstr + ')',
        '-2*sin(pi*t)*x[2]*exp(-' + normstr + ')'
        ), t=0)

class RobinBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0]<0.5

mesh = UnitCubeMesh(2,2,2)
#mesh = Mesh('msh_pan.xml')
hmaxs = []
errors = []
P = range(0,3)
for p in P:
    # mark Dirichlet and Robin boundary condition parts
    boundaries = MeshFunction("sizet", mesh, mesh.topology().dim()-1)
    boundaries.set_all(0)

    # Dirichlet
    dbcs = {
            0: u_ex
            }

    # Robin
    robin = RobinBoundary()
    robin.mark(boundaries, 1)
    n = FacetNormal(mesh)
    alpha = Constant(1.)
    rbcs = {
            1: RobinBC(alpha, u_ex + (1/alpha)*dot(u_ex_grad, n))
            }

    sol = solve_heat(mesh,
            u_init=u_ex,
            lambd=lambd,
            f=f,
            boundaries=boundaries,
            dbcs=dbcs,
            rbcs=rbcs,
            expressions=[u_ex, u_ex_grad],
            scale_dt=0.2,
            u_ex=u_ex)

    errors += [ max(sol["L2errors"]) ]
    hmaxs += [mesh.hmax()]
    mesh = refine(mesh)

hmaxs = np.array(hmaxs)
errors = np.array(errors)
eoc = np.log(errors[1:]/errors[:-1])/np.log(hmaxs[1:]/hmaxs[:-1])
print(errors)
print(eoc)

pp.loglog( hmaxs, errors , '-o')
pp.loglog( [hmaxs[0],hmaxs[-1]], [errors[0], errors[0]*(hmaxs[-1]/hmaxs[0])])
pp.show()
