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

mesh = UnitCube(2,2,2)
# mesh = Mesh('msh_pan.xml')

def solve_heat(mesh, p):
    hmax = mesh.hmax()
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    bc = DirichletBC(V, u_ex, 'on_boundary')

    t = 0
    u_ex.t = 0
    f.t = 0
    dt = 1./(2**p)
    tend = 2

    Avar = u*v*dx + dt*inner(grad(u),grad(v))*dx

    L2errors = []

    u_old = Function(V)
    u_old.interpolate(u_ex)
    u_new = Function(V)

    while t<tend:
        t += dt
        u_ex.t = t
        f.t = t
        print(t)

        bvar = (u_old + dt*f)*v*dx
        solve( Avar == bvar, u_new, bc)

        u_old = u_new.copy()
        L2errors += [errornorm(u_ex,u_new)]

    return hmax, max(L2errors)

hmaxs = []
errors = []
P = range(0,4)
for p in P:
    hmax, err = solve_heat(mesh, p)
    hmaxs += [hmax]
    errors += [err]
    mesh = refine(mesh)

hmaxs = np.array(hmaxs)
errors = np.array(errors)
eoc = np.log(errors[1:]/errors[:-1])/np.log(hmaxs[1:]/hmaxs[:-1])
print(eoc)

pp.loglog( hmaxs, errors , '-o')
pp.loglog( [hmaxs[0],hmaxs[-1]], [errors[0], errors[0]*(hmaxs[-1]/hmaxs[0])])
pp.show()
