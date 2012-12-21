'''
Omega: [0,1]^3
solution: u(x,t) = sin(pi*t)*exp(-norm(x)^2)
=> rhs:   f(x,t) = (pi*cos(pi*t) - sin(pi*t)*(4*norm(x)^2 - 6))*exp(-norm(x)^2)
'''

from dolfin import *
normstr = '(x[0]*x[0]+ x[1]*x[1] + x[2]*x[2])'
u_ex = Expression('sin(pi*t)*exp(-' + normstr + ')', t=0)
f = Expression('(pi*cos(pi*t) - sin(pi*t)*(4*' + normstr + ' - 6))*exp(-' + normstr + ')', t=0)

def solve_heat(n):
    mesh = UnitCube(n, n, n)
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    bc = DirichletBC(V, u_ex, 'on_boundary')

    u_old = Function(V)
    u_old.interpolate(u_ex)
    u_new = Function(V)

    t = 0
    dt = 1./n
    tend = 2

    Avar = u*v*dx + dt*inner(grad(u),grad(v))*dx

    L2errors = []

    while t<tend:
        t += dt
        u_ex.t = t
        f.t = t

        bvar = (u_old + dt*f)*v*dx
        solve( Avar == bvar, u_new, bc)

        u_old = u_new.copy()
        L2errors += [errornorm(u_ex,u_new)]

    return max(L2errors)

errors = []
P = range(1,5)
for p in P:
    errors += [solve_heat(2**p)]

from matplotlib import pyplot as pp
pp.loglog( [1./(2**p) for p in P], errors , '-o')
pp.loglog( [1./(2**P[0]), 1./(2**P[-1])], [errors[0], errors[0]*(1./(2**(P[-1]-P[0])))])
pp.show()
