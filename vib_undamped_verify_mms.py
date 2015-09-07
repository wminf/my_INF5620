import sympy as sym
import numpy as np
V, t, I, w, dt = sym.symbols('V t I w dt')  # global symbols
f = None  # global variable for the source term in the ODE


def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u, t, t) + w**2*u


def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = DtDt(u, dt) + u*w**2 - ode_source_term(u)
    return sym.simplify(R)


def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted.
    Careful to evaluate everything around t=0"""
    f = ode_source_term(u)      # source term
    V = sym.diff(u, t)
    V = V.subs(t, 0)
    # V is the initial condition on first derivative and replaces u(0-dt)
    numer = 2*(u.subs(t, dt) - u.subs(t, 0) - dt*V)
    # numer is the numerator in the discrete second derivative,
    R = numer/(dt**2) + u.subs(t, 0)*w**2 - f.subs(t, 0)
    return sym.simplify(R)


def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u.subs(t, t + dt) - 2*u + u.subs(t, t - dt))/(dt**2)


def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u
    print "Initial conditions u(0)=%s, u'(0)=%s:" % \
          (u.subs(t, 0), sym.diff(u, t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    # f = sym.simplify(ode_lhs(u))
    f = sym.simplify(ode_source_term(u))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)


def my_linear():
    """
    Convenient symbolic definition of straight line
    """
    c, d = sym.symbols('c d')
    return c*t + d


def my_quadratic():
    """
    Convenient symbolic definition of quadratic function
    """
    c, d, e = sym.symbols('c d e')
    return c*t**2 + d*t + e


def my_cubic():
    """
    Convenient symbolic definition of cubic function
    """
    a0, a1, a2, a3 = sym.symbols('a0 a1 a2 a3')
    return a3*t**3 + a2*t**2 + a1*t + a0


if __name__ == '__main__':
    a = my_linear()
    main(a)
    b = my_quadratic()
    main(b)
    c = my_cubic()              # to test polynomial of degree 3
    main(c)

    import nose
    print '==== Running nose test on solver ===='
    # nose.main('vib_undamped_verify_mms')


def solver(u0, up0, omega, deltaT, T, f):
    """Solver for the problem u'' + u*omega**2 = f where u and f are
    functions of time t. Different symbols are used to avoid confusion
    with preceding sympy code.
    u0     = I  = Initial u value
    up0    = V  = u' at t=0
    omega  = w  = frequency
    deltaT = dt = size of time step
         T = End time of simulation
         f = f(t) = forcing function depenent on time in numerical form
    """
    if T < deltaT:              # short check of input variable
        raise ValueError("End time T can't be less than time step size deltaT")
    if T < 0:
        raise ValueError("End time T can't be less than 0")
    if T == 0:
        return u0

    deltaT = float(deltaT)      # to avoid integer division problems
    Nt = int(round(T/deltaT))   # number of time steps
    u = np.zeros(Nt + 1)        # empty vector for u values
    t = np.linspace(0, Nt*dt, Nt+1)  # time mesh

    u[0] = u0                   # first value

    # u[1] = 0.5*(u[0]*(2 - (deltaT**2)*(omega**2)) + 2*deltaT*up0 +
    #             f(0)*deltaT**2)  # first iteration
    u[1] = u[0]*(1 - 0.5*(deltaT*omega)**2) + deltaT*up0 + 0.5*deltaT**2*f(0)
    for n in xrange(1, Nt):     # remaining iterations
        u[n+1] = u[n]*(2 - (deltaT*omega)**2) - u[n-1] + f(t[n])*deltaT**2
    return u, t


def test_solver(tol=1e-14):
    """
    Function to test solver with source term from linear equation
    """
    assert 1 == 1.0
    assert 1 == 2
