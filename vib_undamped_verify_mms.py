import sympy as sym
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
    numer = 2*u.subs(t, 0+dt)*(1 - dt*V) - 2*u.subs(t, 0)
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


def linear():
    main(lambda t: V*t + I)


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
    

# def quadratic():
#     main(lambda t: b*t**2 + c*t + d)

if __name__ == '__main__':
    a = my_linear()
    main(a)
    b = my_quadratic()
    main(b)
    c = my_cubic()
    main(c)
