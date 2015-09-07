import numpy as np
# import scitools.std as plt
import matplotlib.pyplot as plt


def simulate(beta=0.9, Theta=30, epsilon=0, num_periods=6,
             time_steps_per_period=60, plot=True):
    """
    Function to simulate motion of an elastic pendulum. Inputs are:
    beta - dimensionless paramater related to initial length,
mass, gravity and spring constant
    Theta - initial angle in degrees
    epsilon - initial stretch of wire
    num_periods - number of periods to simulate
    time_steps_per_period - time step resolution
    plots - if == True plots are made and shown
    """
    Theta = Theta*np.pi/180     # convert theta to radians
    # initialize vectors for x, y, theta, t
    P = 2*np.pi                 # period
    dt = P/time_steps_per_period  # P is always float because of pi
    T = num_periods*P
    t = np.arange(0, T, dt)  # T+dt because end is exclusive
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    theta = np.zeros_like(t)

    # initial condition for x, y and theta
    x[0] = (1 + epsilon)*np.sin(Theta)
    y[0] = 1 - (1 + epsilon)*np.cos(Theta)
    theta[0] = Theta

    # define convenient functions and constant
    gamma = float(beta)/(1 - beta)  # constant
    omega = 2*np.pi/P

    def length(x, y):
        """Convenient length measure"""
        return float(np.sqrt(x**2 + (y - 1)**2))

    def angle(x, y):
        """Convenient angle, not used much"""
        return np.arctan2(x, -y)
        # return np.arctan(float(x)/y)

    # calculat at time step 1
    x[1] = x[0] + 0.5*gamma*(dt**2)*(1 - beta/length(x[0], y[0]))*x[0]
    y[1] = y[0] + 0.5*gamma*(dt**2)*((1 - beta/length(x[0], y[0]))*(y[0] - 1)
                                     - beta)
    theta[1] = angle(x[1], y[1])
    # solve for remaining time steps
    for n in xrange(1, len(t)-1):
        x[n+1] = 2*x[n] - x[n-1] + (dt**2)*gamma*(1 -
                                                  (beta/length(x[n],
                                                               y[n])))*x[n]
        y[n+1] = 2*y[n] - y[n-1] + (dt**2)*(gamma *
                                            (1 - (beta/length(x[n], y[n])))
                                            * (y[n] - 1) - beta)
        theta[n] = angle(x[n], y[n])
    # start plotting
    if plot:
        plt.figure(0)
        plt.plot(x, y, 'b-')
        plt.title('Pendulum motion')
        plt.xlim(x.min(), x.max())
        plt.ylim(1.3*y.min(), 1)
        plt.axis('equal')
        plt.savefig('tmp_xy.png')
        plt.savefig('tmp_xy.pdf')

        # Plot theta in degrees
        plt.figure(1)
        plt.plot(t, theta*180/np.pi, 'b-')
        plt.title('Angular displacement in degrees')
        plt.ylabel('Degrees')
        plt.xlabel('Time')
        plt.savefig('tmp_theta.png')
        plt.savefig('tmp_theta.pdf')
        if abs(Theta) < 10*np.pi/180:
            # Compare theta and theta_e for small angles (<10 degrees)
            theta_e = Theta*np.cos(omega*t)  # non-elastic scaled sol.
            plt.figure(2)
            plt.plot(t, theta, label='Theta elastic')
            plt.plot(t, theta_e, label='Theta non-elastic')
            plt.title('Elastic vs non-elastic pendulum for $\\beta=$%g' % beta)
            plt.ylabel('Degrees')
            plt.xlabel('Time')
            plt.savefig('tmp_compare.png')
            plt.savefig('tmp_compare.pdf')
        plt.show()            # for testing
        # Plot y vs x (the real physical motion)
        return x, y, theta, t


def test_equiliberium():
    """Test that starting from rest makes x=y=theta=0"""
    x, y, theta, t = simulate(beta=0.9, Theta=0, epsilon=0, num_periods=6,
                              time_steps_per_period=10, plot=False)
    tol = 1e-14
    assert np.abs(x.max()) < tol
    assert np.abs(y.max()) < tol
    assert np.abs(theta.max()) < tol


def test_vertical_motion():
    beta = 0.9
    omega = np.sqrt(beta/(1-beta))
    # Find num_periods. Recall that P=2*pi for scaled pendulum
    # oscillations, while here we don't have gravity driven
    # oscillations, but elastic oscillations with frequency omega.
    period = 2*np.pi/omega
    # We want T = N*period
    N = 5
    # simulate function has T = 2*pi*num_periods
    num_periods = 5/omega
    n = 600
    time_steps_per_period = omega*n
    y_exact = lambda t: -0.1*np.cos(omega*t)
    x, y, theta, t = simulate(beta, 0, 0.1,
                              num_periods,
                              time_steps_per_period,
                              False)
    tol = 0.00055              # ok tolerance for the above resolution
    # No motion in x direction is epxected
    assert np.abs(x.max()) < tol
    # Check motion in y direction
    y_e = y_exact(t)
    diff = np.abs(y_e - y).max()
    if diff > tol:              # plot
        plt.plot(t, y, label='y calculated')
        plt.plot(t, y_e, label='exact')
        plt.legend(loc='best')
        plt.show()
        raw_input('Error in test_vertical_motion; type CR:')
        assert diff < tol, 'diff=%g' % diff


def demo(beta=0.99, Theta=40, num_periods=3):
    """
    Function to demonstrate simulator
    """
    x, y, theta, t = simulate(beta=beta, Theta=Theta, epsilon=0,
                              num_periods=num_periods,
                              time_steps_per_period=600, plot=True)


def text_demo():
    """
    Demonstrate to give plots from the text.
    """
    demo(beta=0.999, Theta=40, num_periods=3)
    demo(beta=0.93, Theta=40, num_periods=1)
