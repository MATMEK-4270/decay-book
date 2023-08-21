import numpy as np

def solver(I, a, T, dt, theta):
    """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
    dt = float(dt)            # avoid integer division
    Nt = int(round(T/dt))     # no of time intervals
    T = Nt*dt                 # adjust T to fit time step dt
    u = np.zeros(Nt+1)           # array of u[n] values
    t = np.linspace(0, T, Nt+1)  # time mesh
    u[0] = I                  # assign initial condition
    for n in range(0, Nt):    # n=0,1,...,Nt-1
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
    return u, t

def u_exact(t, I, a):
    return I*np.exp(-a*t)

import matplotlib.pyplot as plt

def plot_numerical_and_exact(theta, I, a, T, dt):
    """Compare the numerical and exact solution in a plot."""
    u, t = solver(I=I, a=a, T=T, dt=dt, theta=theta)
    t_e = np.linspace(0, T, 1001)        # fine mesh for u_e
    u_e = u_exact(t_e, I, a)
    plt.plot(t, u, 'r--o', t_e, u_e, 'b-')
    plt.legend(['numerical', 'exact'])
    plt.xlabel('t'); plt.ylabel('u')
    plt.title('theta=%g, dt=%g' % (theta, dt))
    plt.savefig('plot_%s_%g.png' % (theta, dt))
plot_numerical_and_exact(I=1, a=2, T=8, dt=0.8, theta=1)
plt.show()
