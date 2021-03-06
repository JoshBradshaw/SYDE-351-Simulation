from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from pprint import pprint
import math

g = 9.8 # m/s^2
M2 = 0.3 # kg
K4 = 136 #N/m
Rdisk = 0.014 # m
J7 = 0.00002
J12 = 0.00002
c9 = 10
c11 = 10
K14 = 15 # n/m

# slip is given by X7 - X9
friction = []

def f(y, t):
    """the system of differential equaitons derived from the bond graph model"""
    P2i = y[0]
    Q4i = y[1]
    P7i = y[2]
    P12i = y[3]
    Q14i = y[4]
    X2 = y[5]

    wr = P7i/J7
    wt = P12i/J12
    #print "u_k(wr, wt) {}".format(u_k(wr, wt))

    P2_ = M2*g - Q4i*K4
    Q4_ = P2i/M2 - Rdisk*P7i/J7
    P7_ = Q4i*K4 - u_k(wr, wt)
    P12_ = u_k(wr, wt) - Rdisk*Q14i*K14 - P12i*c11/J12
    Q14_ = Rdisk*P12i/J12
    X5_ = P2_/M2
    print "c9*(wr-wt): {} u_k: {}".format(c9*(wr-wt), u_k(wr, wt))
    return [P2_, Q4_, P7_, P12_, Q14_, X5_, wr, wt]

def u_k(wr, wt):
    """4th order model of the bushing friction"""
    A = 0.364
    B1 = 0.47
    B2 = 0.5
    return math.copysign(A*(math.tanh(wr/wt) + (B1 * (wr/wt))/(1 + B2*(wr/wt)**4)), wr-wt)


# Initial conditions
P2 = 1. # initial velocity
Q4 = 0. # initial force
P7 = 0.00000001 # initial force
P12 = 0.000000001
Q14 = 0. # initial velocity
X2 = 0.
wr = 0.
wt = 0.
y0 = [P2, Q4, P7, P12, Q14, X2, wr, wt]
t = np.linspace(0, 10, 4000)

# solve the DEs
soln = odeint(f, y0, t)
P2 = soln[:, 0]
Q4 = soln[:, 1]
P7 = soln[:, 2]
P12 = soln[:, 3]
Q13 = soln[:, 4]
X2 = soln[:, 5]
wr = soln[:, 6]
wt = soln[:, 7]

plt.figure()
plt.plot(t, wr, 'r--', t, wt, 'b--', t, wr-wt, 'g--')
plt.figure()
plt.plot(t, X2, 'r--')

plt.xlabel('Time in seconds')
plt.ylabel('cm')
plt.title('Position of the Hanging Mass')
plt.legend(loc=0)
plt.show()