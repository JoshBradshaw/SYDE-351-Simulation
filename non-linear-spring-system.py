from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

g = 9.8 # m/s^2
M2 = 0.3 # kg
K4 = 150 #N/m
Rdisk = 0.014 # m
J7 = 0.00002
J12 = 0.00002
c9 = 1
c11 = 10
K14 = 20 # n/m

# slip is given by X7 - X9
force_values = []

def f(y, t):
    """the system of differential equaitons derived from the bond graph model"""
    P2i = y[0]
    Q4i = y[1]
    P7i = y[2]
    P12i = y[3]
    Q14i = y[4]
    X2 = y[5]
    X7 = y[6]
    X10 = y[7]

    wr = P7i/J7
    wt = P12i/J12

    P2_ = M2*g - nl_k4(Q4i)
    Q4_ = P2i/M2 - Rdisk*P7i/J7
    P7_ = nl_k4(Q4i) - u_k(wr, wt)
    P12_ = u_k(wr, wt) - Rdisk*nl_k14(Q14i) - P12i*c11/J12
    Q14_ = Rdisk*P12i/J12
    X5_ = P2_/M2
    return [P2_, Q4_, P7_, P12_, Q14_, X5_, wr, wt]

def u_k(wr, wt):
    """4th order model of the bushing friction"""
    A = 3.65 # 3.7
    B1 = 0.8
    B2 = 0.5
    return A*(math.tanh(wr/wt) + (B1 * (wr/wt))/(1 + B2*(wr/wt)**4))

def nl_k4(q):
    if q >= 0:
        return K4 * q
    else:
        return 0

def nl_k14(q):
    if q >= 0:
        return K14 * q
    else:
        return 0

# Initial conditions
P2 = 1. # initial velocity
Q4 = 0. # initial force
P7 = 0.00000001 # initial force
P12 = 0.000000001
Q14 = 0. # initial velocity
X2 = 0.
uk = 0.
x7 = 0.
x9 = 0.
y0 = [P2, Q4, P7, P12, Q14, X2, x7, x9]
t = np.linspace(0, 10, 4000)

# solve the DEs
soln = odeint(f, y0, t)
P2 = soln[:, 0]
Q4 = soln[:, 1]
P7 = soln[:, 2]
P12 = soln[:, 3]
Q13 = soln[:, 4]
X2 = soln[:, 5]
x7 = soln[:, 6]
x9 = soln[:, 7]

plt.figure()
plt.plot(t, X2, 'r--')

plt.figure()
plt.plot(t, x7, 'r--', t, x9, 'b--')

plt.xlabel('Time in seconds')
plt.ylabel('m')
plt.title('Position of the Hanging Mass')
plt.legend(loc=0)
plt.show()