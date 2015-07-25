from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

# parameters, set through measurement and making assumptions
SE_1 = 2.94
I_2 = 0.3
C_4 = 136.0
TF_5 = 4092.0
I_7 = 0.00002405 # revolution
R_9 = 0.001 # no slip for now
R_11 = 1.0 # 
TF_12 = 0.00024
C_13 = 15.0

# slip is given by X7 - X9

def f(y, t):
    """the system of differential equaitons derived from the bond graph model"""
    P2i = y[0]
    Q4i = y[1]
    P7i = y[2]
    Q13i = y[3]

    v11 = (1.0/(1.0 + R_9/R_11)) * ((R_9*P7i/(R_11*I_7)) - TF_12*Q13i/(R_11*C_13))
    P2_ = SE_1 - Q4i/C_4
    Q4_ = P2i/I_2 - P7i/(TF_5 * I_7)
    P7_ = Q4i/(C_4*TF_5) - R_9*(P7i/I_7 - v11)
    Q13_ = v11/TF_12

    return [P2_, Q4_, P7_, Q13_]

def u_k(vr, vt):
    """4th order model of the bushing friction"""
    A = 1.
    B1 = 2.
    B2 = 2.
    return math.tanh(vr/vt) + (B1 * (vr/vt))/(1 + B2*(vr/vt)**4)


# Initial conditions
P2 = 0 # initial velocity
Q4 = 0. # initial force
P7 = 0. # initial force
Q13 = 0. # initial velocity
y0 = [P2, Q4, P7, Q13]
t = np.linspace(0, 1000, 2000)

# solve the DEs
soln = odeint(f, y0, t)
P2 = soln[:, 0]
Q4 = soln[:, 1]
P7 = soln[:, 2]
Q13 = soln[:, 3]

plt.figure()
plt.plot(t, Q13, 'r--')

plt.xlabel('Time in seconds')
plt.ylabel('kg m/s')
plt.title('Position of the Hanging Mass')
plt.legend(loc=0)
plt.show()