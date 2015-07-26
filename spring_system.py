from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

# parameters, set through measurement and making assumptions
SE_1 = 2.94 # N
I_2 = 0.3 # kg
C_4 = 1./136. # m/N
TF_5 = 71 # rad/m
I_7 = 0.002 # revolution
R_9 = 1.6 # no slip for now
R_11 = 1.6 # 
TF_12 = 0.014 # m/rad
C_13 = 1/15. # m/N

# slip is given by X7 - X9

def f(y, t):
    """the system of differential equaitons derived from the bond graph model"""
    P2i = y[0]
    Q4i = y[1]
    P7i = y[2]
    Q13i = y[3]
    X1_ = y[4]

    v11 = (R_9*P7i/I_7 - TF_12*Q13i/C_13)/(R_11+R_9)
    P2_ = SE_1 - Q4i/C_4
    Q4_ = P2i/I_2 - P7i/(TF_5 * I_7)
    P7_ = Q4i/(C_4*TF_5) - R_9*(P7i/I_7 - v11)
    Q13_ = v11*TF_12
    X1_ = P2_/I_2

    return [P2_, Q4_, P7_, Q13_, X1_]

def u_k(vr, vt):
    """4th order model of the bushing friction"""
    A = 1.
    B1 = 2.
    B2 = 2.
    return math.tanh(vr/vt) + (B1 * (vr/vt))/(1 + B2*(vr/vt)**4)


# Initial conditions
P2 = 3. # initial velocity
Q4 = 0. # initial force
P7 = 0. # initial force
Q13 = 0. # initial velocity
X1 = 0.
y0 = [P2, Q4, P7, Q13, X1]
t = np.linspace(0, 10, 4000)

# solve the DEs
soln = odeint(f, y0, t)
P2 = soln[:, 0]
Q4 = soln[:, 1]
P7 = soln[:, 2]
Q13 = soln[:, 3]
X1 = soln[:, 4]

plt.figure()
plt.plot(t, X1, 'r--')

plt.xlabel('Time in seconds')
plt.ylabel('kg m/s')
plt.title('Position of the Hanging Mass')
plt.legend(loc=0)
plt.show()