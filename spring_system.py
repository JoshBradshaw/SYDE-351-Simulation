import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

# parameters, set through measurement and making assumptions
U1 = 1.
I2 = 1.
C4 = 1.
TF = 1.
I7 = 1.
R9 = 20
C11 = 1.
R12 = 1.

# slip is given by X7 - X9

def f(y, t):
    """the system of differential equaitons derived from the bond graph model"""
    P2i = y[0]
    Q4i = y[1]
    P7i = y[2]
    Q11i = y[3]
    X2i = y[4]

    P2_ = U1 - Q4i/C4
    Q4_ = P2i/I2 - P7i/(TF*I7)
    P7_ = Q4i/(TF*C4) - R9 *(P7i/I7 + Q11i/(R12*C11))
    Q11_ = P7i/I7 - Q11i/(R12*C11)
    X2_ = P2i/I2
    X7_ = P7i/I7
    X9_ = P7i/I7 + Q11i/(R12*C11)
    return [P2_, Q4_, P7_, Q11_, X2_, X7_, X9_]

def u_k(vr, vt):
    """4th order model of the bushing friction"""
    A = 1.
    B1 = 2.
    B2 = 2.
    return math.tanh(vr/vt) + (B1 * (vr/vt))/(1 + B2*(vr/vt)**4)


# Initial conditions
P2 = 0. # initial velocity
Q4 = 0. # initial force
P7 = 0. # initial force
Q11 = 0. # initial velocity
X2 = 0.
X7 = 0.
X9 = 0.
y0 = [P2, Q4, P7, Q11, X2, X7, X9]
t = np.linspace(0, 100, 2000)

# solve the DEs
soln = odeint(f, y0, t)
P2 = soln[:, 0]
Q4 = soln[:, 1]
P7 = soln[:, 2]
Q11 = soln[:, 3]
X2 = soln[:, 4]
X7 = soln[:, 5]
X9 = soln[:, 6]

plt.figure()
plt.plot(t, X2, 'r--')

plt.xlabel('Time in seconds')
plt.ylabel('kg m/s')
plt.title('Position of the Hanging Mass during initial oscillation')
plt.legend(loc=0)
plt.show()