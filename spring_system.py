from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

g = 9.8 # m/s^2
M2 = 0.3 # kg
K4 = 136 #N/m
Rdisk = 0.014 # m
J7 = 0.00002
c9 = 0.02
c11 = 0.002
K13 = 15 # n/m


# slip is given by X7 - X9
actual_ratios = []


def f(y, t):
    """the system of differential equaitons derived from the bond graph model"""
    P2i = y[0]
    Q4i = y[1]
    P7i = y[2]
    Q13i = y[3]
    X1i = y[4]
    W8i = y[5]
    W11i = y[6]


    w11_ = (c9*P7i/J7 - Rdisk*Q13i*K13)/(c11 + c9)
    P2_ = M2*g - Q4i*K4
    Q4_ = P2i/M2 - Rdisk*P7i/J7
    P7_ = Rdisk*Q4i*K4 - c9*(P7i/J7 - w11_)
    Q13_ = Rdisk*w11_
    X1_ = P2_/M2
    W8_ = P7i/J7
    print "c9*P7i/J7 {} c9*(P7i/J7 - w11_) {}".format(c9*P7i/J7, c9*(P7i/J7 - w11_))
    return [P2_, Q4_, P7_, Q13_, X1_, W8_, w11_]

def u_k(wr_wt):
    """4th order model of the bushing friction"""
    A = 0.02
    B1 = 0.47
    B2 = 8.
    return A*(math.tanh(wr_wt) + (B1 * (wr_wt))/(1 + B2*(wr_wt)**4))



# Initial conditions
P2 = 1. # initial velocity
Q4 = 0. # initial force
P7 = 0. # initial force
Q13 = 0. # initial velocity
X1 = 0.
W8 = 0.
W11 = 0.
y0 = [P2, Q4, P7, Q13, X1, W8, W11]
t = np.linspace(0, 10, 4000)

# solve the DEs
soln = odeint(f, y0, t)
P2 = soln[:, 0]
Q4 = soln[:, 1]
P7 = soln[:, 2]
Q13 = soln[:, 3]
X1 = soln[:, 4]
W8 = soln[:, 5]
W11 = soln[:, 6]

plt.figure()
actual_ratios = np.linspace(-50, 50, 2000)
uk_vals = [u_k(ratio) for ratio in actual_ratios]

plt.plot(actual_ratios, uk_vals, 'r--')
#plt.show()

plt.figure()
plt.plot(t, W8, 'r--', t, W11, 'b--')

plt.xlabel('Time in seconds')
plt.ylabel('cm')
plt.title('Position of the Hanging Mass')
plt.legend(loc=0)
plt.show()