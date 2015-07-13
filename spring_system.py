import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# parameters, set through measurement and making assumptions
U1 = 1 # mg
L2 = 1. # mass
C5 = 3. # spring const
R6 = 1.2 # friction relation
C8 = 5. # spring const
L9 = 1. # mass

def f(y, t):
	P2i = y[0]
	Q5i = y[1]
	Q8i = y[2]
	P9i = y[3]

	f0 = U1 - Q5i/C5 - R6*P2i/C5 + R6*P9i/L9
	f1 = P2i/C5 - P9i/L9
	f2 = P9i/L9
	f3 = -Q5i/C5 - R6*P2i/L2 + R6*P9/L9 - Q8i/C8
	return [f0, f1, f2, f3]

# Initial conditions
P2 = 0. # initial velocity
Q5 = 0. # initial force
Q8 = 0. # initial force
P9 = 0. # initial velocity
y0 = [P2, Q5, Q8, P9]
t = np.linspace(0, 100, 2000)

# solve the DEs
soln = odeint(f, y0, t)
P2 = soln[:, 0]
Q2 = soln[:, 1]
Q8 = soln[:, 2]
P9 = soln[:, 3]

plt.figure()
plt.plot(t, -Q2, label='Inertia 2')
plt.xlabel('Time in seconds')
plt.ylabel('kg m/s')
plt.title('Momentum of Hanging Mass during initial oscillation')
plt.legend(loc=0)
plt.show()