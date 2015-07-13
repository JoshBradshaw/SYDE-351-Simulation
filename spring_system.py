import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# parameters, set through measurement and making assumptions
U1 = 1 # mg
L2 = 1. # mass
C3 = 3. # spring const
R4 = 0.95 # friction relation
C5 = 5. # sprint const
L6 = 1. # mass

def f(y, t):
	P2i = y[0]
	Q3i = y[1]
	Q5i = y[2]
	P6i = y[3]

	f0 = U1 - Q3i/C3 - R4*P2i/L2 + R4*P6i/L6
	f1 = -Q3i/C3 + R4*P2i/L2 - R4*P6i/L6 - Q5/C5
	f2 = Q3i/L2 - P6i/L6
	f3 = P6i/L6
	return [f0, f1, f2, f3]

# Initial conditions
P2 = 0. # initial velocity
Q3 = 0. # initial force
Q5 = 0. # initial force
P6 = 0. # initial velocity
y0 = [P2, Q3, Q5, P6]
t = np.linspace(0, 20, 2000)

# solve the DEs
soln = odeint(f, y0, t)
I2 = soln[:, 0]
U2 = soln[:, 1]
U5 = soln[:, 2]
I6 = soln[:, 3]

plt.figure()
plt.plot(t, U5, label='Inertia 2')
plt.xlabel('Time in seconds')
plt.ylabel('kg m/s')
plt.title('Momentum of Hanging Mass during initial oscillation')
plt.legend(loc=0)
plt.show()