'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 3
Problem 2
'''
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math

def f(thetas, t, b, gamma, omega):
    theta=thetas[0]
    thetaDot=thetas[1]
    thetaDouble=-b*thetaDot - math.sin(theta) + gamma*math.cos(omega*t)
    return thetaDot, thetaDouble
    

theta0=-0.1
thetaDot0=0.1
thetas=[theta0,thetaDot0]
b=0.05
omega=0.7
gamma=0.4
steps=101
boost=5
t = np.linspace(0, 10*boost, steps*boost)
sol = odeint(f, thetas, t, args=(b, gamma, omega))
plt.plot(t, sol[:, 0], 'b', label='theta(t)')
#plt.plot(t, sol[:, 1], 'g', label='thetaDot(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('theta')
plt.grid()
plt.show()


plt.plot(sol[:,0], sol[:, 1], 'g', label='theta-Dot(theta)')
plt.legend(loc='best')
plt.xlabel('theta')
plt.ylabel('theta-Dot')
plt.grid()
plt.show()

#plt.plot(sol[:,0], sol[:, 1], 'g', label='theta-Dot(theta)')
#plt.legend(loc='best')
#plt.xlabel('theta')
#plt.ylabel('theta-Dot')
#plt.grid()
#plt.show()
