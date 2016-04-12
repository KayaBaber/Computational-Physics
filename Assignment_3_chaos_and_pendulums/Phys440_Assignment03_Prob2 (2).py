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
    #pendulum driven-damped function
    theta=thetas[0]
    thetaDot=thetas[1]
    thetaDouble=-b*thetaDot - math.sin(theta) + gamma*math.cos(omega*t)
    return thetaDot, thetaDouble
    

#initial conditions
theta0=-0.0
thetaDot0=0.0
thetas=[theta0,thetaDot0]

#generating loop
for i in range(7):
    #constants
    b=0.05
    omega=0.7
    gamma=0.4+(i*0.1)
    #FIX YO OMEGA STUFF

    #computation parameters
    steps=100
    periods=10
    t = np.linspace(0, periods*(math.pi*2.0/omega), steps*(math.pi*2.0/omega)+1)

    #ODE solution
    sol = odeint(f, thetas, t, args=(b, gamma, omega))

    #TAKE THE STROBE

    plt.plot(t, sol[:, 1], 'b', label='thetaDot(t)')
    plt.xlabel('time')
    plt.ylabel('theta')
    plt.grid()
    #plt.savefig('/Users/student/kbaber/Desktop/Phys440/Assignment 3/plots//gamma'+str(gamma)+'_thetaDot_t.png',bbox_inches='tight')
    plt.savefig('\Users\Kaya\Google Drive\School\Phys 440\Assignments\Assignment 3\plots\\gamma'+str(gamma)+'_thetaDot_t.png',bbox_inches='tight')
    #plt.show()
    plt.clf()

    plt.plot(sol[:,0], sol[:, 1], 'g', label='theta-Dot(theta)')
    plt.xlabel('theta')
    plt.ylabel('theta-Dot')
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.savefig('/Users/student/kbaber/Desktop/Phys440/Assignment 3/plots//gamma'+str(gamma)+'_thetaDot_theta.png',bbox_inches='tight')
    plt.savefig('\Users\Kaya\Google Drive\School\Phys 440\Assignments\Assignment 3\plots\\gamma'+str(gamma)+'_thetaDot_theta.png',bbox_inches='tight')
    #plt.show()
    plt.clf()


print t