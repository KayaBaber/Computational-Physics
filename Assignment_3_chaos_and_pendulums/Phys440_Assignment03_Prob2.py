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

#constants
b=0.05
omega=0.7

#computation parameters
steps=100
periods=310
t = np.linspace(0, periods*(math.pi*2.0*omega), steps*periods+1)

#generating loop
for i in range(7):

    gamma=0.4+(i*0.1)

    #ODE solution
    sol = odeint(f, thetas, t, args=(b, gamma, omega))

    #Cut off data from before 200 driving periods
  
    #plot theta vs time
    plt.plot(t[210*steps:], sol[:, 1][210*steps:], 'b', label='thetaDot(t)')
    plt.xlabel('time')
    plt.ylabel('theta-Dot')
    plt.grid()
    plt.savefig('plots/gamma'+str(gamma)+'_thetaDot_t.png',bbox_inches='tight')
    #plt.savefig('plots\\gamma'+str(gamma)+'_thetaDot_t.png',bbox_inches='tight')
    #plt.show()
    plt.clf()

    #clips the plot to keep theta between -pi and +pi
    thetaLog=((np.array(sol[:,0][210*steps:])+math.pi)%(2*math.pi))-math.pi
    #plot phase space plot
    plt.plot(thetaLog, sol[:, 1][210*steps:], 'g.', label='theta-Dot(theta)')
    plt.xlabel('theta')
    plt.ylabel('theta-Dot')
    plt.title('Phase Space Plot')
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('plots/gamma'+str(gamma)+'_thetaDot_theta.png',bbox_inches='tight')
    #plt.savefig('plots\\gamma'+str(gamma)+'_thetaDot_theta.png',bbox_inches='tight')
    #plt.show()
    plt.clf()
    
    
    
    #selects only points that coincide with the period omega
    strobedTheta=sol[:,0][210*steps:-1:steps]
    strobedThetaDot=sol[:,1][210*steps:-1:steps]
    strobedTheta=((strobedTheta+math.pi)%(2*math.pi))-math.pi
    #plot strobed phase space plot
    plt.plot(strobedTheta, strobedThetaDot, 'r.', label='theta-Dot(theta)')
    plt.xlabel('theta')
    plt.ylabel('theta-Dot')
    plt.title('Strobed Phase Space Plot')
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('plots/gamma'+str(gamma)+'_thetaDot_theta_strobed.png',bbox_inches='tight')
    #plt.savefig('plots\\gamma'+str(gamma)+'_thetaDot_theta.png',bbox_inches='tight')
    #plt.show()
    plt.clf()
    
    

