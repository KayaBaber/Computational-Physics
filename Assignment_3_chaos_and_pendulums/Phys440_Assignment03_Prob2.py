'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 3
Problem 2
'''
import numpy as np
import matplotlib.pyplot as plt
import math


def forced_damped_pendulum_sim(theta0, thetaDot0, b, omega, gamma, steps, stepSize): 
    thetaLog=[theta0] 
    thetaDotLog=[thetaDot0]
    timeLog=[0.0]

    #Forward Euler step
    thetaDouble0=-b*thetaDot0 - math.sin(theta0) + gamma*math.cos(omega*stepSize)
    theta = theta0 + stepSize*thetaDot0
    thetaDot = thetaDot0 + stepSize*thetaDouble0
    
    thetaLog.append(theta)
    thetaDotLog.append(thetaDot)

    #LeapFrog Loop
    for i in range(1,steps-1):
        theta = thetaLog[-2] + 2.0*stepSize*thetaDotLog[-1]
        thetaDouble=-b*thetaDot - math.sin(theta) + gamma*math.cos(omega*i*stepSize)
        thetaDot = thetaDotLog[-2] + 2.0*stepSize*thetaDouble
        
        timeLog.append(i*stepSize)
        thetaLog.append(theta)
        thetaDotLog.append(thetaDot)
    timeLog.append(steps*stepSize)    
    return thetaLog, thetaDotLog, timeLog





theta0=-0.1
thetaDot0=0.1
b=0.0
omega=0.7
gamma=0.0
steps=100*100
stepSize=1.0/1000.0


thetaLog, thetaDotLog, timeLog = forced_damped_pendulum_sim(theta0, thetaDot0, b, omega, gamma, steps, stepSize)
plt.plot(timeLog,thetaLog)
plt.show()

