'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 3
Problem 1b
'''
import numpy as np
import matplotlib.pyplot as plt
import math


def pendulum_sim(theta0, momen0, steps, stepSize): 
    thetaLog=[theta0] 
    momenLog=[momen0]

    #Forward Euler step
    theta = theta0 + stepSize*momen0
    momen = momen0 - stepSize*math.sin(theta0)
    thetaLog.append(theta)
    momenLog.append(momen)

    #LeapFrog Loop
    for i in range(steps-2):
        theta = thetaLog[-2] + 2.0*stepSize*momenLog[-1]
        momen = momenLog[-2] - 2.0*stepSize*math.sin(thetaLog[-1])
        thetaLog.append(theta)
        momenLog.append(momen)
    return thetaLog, momenLog
    
    
    
    

stepSize=1.0/1000.0
steps=100*100
resolution=10
for i in range(resolution):
    for j in range(resolution): 
        theta0=-3.0 + i*6.0/float(resolution)
        momen0=-3.0 + j*6.0/float(resolution)
        thetaLog, momenLog = pendulum_sim(theta0, momen0, steps, stepSize)
        thetaLog=((np.array(thetaLog)+math.pi)%(2*math.pi))-math.pi
        plt.plot(thetaLog,momenLog,'.b')
plt.title('Phase Space Plot')
plt.ylabel('Momentum')
plt.xlabel('Theta (radians)')
plt.show()
#plt.plot(range(steps),thetaLog)
#plt.title('Theta')
#plt.show()
#plt.plot(range(steps),momenLog)
#plt.title('Momentum')
#plt.show()
