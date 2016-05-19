'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 5 - PDEs
Exploration 2 - Parabolic PDEs: The Wave Equation
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math

def make_quad_ops(N,stepSize):
    iden = np.identity(N)
    #make derivative operator matrix
    #i(diag (0,1,2,3,4,.. -4, -3, -2, -1)), where i is imaginary i
    opFlat = np.append(np.linspace(0,N/2, (N/2)+1)*1J,
                       np.linspace((-N/2)+1,0,(N/2))[:-1]*1J)
    opDiag = np.diag(opFlat*(stepSize/2.))
    top = np.append(iden, opDiag, axis=1)
    bot = np.append(opDiag, iden, axis=1)
    posQuad = np.append(top, bot, axis=0)
    negTop = np.append(iden, -opDiag, axis=1)
    negBot = np.append(-opDiag, iden, axis=1)
    negQuad = np.append(negTop, negBot, axis=0)
    return posQuad, negQuad


def step_forward(velDen,posQuad,negQuad):
    #matrix multiply the negative quad by the FFT column vector
    RHS = np.dot(negQuad,velDen)
    #linear algebra solve the pos_quad*newFFTvector = above-result for newFFTvector
    solution = LA.tensorsolve(posQuad,RHS)
    return solution



L = 2.*math.pi  #set the x range to (0->2pi)
N = 100        #number of spatial intervals and points (since it loops)
steps = 629    #number of timesteps
stepSize = 0.01  #temporal step size

#make initial velocity array in physical space
velPhysFlat = [0 + 0J]*N
velPhysGauss = np.exp(-10 * ( (np.linspace(0,L,N+1)[:-1]-math.pi) ** 2)) + [0J]*N
velPhysGaussL = np.exp(-10 * ( (np.linspace(0,L,N+1)[:-1]- 1.5*math.pi) ** 2)) + [0J]*N
#velPhysSinGauss = np.multiply(np.sin((0.2*np.linspace(0,L,N+1)[:-1])*4),np.exp(-10*( (np.linspace(0,L,N+1)[:-1]-*math.pi) ** 2)))

#make initial density array in physical space
denPhysFlat = [0 + 0J]*N
denPhysGauss = np.exp(-10 * ( (np.linspace(0,L,N+1)[:-1]-math.pi) ** 2)) + [0J]*N
denPhysGaussR = np.exp(-10 * ( (np.linspace(0,L,N+1)[:-1]- 0.5*math.pi) ** 2)) + [0J]*N
denPhysSinGauss = np.multiply(np.sin((np.linspace(0,L,N+1)[:-1])*4),np.exp(-1*((np.linspace(0,L,N+1)[:-1]-math.pi) ** 2)))

#fft both to fourier space
velF = np.fft.fft(denPhysFlat)
denF = np.fft.fft(denPhysSinGauss)

#make a column vector of density_f appended to velocity_f
velDen = velF
velDen = np.append(velDen,denF)

#make quad matrix [[I][op],[op][I]] and negative quad [[I][-op],[-op][I]]
posQuad, negQuad = make_quad_ops(N,stepSize)

#simulation loop
velLog = []
denLog = []
for i in range(steps):
    #make seperate copies of the velocity_f and density_f components of newFFTvector
    velF = velDen[:N]
    denF = velDen[N:]
    #inverse FFT the components and append the real parts of velocity and density to a log
    vel = np.fft.ifft(velF)
    den = np.fft.ifft(denF)
    velLog.append(vel.real)
    
    plt.plot(np.linspace(0,L,N),vel.real)
#    plt.ylim([-1,1])
#    plt.xlim([0,L])
#    plt.ylabel("Velocity")
#    plt.xlabel("Position")    
#    plt.savefig('frames/velocityImage'+str(i+1)+'.png',bbox_inches='tight')    
#    plt.clf()
    
    plt.plot(np.linspace(0,L,N),den.real)
    plt.plot(np.linspace(0,L,N),den.real-vel.real)
    plt.ylim([-1,1])
    plt.xlim([0,L])
    plt.ylabel("Density, Velocity, Difference")
    plt.xlabel("Position")     
    plt.savefig('frames/comboImage'+str(i+1)+'.png',bbox_inches='tight')    
    plt.clf()
    print i
    denLog.append(den.real)
    #step forward
    velDen = step_forward(velDen,posQuad,negQuad)



#plot the velocity and density logs in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.linspace(0,L,N)
y = np.linspace(0,steps*stepSize,steps)
XX, YY = np.meshgrid(x,y)
ZZ = velLog
ax.plot_surface(XX, YY, ZZ, 
                #cmap=cm.spectral
               )
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Velocity')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.linspace(0,L,N)
y = np.linspace(0,steps*stepSize,steps)
XX, YY = np.meshgrid(x,y)
ZZ = denLog
ax.plot_surface(XX, YY, ZZ, 
                #cmap=cm.spectral
               )
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Density')
plt.show()






 
