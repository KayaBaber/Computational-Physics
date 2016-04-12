import math
import numpy as np
import matplotlib.pyplot as plt


def pythag_pos_diff(objPosArray1,objPosArray2):
    
    #calculates the pythagorian position difference between two n-dimensional position arrays
    pythagSum=0
    for l in range(len(objPosArray1)):
        pythagSum+=math.pow(objPosArray1[l]-objPosArray2[l],2)
    return math.pow(pythagSum,1.0/2.0)


def nbody_sim_euler(posArray,velArray, massArray,stepSize,stepCount):
    '''Returns an array of computed position steps for the inputed objects
    '''
    numObjects=len(massArray)
    numDimensions=len(posArray[0])
    #simulation loop using forward Euler
    #1. update position
    #2. compute acceleration
    #3. update velocity
    G=math.pow(2.0*math.pi,2.0)
    historicPosArray=[[[]for j in range(numDimensions)] for j in range(numObjects)]
    historicVelArray=[[[]for j in range(numDimensions)] for j in range(numObjects)]
    
    #start with forward Euler
    for i in range(stepCount):
        #update positon
        for j in range(numObjects):
            for l in range(numDimensions):
                historicPosArray[j][l].append(posArray[j][l])
                posArray[j][l]+=stepSize*velArray[j][l]
        #compute acceleration
        accelArray=[[0]*numDimensions for j in range(numObjects)]
        for j in range(numObjects):
            for l in range(numDimensions):
                for k in range(j+1,numObjects):
                    #print posArray[j],posArray[k]
                    newAccel=G*massArray[k]*(posArray[k][l]-posArray[j][l]) / math.pow(pythag_pos_diff(posArray[j],posArray[k]),3.0)
                    oppositeNewAccel=-newAccel*massArray[j]/massArray[k]
                    accelArray[j][l]+=newAccel
                    accelArray[k][l]+=oppositeNewAccel
          
        #update velocity
        for j in range(numObjects):
            for l in range(numDimensions):
                velArray[j][l]+=stepSize*accelArray[j][l]
        #ang=angular_momentum(posArray, velArray, massArray)
        energy=total_energy(posArray, velArray, massArray)
        #print ang[0]
        print energy
    return historicPosArray
    

def angular_momentum(posArray, velArray, massArray):
    '''calculates the total angular momentum for a system of objects
    '''
    numObjects=len(massArray)
    numDimensions=len(posArray[0])
    momentum=[[0.0]*numDimensions for j in range(numObjects)]
    angMomen=[[0.0]*numDimensions for j in range(numObjects)]
    angMomenTot=[0.0]*numDimensions
    
    for j in range(numObjects):
        for i in range(numDimensions):
            momentum[j][i]=massArray[j]*velArray[j][i]
        angMomen[j]=np.cross(posArray[j],momentum[j])
    for j in range(numObjects):
        for i in range(numDimensions):
            angMomenTot[i]+=angMomen[j][i]
    return angMomenTot


def total_energy(posArray, velArray, massArray):
    '''calculates the total energy (gravitational and kinetic) for a system of objects
    '''
    numObjects=len(massArray)
    numDimensions=len(posArray[0])
    G=math.pow(2.0*math.pi,2.0)
    gravPotential=0.0
    kineticEnergy=0.0
    speedArray=np.linalg.norm(velArray,axis=1)
    
    for j in range(numObjects):
        for k in range(j+1,numObjects):
            r = pythag_pos_diff(posArray[j],posArray[k])
            gravPotential+=(-G*massArray[j]*massArray[k])/r
    for j in range(numObjects):
        kineticEnergy+=(1.0/2.0)*massArray[j]*(speedArray[j]**2.0)
    
    energyTot=gravPotential+kineticEnergy
    return energyTot
            

numObjects=2
numDimensions=3
#Quantiy      Units
#time      |   Year
#distance  |   AU
#velocity  |   AU/Year
#mass      |   Solar Mass

#Number       Object
#0         |   Sun
#1         |   Jupiter
#2         |   Saturn

#generate velocity and position vectors
#the arrays are defined as ____Array[jth object][dimensions x,y,z: 0,1,2]
#velocity array is half a timestep ahead of the position and acceleration arrays
posArray=[[0.0]*numDimensions for j in range(numObjects)]
velArray=[[0.0]*numDimensions for j in range(numObjects)]

#positions and velocitys of solar system bodies in baryocentric coords
#all starting values at A.D. 2016-Feb-29 00:00:00.0000
posArray[0]=[3.769360061446780E-03,  1.812773374712511E-03, -1.624014729877611E-04]
posArray[1]=[-5.291527066741596E+00,  1.182309344528905E+00,  1.134242057706169E-01]
#posArray[2]=[-3.421132666518896E+00, -9.407865807544542E+00,  2.997343479755205E-01]

velArray[0]=[2.116907618849916E-07*365,  6.980633760687217E-06*365, -1.203681274902981E-08*365]
velArray[1]=[-1.733768918092524E-03*365, -7.008096404329502E-03*365,  6.792695786000429E-05*365]
#velArray[2]=[4.937187945581996E-03*365, -1.922883900630874E-03*365, -1.633640215930448E-04*365]

#generate massArray[jth object]
massArray=[1.0, 0.0009543]#,0.0002857

print posArray
print velArray
print massArray

stepSize=25/365.25 #1 day
stepCount=250*3


historicPosArray = nbody_sim_euler(posArray,velArray,massArray,stepSize,stepCount)
plt.plot(historicPosArray[0][0],historicPosArray[0][1],'.') 
plt.plot(historicPosArray[1][0],historicPosArray[1][1],'.--')
#plt.plot(historicPosArray[2][0],historicPosArray[1][1],'.--')
#plt.show()
#print historicPosArray[1][0]

