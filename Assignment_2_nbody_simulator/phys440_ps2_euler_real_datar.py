import math
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()


def angular_momentum(posArray, velArray, massArray):
    '''calculates the total angular momentum for a system of objects
    '''
    numObjects=len(massArray)
    numDimensions=len(posArray[0])
    momentum=[[0.0]*numDimensions for j in range(numObjects)]
    angMomen=[[0.0]*numDimensions for j in range(numObjects)]
    angMomenTot=[0.0]*numDimensions
    
    for j in range(numObjects):
        #take cross product of postion and velocity
        angMomen[j]=np.cross(posArray[j],velArray[j])
        #multiply each component by mass
        for i in range(numDimensions):
            angMomen[j][i]=massArray[j]*angMomen[j][i]
        #print angMomen[j]
    #sum up all angular momenta
    for j in range(numObjects):
        for i in range(numDimensions):
            angMomenTot[i]=angMomen[j][i]+angMomenTot[i]
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
    angStart=angular_momentum(posArray, velArray, massArray)
    energyStart=total_energy(posArray, velArray, massArray)
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
        #print accelArray[1]   
        #update velocity
        for j in range(numObjects):
            for l in range(numDimensions):
                velArray[j][l]+=stepSize*accelArray[j][l]
                historicVelArray[j][l].append(velArray[j][l])
        #print ang, energy
        #print velArray[1]
    ang=angular_momentum(posArray, velArray, massArray)
    energy=total_energy(posArray, velArray, massArray)
    print 'total drift\nTotal Angular Momentum\n       x                  y                  z         Total Energy'
    print np.subtract(ang,angStart), (energy-energyStart)
    return historicPosArray, historicVelArray
    
    
def plot_planets(historicPosArray,numObjects, numDimensions):
    titles=["X-Y Top-Down View","Y-Z"]
    #for j in range(numDimensions-1): 
    for i in range(numObjects):
        plt.plot(historicPosArray[i][0],historicPosArray[i][1],'.-')
            
    

numObjects=9
numDimensions=3
#Quantiy      Units
#time      |   Year
#distance  |   AU
#velocity  |   AU/Year
#mass      |   Solar Mass

#Number       Object
#0         |   Sun
#1         |   Mercury
#2         |   Venus
#3         |   Earth + Moon
#4         |   Mars
#5         |   Jupiter
#6         |   Saturn
#7         |   Uranus
#8         |   Neptune

#generate velocity and position vectors
#the arrays are defined as ____Array[jth object][dimensions x,y,z: 0,1,2]
#velocity array is half a timestep ahead of the position and acceleration arrays
posArray=[[0.0]*numDimensions for j in range(numObjects)]
velArray=[[0.0]*numDimensions for j in range(numObjects)]


#generate massArray[jth object]
massArray=[1.0, 1.661E-07, 2.4484E-06, 3.00346E-06, 3.2279E-07, 0.0009543, 0.0002857, 4.36576E-05, 5.15E-05]

#positions and velocitys of solar system bodies in baryocentric coords
#all starting values at A.D. 2016-Feb-29 00:00:00.0000
#SUN
posArray[0]=[3.769360061446780E-03,  1.812773374712511E-03, -1.624014729877611E-04]
velArray[0]=[2.116907618849916E-07*365,  6.980633760687217E-06*365, -1.203681274902981E-08*365]

#MERCURY
posArray[1]=[6.384933295823990E-02, -4.515747976027313E-01, -4.272075514855849E-02]
velArray[1]=[2.224830448309151E-02*365,  5.147225564402190E-03*365, -1.621110695570192E-03*365]

#VENUS
posArray[2]=[1.167830458441351E-01, -7.165914711359284E-01, -1.653385613745591E-02]
velArray[2]=[1.984550800517364E-02*365,  3.079327873621045E-03*365, -1.103086914976536E-03*365]

#EARTH
posArray[3]=[-9.247648258319329E-01,  3.469717738912940E-01, -1.730700077332298E-04]
velArray[3]=[-6.279271560859209E-03*365, -1.617967468331688E-02*365,  2.694319332791631E-07*365]

#MARS
posArray[4]=[-1.510114581167425E+00, -5.693088681586485E-01,  2.502401748406479E-02]
velArray[4]=[5.461536539056124E-03*365, -1.189036297311848E-02*365, -3.833644710566484E-04*365]

#JUPITER
posArray[5]=[-5.291527066741596E+00,  1.182309344528905E+00,  1.134242057706169E-01]
velArray[5]=[-1.733768918092524E-03*365, -7.008096404329502E-03*365,  6.792695786000429E-05*365]

#SATURN
posArray[6]=[-3.421132666518896E+00, -9.407865807544542E+00,  2.997343479755205E-01]
velArray[6]=[4.937187945581996E-03*365, -1.922883900630874E-03*365, -1.633640215930448E-04*365]

#URANUS
posArray[7]=[1.879301283705500E+01,  6.762792046709805E+00, -2.183506612863105E-01]
velArray[7]=[-1.360455775595364E-03*365,  3.517436778076412E-03*365,  3.073947954095653E-05*365]

#NEPTUNE
posArray[8]=[ 2.801943446698039E+01, -1.060190333542911E+01, -4.274100493913166E-01]
velArray[8]=[1.090032425561670E-03*365,  2.955023769315658E-03*365, -8.560834225044430E-05*365]

#print 'Start Postions',posArray
#print 'Start Velocities',velArray
#print 'Masses',massArray

stepSize=1.0/(365.25)#*100000000*0.928514127163)
stepCount=200
print 'Step Size: ' +str(stepSize*365.25) +' days'
print 'Simulation Duration: '+str(stepSize*stepCount) +' years OR '+str(stepSize*stepCount*365.25)+' days'



historicPosArray, historicVelArray = nbody_sim_euler(posArray,velArray,massArray,stepSize,stepCount)
#plot_planets(historicPosArray, numObjects, numDimensions)
#plt.plot(historicPosArray[0][0],historicPosArray[0][1],'.') 
#plt.show()
#plt.plot(historicVelArray[0][0],'.')
#plt.show()
#plt.plot(historicPosArray[1][0],historicPosArray[1][1],'.-')
#plt.plot(historicPosArray[2][0],historicPosArray[2][1],'.-')
#plt.gca().set_aspect('equal', adjustable='box')
print 'Execution Duration: '+str(datetime.now() - startTime)
speedup=((3.154E+07)*stepSize*stepCount)/((datetime.now() - startTime).total_seconds())
print 'Speedup: '+str(speedup)