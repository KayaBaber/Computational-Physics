import math
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()

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
    historicAngMomenDrift=[]
    historicEngDrift=[]
    
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
                
        #compute drift
        ang=angular_momentum(posArray, velArray, massArray)
        energy=total_energy(posArray, velArray, massArray)
        angMomenDrift=np.subtract(ang,angStart)
        historicAngMomenDrift.append((angMomenDrift[0]**2 + angMomenDrift[1]**2 + angMomenDrift[2]**2)**(1.0/2.0))
        historicEngDrift.append((energy-energyStart))
    print 'Execution Duration: '+str(datetime.now() - startTime)
    speedup=((3.154E+07)*stepSize*stepCount)/((datetime.now() - startTime).total_seconds())
    print 'Speedup: '+str(speedup)
    print 'Total drift'
    print np.subtract(ang,angStart), (energy-energyStart)
    return historicPosArray, historicAngMomenDrift, historicEngDrift
    

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
    speedArray=np.linalg.norm(velArray, axis=1)
    
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

#generate massArray[jth object]
massArray=[1.0, 0.0009543]#,0.0002857



#test center of mass coords
# x_sun = -M_jup*x_jup
# T = ((x_jup - x_sun) / (M_jup+M_sun))^(1/2)
# v_jup=2*pi*x_jup / T
# v_sun = -M_jup*v_jup
jupOrbit=5.0
sunOrbit=-jupOrbit*massArray[1]
period= (((jupOrbit-sunOrbit)**3)/(1.0+massArray[1]))**(1.0/2.0)
jupSpeed=2.0*math.pi*jupOrbit / period
sunSpeed=-massArray[1]*jupSpeed

posArray[1]=[jupOrbit,0.0,0.0]
velArray[1]=[0.0,jupSpeed,0.0]
posArray[0]=[sunOrbit,0.0,0.0]
velArray[0]=[0.0,sunSpeed,0.0]

print '\n\n--------------------------------------------------\nEULER CIRCULAR\n--------------------------------------------------'
#print 'Start Postions',posArray
#print 'Start Velocities',velArray
#print 'Masses',massArray

stepSize=5/365.25 #1 day
stepCount=20000
print 'Step Size: ' +str(stepSize*365.25) +' days'
print 'Steps: '+str(stepCount)
print 'Simulation Duration: '+str(stepSize*stepCount) +' years OR '+str(stepSize*stepCount*365.25)+' days'

for i in range(1,10):
    historicPosArray,historicAngMomenDrift, historicEngDrift = nbody_sim_euler(posArray,velArray,massArray,stepSize*i,stepCount)
    #plt.plot(historicPosArray[0][0],historicPosArray[0][1],'.') 
    #plt.plot(historicPosArray[1][0],historicPosArray[1][1],'.--')
    #plt.plot(historicPosArray[2][0],historicPosArray[1][1],'.--'
    #plt.gca().set_aspect('equal', adjustable='box')

    #plt.show()
    plt.clf()
    plt.plot(historicAngMomenDrift)
    plt.title('Step Size: ' +str(stepSize*i*365.25) +' days')
    plt.savefig('/Users/student/kbaber/Desktop/Phys440//angMomenDiffPlot_stepsize_'+str(stepSize*i*365.25)+'.png',bbox_inches='tight')
    #plt.plot(historicEngDrift)
    #plt.show()


