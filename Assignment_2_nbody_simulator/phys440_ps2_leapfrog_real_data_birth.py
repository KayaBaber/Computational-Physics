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


def nbody_sim_leapfrog(posArray,velArray, massArray,stepSize,stepCount):
    '''Returns an array of computed position steps for the inputed objects
    '''
    numObjects=len(massArray)
    numDimensions=len(posArray[0])
    angStart=angular_momentum(posArray, velArray, massArray)
    energyStart=total_energy(posArray, velArray, massArray)
    #0. calculate first velocity step with forward euler
    #1. update position
    #2. compute acceleration
    #3. update velocity
    G=math.pow(2.0*math.pi,2.0)
    historicPosArray=[[[]for j in range(numDimensions)] for j in range(numObjects)]
    historicVelArray=[[[]for j in range(numDimensions)] for j in range(numObjects)]
    historicAngMomenDrift=[]
    historicEngDrift=[]
    
    #start with forward Euler
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
            velArray[j][l]+=(1.0/2.0)*stepSize*accelArray[j][l]
            historicVelArray[j][l].append(velArray[j][l])
            
    #simulation loop using leap frog
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
                
        #compute drift
        ang=angular_momentum(posArray, velArray, massArray)
        energy=total_energy(posArray, velArray, massArray)
        angMomenDrift=np.subtract(ang,angStart)
        historicAngMomenDrift.append((angMomenDrift[0]**2 + angMomenDrift[1]**2 + angMomenDrift[2]**2)**(1.0/2.0))
        historicEngDrift.append((energy-energyStart))

    print 'Execution Duration: '+str(datetime.now() - startTime)
    speedup=((3.154E+07)*stepSize*stepCount)/((datetime.now() - startTime).total_seconds())
    print 'Speedup: '+str(speedup)
    print 'total drift'
    print np.subtract(ang,angStart), (energy-energyStart)
    return historicPosArray, historicAngMomenDrift, historicEngDrift
    
    
def plot_planets(historicPosArray,numObjects, numDimensions):

    #for j in range(numDimensions-1): 
    for i in range(1):
        plt.plot(historicPosArray[i][0],historicPosArray[i][1],'.-')
    plt.gca().set_aspect('equal', adjustable='box')
    

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
#all starting values at A.D. 1994-Mar-12 00:00:00.0000
#SUN
posArray[0]=[7.971719296098350E-04,  6.530924435763692E-03, -3.840865188581198E-05]
velArray[0]=[-5.057483005901317E-06*365,  3.950592533314976E-06*365,  1.250632236911193E-07*365]

#MERCURY
posArray[1]=[-3.168209006763542E-01, -3.100284070655648E-01,  3.261781129133621E-03] 
velArray[1]=[1.413036973396787E-02*365, -1.863818179319557E-02*365, -2.820247511064388E-03*365]

#VENUS
posArray[2]=[ 6.703870134057690E-01,  2.820919488225764E-01, -3.492924483480292E-02] 
velArray[2]=[-7.770139352355757E-03*365,  1.861901186591109E-02*365,  7.024925476467198E-04*365]

#EARTH
posArray[3]=[-9.813301325021555E-01,  1.574386351629199E-01, -3.957367973686485E-05] 
velArray[3]=[-2.899220177422116E-03*365, -1.707155547951476E-02*365,  1.254070235790866E-07*365]

#MARS
posArray[4]=[ 1.078786065897782E+00, -8.625912141839485E-01, -4.475462451602934E-02] 
velArray[4]=[ 9.309976370119377E-03*365,  1.209485051343706E-02*365,  2.432439495215345E-05*365]

#JUPITER
posArray[5]=[ -4.391049819964502E+00, -3.196505626109847E+00,  1.115744883801592E-01] 
velArray[5]=[ 4.346257193870669E-03*365, -5.749413180499637E-03*365, -7.354330890512936E-05*365]

#SATURN
posArray[6]=[ 8.733336994200078E+00, -4.345364839040485E+00, -2.714746726539586E-01]
velArray[6]=[ 2.183425528529507E-03*365,  4.983355546867565E-03*365, -1.738662890388426E-04*365]

#URANUS
posArray[7]=[ 7.679161474619261E+00, -1.807440948712653E+01, -1.666412541349566E-01] 
velArray[7]=[ 3.590748782173861E-03*365,  1.355454603714086E-03*365, -4.154112618083728E-05*365]

#NEPTUNE
posArray[8]=[ 1.095828259071641E+01, -2.811168975971997E+01,  3.263530816155578E-01]
velArray[8]=[  2.904104193585754E-03*365,  1.157377132972549E-03*365, -9.086645981465949E-05*365]



posArrayF=[[0.0]*numDimensions for j in range(numObjects)]
velArrayF=[[0.0]*numDimensions for j in range(numObjects)]
#positions and velocitys of solar system bodies in baryocentric coords
#all final values at A.D. 2016-Feb-29 00:00:00.0000
#SUN
posArrayF[0]=[3.769360061446780E-03,  1.812773374712511E-03, -1.624014729877611E-04]
velArrayF[0]=[2.116907618849916E-07*365,  6.980633760687217E-06*365, -1.203681274902981E-08*365]

#MERCURY
posArrayF[1]=[6.384933295823990E-02, -4.515747976027313E-01, -4.272075514855849E-02]
velArrayF[1]=[2.224830448309151E-02*365,  5.147225564402190E-03*365, -1.621110695570192E-03*365]

#VENUS
posArrayF[2]=[1.167830458441351E-01, -7.165914711359284E-01, -1.653385613745591E-02]
velArrayF[2]=[1.984550800517364E-02*365,  3.079327873621045E-03*365, -1.103086914976536E-03*365]

#EARTH
posArrayF[3]=[-9.247648258319329E-01,  3.469717738912940E-01, -1.730700077332298E-04]
velArrayF[3]=[-6.279271560859209E-03*365, -1.617967468331688E-02*365,  2.694319332791631E-07*365]

#MARS
posArrayF[4]=[-1.510114581167425E+00, -5.693088681586485E-01,  2.502401748406479E-02]
velArrayF[4]=[5.461536539056124E-03*365, -1.189036297311848E-02*365, -3.833644710566484E-04*365]

#JUPITER
posArrayF[5]=[-5.291527066741596E+00,  1.182309344528905E+00,  1.134242057706169E-01]
velArrayF[5]=[-1.733768918092524E-03*365, -7.008096404329502E-03*365,  6.792695786000429E-05*365]

#SATURN
posArrayF[6]=[-3.421132666518896E+00, -9.407865807544542E+00,  2.997343479755205E-01]
velArrayF[6]=[4.937187945581996E-03*365, -1.922883900630874E-03*365, -1.633640215930448E-04*365]

#URANUS
posArrayF[7]=[1.879301283705500E+01,  6.762792046709805E+00, -2.183506612863105E-01]
velArrayF[7]=[-1.360455775595364E-03*365,  3.517436778076412E-03*365,  3.073947954095653E-05*365]

#NEPTUNE
posArrayF[8]=[ 2.801943446698039E+01, -1.060190333542911E+01, -4.274100493913166E-01]
velArrayF[8]=[1.090032425561670E-03*365,  2.955023769315658E-03*365, -8.560834225044430E-05*365]

print '\n\n--------------------------------------------------\nLEAPFROG REAL DATA\n--------------------------------------------------'
#print 'Start Postions',posArray
#print 'Start Velocities',velArray
#print 'Masses',massArray

stepSize=1/(365.25)#*100000000*0.928514127163)
stepCount=8025
print 'Step Size: ' +str(stepSize*365.25) +' days'
print 'Simulation Duration: '+str(stepSize*stepCount) +' years OR '+str(stepSize*stepCount*365.25)+' days'



historicPosArray,historicAngMomenDrift, historicEngDrift = nbody_sim_leapfrog(posArray,velArray,massArray,stepSize,stepCount)


plot_planets(historicPosArray, numObjects, numDimensions)
plt.title('Leapfrog Method\nReal Solar System Data from 3/12/1994\n'+str(stepCount)+' steps of '+str(stepSize*365.25)+' days',fontsize=15)
plt.xlabel('AU',fontsize=20)
plt.ylabel('AU',fontsize=20)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

#plt.plot(historicAngMomenDrift)
#plt.title('Real Solar System Data from 3/12/1994\n'+str(stepCount)+' steps of '+str(stepSize*365.25)+' days',fontsize=15)
#plt.xlabel('Steps',fontsize=20)
#plt.ylabel('Total Angular Momentum Drift',fontsize=15)
#plt.show()
#
#plt.plot(historicEngDrift)
#plt.title('Real Solar System Data from 3/12/1994\n'+str(stepCount)+' steps of '+str(stepSize*365.25)+' days',fontsize=15)
#plt.xlabel('Steps',fontsize=20)
#plt.ylabel('Total Energy Drift',fontsize=15)
#plt.show()
#
#diffArray = np.absolute(np.subtract(posArray,posArrayF))
#
#print 'Initial\n', np.array(posArray)
#print '\nActual\n',np.array(posArrayF)
#print '\nAbsolute Difference\n',diffArray
#
#plt.plot(diffArray)
#plt.title('Real Solar System Data from 3/12/1994\n'+str(stepCount)+' steps of '+str(stepSize*365.25)+' days',fontsize=15)
#plt.xlabel('Object',fontsize=20)
#plt.ylabel('TDiffence in position - AU',fontsize=15)
#plt.show()
#for i in range(numObjects):
#    plt.plot(posArray[i][0],posArray[i][1],'.-r')
#    plt.plot(posArrayF[i][0],posArrayF[i][1],'.-b')
#plt.title('Object Position Difference\n'+str(stepCount)+' steps of '+str(stepSize*365.25)+' days',fontsize=15)
#plt.xlabel('AU',fontsize=20)
#plt.ylabel('AU',fontsize=20)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.show()