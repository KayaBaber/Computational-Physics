# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 23:04:31 2016

@author: Kaya
"""

#Phys 440 Assignment 1
#Kaya Baber

#Problem 2

import math
import matplotlib.pyplot as plt

def trap_rule(arrayIn, l, initialCond):
  #takes arrayIn of equally spaced samples of a function, a range l, and a initial condition for the first point
  #computes the indefinite integral of the array with the trap rule
  #returns arrayIntegrate, an array of the integral up to that point

  length=len(arrayIn)
  stepSize=l/(length)

  arrayIntegral=[initialCond]
  trapSum=initialCond
  for i in range(1, length):
    trapSum = (trapSum + (stepSize/2)*(arrayIn[i]+arrayIn[i-1]))
    arrayIntegral.append(trapSum)
  
  return arrayIntegral

  
  
def simpson(arrayIn,length):
  #computes the definite intergral of the array with simpson's rule
  #returns definite Integral out, a number equal to the total area under the function
  
    h = float(length) / float(len(arrayIn))
    s = arrayIn[0] + arrayIn[len(arrayIn)-1]

    for i in range(1, len(arrayIn), 2):
        s += 4 * arrayIn[i]
    for i in range(2, len(arrayIn)-1, 2):
        s += 2 * arrayIn[i]

    return s * h / 3

#nth derivative function using pth order method
  #for n=(1,2) and p=(2,4)
  #takes input of an array arrayIn, n, p, and a range l
  #outputs an array arrayDeriv that is the nth order derivative of arrayIn
def numericDeriv(arrayIn, n, p, l):
  length=len(arrayIn)
  stepSize = l/(float(length))
  arrayDeriv=[]
  
  
  
  
  #define coefficients and calculate differences
  if p==2:
    if n==1:
      centralCoefficients=[-1.0/2.0, 0, 1.0/2.0]
      forBackCoefficients=[-3.0/2.0, 2.0, -1.0/2.0]
      
      #first do start and end points with forward/backward-difference
      localDeriv=( (arrayIn[0]*forBackCoefficients[0]) + (arrayIn[1]*forBackCoefficients[1]) + (arrayIn[2]*forBackCoefficients[2]) )/stepSize
      lastDeriv=( -(arrayIn[length-1]*forBackCoefficients[0]) - (arrayIn[length-2]*forBackCoefficients[1]) - (arrayIn[length-3]*forBackCoefficients[2]) )/stepSize
      arrayDeriv.append(localDeriv)
      
      #next do all middle points with central-difference
      for i in range(1,length-2):
        localDeriv=( (arrayIn[i-1]*centralCoefficients[0]) + (arrayIn[i]*centralCoefficients[1]) + (arrayIn[i+1]*centralCoefficients[2]) )/stepSize
        arrayDeriv.append(localDeriv)
      
    elif n==2:
      centralCoefficients=[1.0, -2.0, 1.0]
      forBackCoefficients=[2.0, -5.0, 4.0, -1.0]
      
      #first do start and end points with forward/backward-difference
      localDeriv=( (arrayIn[0]*forBackCoefficients[0]) + (arrayIn[1]*forBackCoefficients[1]) + (arrayIn[2]*forBackCoefficients[2]) + (arrayIn[3]*forBackCoefficients[3]) )/(stepSize*stepSize)
      lastDeriv=-( (arrayIn[length-1]*forBackCoefficients[0]) + (arrayIn[length-2]*forBackCoefficients[1]) + (arrayIn[length-3]*forBackCoefficients[2]) + (arrayIn[length-4]*forBackCoefficients[3]) )/(stepSize*stepSize)
      arrayDeriv.append(localDeriv)
      
      #next do all middle points with central-difference
      for i in range(1,length-2):
        localDeriv=( (arrayIn[i-1]*centralCoefficients[0]) + (arrayIn[i]*centralCoefficients[1]) + (arrayIn[i+1]*centralCoefficients[2]) )/(stepSize*stepSize)
        arrayDeriv.append(localDeriv)
 
    #add end points(s) with backward-difference to arrayDeriv
    arrayDeriv.append(lastDeriv)
    
    
    
     
  elif p==4:
    if n==1:
      centralCoefficients=[1.0/12.0, -2.0/3.0, 0, 2.0/3.0, -1.0/12.0]
      forBackCoefficients=[-25.0/12.0, 4.0,	-3.0,	4.0/3.0, -1.0/4.0]
      
      #first do start and end points with forward/backward-difference
      firstDeriv=( (arrayIn[0]*forBackCoefficients[0]) + (arrayIn[1]*forBackCoefficients[1]) + (arrayIn[2]*forBackCoefficients[2]) + (arrayIn[3]*forBackCoefficients[3]) + (arrayIn[4]*forBackCoefficients[4]))/stepSize
      localDeriv=( (arrayIn[1]*forBackCoefficients[0]) + (arrayIn[2]*forBackCoefficients[1]) + (arrayIn[3]*forBackCoefficients[2]) + (arrayIn[4]*forBackCoefficients[3]) + (arrayIn[5]*forBackCoefficients[4]))/stepSize
      secondLastDeriv= ( -(arrayIn[length-2]*forBackCoefficients[0]) - (arrayIn[length-3]*forBackCoefficients[1]) - (arrayIn[length-4]*forBackCoefficients[2]) - (arrayIn[length-5]*forBackCoefficients[3]) - (arrayIn[length-6]*forBackCoefficients[4])  )/stepSize
      lastDeriv=( -(arrayIn[length-1]*forBackCoefficients[0]) - (arrayIn[length-2]*forBackCoefficients[1]) - (arrayIn[length-3]*forBackCoefficients[2]) - (arrayIn[length-4]*forBackCoefficients[3]) - (arrayIn[length-5]*forBackCoefficients[4]) )/stepSize
      arrayDeriv.append(firstDeriv)
      arrayDeriv.append(localDeriv)
      
      #next do all middle points with central
      for i in range(2,length-3):
        localDeriv=( (arrayIn[i-2]*centralCoefficients[0]) + (arrayIn[i-1]*centralCoefficients[1]) + (arrayIn[i]*centralCoefficients[2]) + (arrayIn[i+1]*centralCoefficients[3]) + (arrayIn[i+2]*centralCoefficients[4]) )/stepSize
        arrayDeriv.append(localDeriv)
      
      
    elif n==2:
      centralCoefficients=[-1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0]
      forBackCoefficients=[15.0/4.0, -77.0/6.0, 107.0/6.0, -13.0, 61.0/12.0, -5.0/6.0]
      #first do start and end points with forward/backward-difference
      firstDeriv=( (arrayIn[0]*forBackCoefficients[0]) + (arrayIn[1]*forBackCoefficients[1]) + (arrayIn[2]*forBackCoefficients[2]) + (arrayIn[3]*forBackCoefficients[3]) + (arrayIn[4]*forBackCoefficients[4]) + (arrayIn[5]*forBackCoefficients[5]) )/(stepSize*stepSize)
      localDeriv=( (arrayIn[1]*forBackCoefficients[0]) + (arrayIn[2]*forBackCoefficients[1]) + (arrayIn[3]*forBackCoefficients[2]) + (arrayIn[4]*forBackCoefficients[3]) + (arrayIn[5]*forBackCoefficients[4]) + (arrayIn[6]*forBackCoefficients[5]))/(stepSize*stepSize)
      secondLastDeriv= ( (arrayIn[length-2]*forBackCoefficients[0]) + (arrayIn[length-3]*forBackCoefficients[1]) + (arrayIn[length-4]*forBackCoefficients[2]) + (arrayIn[length-5]*forBackCoefficients[3]) + (arrayIn[length-6]*forBackCoefficients[4]) + (arrayIn[length-7]*forBackCoefficients[5])  )/(stepSize*stepSize)
      lastDeriv= ( (arrayIn[length-1]*forBackCoefficients[0]) + (arrayIn[length-2]*forBackCoefficients[1]) + (arrayIn[length-3]*forBackCoefficients[2]) + (arrayIn[length-4]*forBackCoefficients[3]) + (arrayIn[length-5]*forBackCoefficients[4]) + (arrayIn[length-6]*forBackCoefficients[5]) )/(stepSize*stepSize)
      arrayDeriv.append(firstDeriv)
      arrayDeriv.append(localDeriv)
    
      #next do all middle points with central
      for i in range(2,length-3):
        localDeriv=( (arrayIn[i-2]*centralCoefficients[0]) + (arrayIn[i-1]*centralCoefficients[1]) + (arrayIn[i]*centralCoefficients[2]) +(arrayIn[i+1]*centralCoefficients[3]) + (arrayIn[i+2]*centralCoefficients[4]))/stepSize
        arrayDeriv.append(localDeriv)
    
    #add end points(s) with backward-difference to arrayDeriv
    arrayDeriv.append(secondLastDeriv)
    arrayDeriv.append(lastDeriv)
    
    
  return arrayDeriv
  
  
  



def lineChargeZpot(z,Q,a):
  l=20.0
  N=100.0
  chargeArray=[]
  chargeDensArray=[]
  xArray=[]
  stepSize=l/(N-1)
    
  for i in range(-int(N/2),int(N/2)):
    xArray.append(i*stepSize)
  for i in range(int(N)):
    chargeArray.append(math.pow(math.e,-math.pow(xArray[i]/a,2))/math.pow(math.pow(xArray[i],2)+math.pow(z,2),(1.0/2.0)))
    #chargeArray.append(1.0/ math.pow(math.pow(xArray[i],2)+math.pow(z,2),(1.0/2.0)))
    #chargeDensArray.append((Q/a))
    chargeDensArray.append((Q/a)*math.pow(math.e,-math.pow(xArray[i]/a,2)))
  potential=(Q/a)*simpson(chargeArray,l)
  
  
  return (potential,chargeDensArray,xArray)


H=100
stepSize=1.0
zPotential=[]
zEfield=[]
expectLog=[]
for i in range(1,H):
  zPotential.append(lineChargeZpot(i*stepSize,1.0,1.0)[0])
  #expectLog.append(4*math.log(math.e,((i+6)/4.0))-1)
zEfield=numericDeriv(zPotential, 1, 2, 20)
print zPotential
#print expectLog

#plt.plot(zPotential)
#plt.plot(expectLog)
plt.plot(lineChargeZpot(1.0,1.0,1.0)[2],lineChargeZpot(1.0,1.0,1.0)[1])
#plt.plot(zEfield)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':20})
plt.xlabel('X Position on Charged Line',fontsize=20)
plt.ylabel('Linear Charge Density',fontsize=20)

#plt.ylabel('Potential',fontsize=20)
#plt.ylabel('E-Field_z',fontsize=20)
plt.grid()
plt.show()