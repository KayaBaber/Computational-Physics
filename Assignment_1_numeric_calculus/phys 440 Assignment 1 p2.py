# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 22:07:57 2016

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
  stepSize=l/(float(length))

  arrayIntegral=[initialCond]
  trapSum=initialCond
  for i in range(1, length):
    trapSum = (trapSum + (stepSize/2)*(arrayIn[i]+arrayIn[i-1]))
    arrayIntegral.append(trapSum)
  
  return arrayIntegral

  
  
def simpson(arrayIn,length):
  #computes the definite intergral of the array with simpson's rule
  #returns definiteIntegral out, a number equal to the total area under the function
  
    h = float(length) / float(len(arrayIn))
    s = arrayIn[0] + arrayIn[len(arrayIn)-1]

    for i in range(1, len(arrayIn), 2):
        s += 4 * arrayIn[i]
    for i in range(2, len(arrayIn)-1, 2):
        s += 2 * arrayIn[i]

    return s * h / 3




l=5.0
  
maxErrorScaled=[]
Narray=[]
maxErrorScaledLog=[]
NarrayLog=[]
simpsonErrorScaled=[]
for p in range(2,8):
  N=math.pow(4,p)
  Narray.append(N)
  NarrayLog.append(math.log10(N))
  testArray=[]
  knownArray=[]
  xArray=[]
  stepSize=l/(N-1)
  for i in range(int(N)):
    xArray.append(i*stepSize)
    testArray.append(math.pow(math.e,xArray[i]))
    knownArray.append(math.pow(math.e,xArray[i]))

    computeArray=trap_rule(testArray,l,0)

    simpsonSum=simpson(testArray,l)
  errorArray=[]

  for i in range(int(N)-1):
    errorArray.append(abs(knownArray[i]-computeArray[i]))
  simpsonErrorScaled.append(abs(knownArray[-1]-simpsonSum)/max(knownArray))
  maxErrorScaled.append(max(errorArray)/max(knownArray))
  maxErrorScaledLog.append(math.log10(max(errorArray)/max(knownArray)))

#plt.plot(xArray,computeArray,'.')
#plt.plot(xArray,knownArray)
plt.plot(Narray,maxErrorScaled,marker='o',linestyle='--',label='Trap Rule Method')
plt.plot(Narray,simpsonErrorScaled,marker='s',linestyle='--',label='Simpsons Rule Method')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':20})
plt.loglog()
plt.xlabel('Number of Subintervals Used',fontsize=20)
plt.ylabel('Scaled Maximum Error',fontsize=20)
plt.grid()
plt.show()  