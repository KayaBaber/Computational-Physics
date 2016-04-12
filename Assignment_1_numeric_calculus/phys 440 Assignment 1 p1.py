# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 20:34:46 2016

@author: Kaya
"""

#Phys 440 Assignment 1
#Kaya Baber

#Problem 1

import math
import matplotlib.pyplot as plt



#nth derivative function using pth order method
  #for n=(1,2) and p=(2,4)
  #takes input of an array arrayIn, n, p, and a range l
  #outputs an array arrayDeriv that is the nth order derivative of arrayIn
def numericDeriv(arrayIn, n, p, l):
  length=len(arrayIn)
  stepSize = l/(length)
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
  
  
#error function
maxErrorScaled12=[]
maxErrorScaled22=[]
maxErrorScaled14=[]
maxErrorScaled24=[]
Narray=[]
maxErrorScaledLog12=[]
maxErrorScaledLog22=[]
maxErrorScaledLog14=[]
maxErrorScaledLog24=[]
NarrayLog=[]
for p in range(2,10):
  N=math.pow(4,p)
  NarrayLog.append(math.log10(N))
  Narray.append(N)
  testArray=[]
  testArrayDeriv2=[]
  testArrayDeriv1=[]
  for i in range(0,int(N)):
    testArray.append(math.sin(2*math.pi*i/N))
    testArrayDeriv2.append(-math.pow((2*math.pi/N),2)*math.sin(2*math.pi*i/N))
    testArrayDeriv1.append((2*math.pi/N)*math.cos(2*math.pi*i/N))
  realArrayDeriv12=numericDeriv(testArray,1,2,int(N))
  realArrayDeriv22=numericDeriv(testArray,2,2,int(N))
  realArrayDeriv14=numericDeriv(testArray,1,4,int(N))
  realArrayDeriv24=numericDeriv(testArray,2,4,int(N))
  errorArray12=[]
  errorArray22=[]
  errorArray14=[]
  errorArray24=[]
  for i in range(int(N)-1):
    errorArray12.append(abs(testArrayDeriv1[i]-realArrayDeriv12[i]))
    errorArray22.append(abs(testArrayDeriv2[i]-realArrayDeriv22[i]))
    errorArray14.append(abs(testArrayDeriv1[i]-realArrayDeriv14[i]))
    errorArray24.append(abs(testArrayDeriv2[i]-realArrayDeriv24[i]))
  maxErrorScaled12.append(max(errorArray12)/max(testArrayDeriv1))
  maxErrorScaled22.append(max(errorArray22)/max(testArrayDeriv2))
  maxErrorScaled14.append(max(errorArray14)/max(testArrayDeriv1))
  maxErrorScaled24.append(max(errorArray24)/max(testArrayDeriv2))
  maxErrorScaledLog12.append(math.log10(max(errorArray12)/max(testArrayDeriv1)))
  maxErrorScaledLog22.append(math.log10(max(errorArray22)/max(testArrayDeriv2)))
  maxErrorScaledLog14.append(math.log10(max(errorArray14)/max(testArrayDeriv1)))
  maxErrorScaledLog24.append(math.log10(max(errorArray24)/max(testArrayDeriv2)))
print '\nmax error scaled:\n'+str(maxErrorScaled22)+'\n'+str(Narray)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(Narray,maxErrorScaled24,marker='o',linestyle='--', label='1st Deriv 2nd Order Accurate')
plt.plot(Narray,maxErrorScaled14,marker='^',linestyle='-.', label='1st Deriv 4th Order Accurate')
plt.plot(Narray,maxErrorScaled22,marker='s',linestyle='--', label='2nd Deriv 2nd Order Accurate')
plt.plot(Narray,maxErrorScaled12,marker='x',linestyle='-.', label='2nd Deriv 4th Order Accurate')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':20})
plt.loglog()
plt.xlabel('Number of Subintervals Used',fontsize=20)
plt.ylabel('Scaled Maximum Error',fontsize=20)
plt.grid()
plt.show()