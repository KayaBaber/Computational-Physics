'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 5 - PDEs
Exploration 2 - Parabolic PDEs: The Wave Equation
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math



def make_banded(N,M):
    bandTopBot = [-1.]*(N-1)
    bandMid = [2. + (4.*M)/(N**2) ]*N
    banded = np.diag(bandMid)
    banded = np.add(banded,np.diag(bandTopBot,1))
    banded = np.add(banded,np.diag(bandTopBot,-1))
    return banded
    
    
def make_operator(N,M):
    bandedCrank = make_banded(N,M)
    negativeCrank = bandedCrank * (-1)
    bandedCrank[0] = [1] + [0]*(N-1)
    bandedCrank[-1] = [0]*(N-1) + [1]
    invertedCrank = LA.inv(bandedCrank)
    operatorCrank = invertedCrank.dot(negativeCrank)
    return operatorCrank

