#!/usr/bin/python
## example2_18

import numpy
from conjGrad import *

def Ax(v):
    n = len(v)
    Ax = numpy.zeros((n),dtype=float)
    Ax[0] = 2.0*v[0] - v[1]+v[n-1]
    Ax[1:n-1] = -v[0:n-2] + 2.0*v[1:n-1] -v [2:n]
    Ax[n-1] = -v[n-2] + 2.0*v[n-1] + v[0]
    return Ax

n = eval(raw_input("Number of equations ==> "))
b = numpy.zeros((n),dtype=float)
b[n-1] = 1.0
x = numpy.zeros((n),dtype=float)
x,numIter = conjGrad(Ax,x,b)
print "\nThe solution is:\n",x
print "\nNumber of iterations =",numIter
    
