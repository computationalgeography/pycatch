#!/usr/bin/python
## example2_13
import numpy
from LUpivot import *

def matInv(a):
    n = len(a[0])
    aInv = numpy.identity(n)*1.0
    a,seq = LUdecomp(a)
    for i in range(n):
        aInv[:,i] = LUsolve(a,aInv[:,i],seq)
    return aInv    
    
a = numpy.array([[ 0.6, -0.4,  1.0],\
           [-0.3,  0.2,  0.5],\
           [ 0.6, -1.0,  0.5]])
aOrig = a.copy()  # Save original [a]
aInv = matInv(a)  # Invert [a] (original [a] is destroyed)
print "\naInv =\n",aInv
print "\nCheck: a*aInv =\n", numpy.dot(aOrig,aInv)
           
    
