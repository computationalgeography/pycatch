#!/usr/bin/python
## example2_14
import numpy
from LUdecomp3 import *

n = 6
d = numpy.ones((n))*2.0
e = numpy.ones((n-1))*(-1.0)
c = e.copy()
d[n-1] = 5.0
aInv = numpy.identity(n)*1.0
c,d,e = LUdecomp3(c,d,e)
for i in range(n):
    aInv[:,i] = LUsolve3(c,d,e,aInv[:,i])
print "\nThe inverse matrix is:\n",aInv
