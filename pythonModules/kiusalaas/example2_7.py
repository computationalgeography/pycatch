#!/usr/bin/python
## example2_7
#from numarray import zeros,array,Float64,product,diagonal
import numpy
from LUdecomp import *

a = numpy.array([[ 3.0, -1.0,  4.0], \
           [-2.0,  0.0,  5.0], \
           [ 7.0,  2.0, -2.0]])
b = numpy.array([[ 6.0,  3.0,  7.0], \
           [-4.0,  2.0, -5.0]])
a = LUdecomp(a)
det = numpy.prod(numpy.diag(a))
print "\nDeterminant =",det
for i in range(len(b)):
   x = LUsolve(a,b[i])
   print "x",i+1,"=",x
raw_input("\nPress return to exit")

