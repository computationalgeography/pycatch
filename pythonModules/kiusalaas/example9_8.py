## example9_8
import numpy
from householder import *

a = numpy.array([[ 7.0, 2.0,  3.0, -1.0],  \
           [ 2.0, 8.0,  5.0,  1.0],  \
           [ 3.0, 5.0, 12.0,  9.0],  \
           [-1.0, 1.0,  9.0,  7.0]])
d,c = householder(a)
print "Principal diagonal {d}:\n", d
print "\nSubdiagonal {c}:\n",c
print "\nTransformation matrix [P]:"
print computeP(a)
raw_input("\nPress return to exit")


