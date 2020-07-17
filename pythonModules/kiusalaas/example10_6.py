## example10_6
from fletcherReeves import *
import numpy
import math

def F(x):
    return  8.0/x[0] - x[0]*(math.tan(x[1]) - 2.0/math.cos(x[1]))

def gradF(x):
    g = numpy.zeros((2),dtype=float)
    g[0] = -8.0/(x[0]**2) - math.tan(x[1]) + 2.0/math.cos(x[1])
    g[1] = x[0]*(-1.0/math.cos(x[1]) + 2.0*math.tan(x[1]))/math.cos(x[1])
    return g

x = numpy.array([2.0, 0.0])
x,nIter = optimize(F,gradF,x)
b = 8.0/x[0] - x[0]*math.tan(x[1])
print "h =",x[0],"m"
print "b =",b,"m"
print "theta =",x[1]*180.0/math.pi,"deg"
print "perimeter =", F(x),"m"
print "Number of iterations =",nIter
raw_input("Press return to exit")
