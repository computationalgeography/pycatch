## example8_6
import numpy
from LUdecomp3 import *
import math

def equations(x,h,m): # Set up finite difference eqs.
    h2 = h*h
    d = numpy.ones((m + 1))*(-2.0 + 4.0*h2)
    c = numpy.ones((m),dtype = float)
    e = numpy.ones((m),dtype = float)
    b = numpy.ones((m+1))*4.0*h2*x
    d[0] = 1.0
    e[0] = 0.0
    b[0] = 0.0
    c[m-1] = 2.0
    return c,d,e,b

xStart = 0.0        # x at left end
xStop = math.pi/2.0      # x at right end
m = 10              # Number of mesh spaces
h = (xStop - xStart)/m
x = numpy.arange(xStart,xStop + h,h)
c,d,e,b = equations(x,h,m)
c,d,e = LUdecomp3(c,d,e)
y = LUsolve3(c,d,e,b)
print "\n        x              y"
for i in range(m + 1):
    print "%14.5e %14.5e" %(x[i],y[i])
raw_input("\nPress return to exit")


