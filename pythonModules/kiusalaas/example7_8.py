## example7_8

import numpy
import run_kut5
import printSoln
import math

def F(x,y):
    F = numpy.zeros((2),dtype=float) 
    F[0] = y[1]
    F[1] = -9.80665 + 65.351e-3 * y[1]**2 * math.exp(-10.53e-5*y[0])
    return F

x = 0.0
xStop = 10.0
y = numpy.array([9000, 0.0])
h = 0.5
freq = 1
X,Y = run_kut5.integrate(F,x,y,xStop,h,1.0e-2)
printSoln.printSoln(X,Y,freq)
raw_input("\nPress return to exit")
