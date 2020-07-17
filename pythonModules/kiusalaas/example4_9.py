#!/usr/bin/python
## example4_9
import numpy
import math
import newtonRaphson2

def f(x):
    f = numpy.zeros((len(x)),dtype=float)
    f[0] = math.sin(x[0]) + x[1]**2 + math.log(x[2]) - 7.0
    f[1] = 3.0*x[0] + 2.0**x[1] - x[2]**3 + 1.0
    f[2] = x[0] + x[1] + x[2] - 5.0
    return f

x = numpy.array([1.0, 1.0, 1.0])
print newtonRaphson2.newtonRaphson2(f,x)
raw_input ("\nPress return to exit")
