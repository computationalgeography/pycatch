#!/usr/bin/python
## example3_4
import numpy
import newtonPoly
import math

xData = numpy.array([0.15,2.3,3.15,4.85,6.25,7.95])
yData = numpy.array([4.79867,4.49013,4.2243,3.47313,2.66674,1.51909])


a = newtonPoly.coeffts(xData,yData)
print " x    yInterp   yExact"
print "-----------------------"
for x in numpy.arange(0.0,8.1,0.5):
    y = newtonPoly.evalPoly(a,xData,x)
    yExact = 4.8*math.cos(math.pi*x/20.0)
    print "%3.1f %9.5f %9.5f"% (x,y,yExact)
raw_input("\nPress return to exit")


