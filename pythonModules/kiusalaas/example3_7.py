#!/usr/bin/python
## example3_7
import numpy
import cubicSpline

xData = numpy.array([1,2,3,4,5],dtype=float)
yData = numpy.array([0,1,0,1,0],dtype=float)
k = cubicSpline.curvatures(xData,yData)
while 1:
    try: x = eval(raw_input("\nx ==> "))
    except SyntaxError: break
    print "y =",cubicSpline.evalSpline(xData,yData,k,x)
raw_input("Done. Press return to exit")
