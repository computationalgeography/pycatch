#!/usr/bin/python
## example3_10
import numpy

xData = numpy.array([-0.04,0.93,1.95,2.90,3.83,5.0,      \
                5.98,7.05,8.21,9.08,10.09])
yData = numpy.array([-8.66,-6.44,-4.36,-3.27,-0.88,0.87, \
                3.31,4.63,6.19,7.4,8.85])

while 1:
    try:
        m = eval(raw_input("\nDegree of polynomial ==> "))
        coeff = numpy.polyfit(xData,yData,m)
        print "Coefficients are:\n",coeff
        ypred = numpy.polyval(coeff,xData)       
        delta = yData - ypred
        print "Std. deviation =",delta.std()
    except SyntaxError: break
    raw_input("Finished. Press return to exit")
