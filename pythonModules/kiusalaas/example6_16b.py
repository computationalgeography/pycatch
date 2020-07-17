## problem6_3_13
import numpy
import math
import triangleQuad
    
def f(x,y):
    return (x**2 + y**2)/2.0        \
           -(x**3 - 3.0*x*y**2)/6.0 \
           -2.0/3.0

xCorner = numpy.array([-1.0, -1.0, 2.0])
yCorner = numpy.array([math.sqrt(3.0), -math.sqrt(3.0), 0.0])
print "Integral =",triangleQuad.triangleQuad(f,xCorner,yCorner)
raw_input("Press return to  exit")
