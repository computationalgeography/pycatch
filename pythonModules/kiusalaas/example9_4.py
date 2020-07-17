## example9_4
import numpy
import math

s = numpy.array([[-30.0,  10.0,  20.0], \
           [ 10.0,  40.0, -50.0], \
           [ 20.0, -50.0, -10.0]])
v = numpy.array([1.0, 0.0, 0.0])
for i in range(100):
    vOld = v.copy()
    z = numpy.dot(s,v)
    zMag = math.sqrt(numpy.dot(z,z))
    v = z/zMag
    if numpy.dot(vOld,v) < 0.0:   
        sign = -1.0
        v = -v
    else: sign = 1.0
    if math.sqrt(numpy.dot(vOld - v,vOld - v)) < 1.0e-6: break
lam = sign*zMag
print "Number of iterations =",i
print "Eigenvalue =",lam
raw_input("Press return to exit")
