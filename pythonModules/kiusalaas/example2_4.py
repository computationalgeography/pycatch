## example2_4
import numpy
from gaussElimin import *

def vandermode(v):
    n = len(v)
    a = numpy.zeros((n,n),dtype=float)
    for j in range(n):
        a[:,j] = v**(n-j-1)
    return a

v = numpy.array([1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
b = numpy.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
a = vandermode(v)
aOrig = a.copy()    # Save original matrix
bOrig = b.copy()    # and the constant vector

x = gaussElimin(a,b)
det = numpy.prod(numpy.diag(a))
print 'x =\n',x
print '\ndet =',det
print '\nCheck result: [a]{x} - b =\n', \
      numpy.dot(aOrig,x) - bOrig
