## module householder
''' d,c = householder(a).
    Householder similarity transformation of matrix [a] to 
    the tridiagonal form [c\d\c].

    p = computeP(a).
    Computes the acccumulated transformation matrix [p]
    after calling householder(a).
'''    
import numpy
import math


def householder(a): 
    n = len(a)
    for k in range(n-2):
        u = a[k+1:n,k]
        uMag = math.sqrt(numpy.dot(u,u))
        if u[0] < 0.0: uMag = -uMag
        u[0] = u[0] + uMag
        h = numpy.dot(u,u)/2.0
        v = numpy.dot(a[k+1:n,k+1:n],u)/h
        g = numpy.dot(u,v)/(2.0*h)
        v = v - g*u
        a[k+1:n,k+1:n] = a[k+1:n,k+1:n] - numpy.outer(v,u) \
                         - numpy.outer(u,v)
        a[k,k+1] = -uMag
    return numpy.diagonal(a),numpy.diagonal(a,1)

def computeP(a): 
    n = len(a)
    p = numpy.identity(n)*1.0
    for k in range(n-2):
        u = a[k+1:n,k]
        h = numpy.dot(u,u)/2.0
        v = numpy.dot(p[1:n,k+1:n],u)/h           
        p[1:n,k+1:n] = p[1:n,k+1:n] - numpy.outer(v,u)
    return p
      
                

