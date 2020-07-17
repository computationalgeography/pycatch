## module conjGrad
''' x, numIter = conjGrad(Av,x,b,tol=1.0e-9)
    Conjugate gradient method for solving [A]{x} = {b}.
    The matrix [A] should be sparse. User must supply
    the function Av(v) that returns the vector [A]{v}.
'''    
import numpy
import math

def conjGrad(Av,x,b,tol=1.0e-9):
    n = len(b)
    r = b - Av(x)
    s = r.copy()
    for i in range(n):
        u = Av(s)
        alpha = numpy.dot(s,r)/numpy.dot(s,u)
        x = x + alpha*s
        r = b - Av(x)
        if(math.sqrt(numpy.dot(r,r))) < tol:
            break
        else:
            beta = -numpy.dot(r,u)/numpy.dot(s,u)
            s = r + beta*s
    return x,i


    
