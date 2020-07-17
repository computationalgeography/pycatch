## example8_1
import numpy
import run_kut4
import brent
import printSoln

def initCond(u):  # Init. values of [y,y']; use 'u' if unknown
    return numpy.array([0.0, u])  

def r(u):         # Boundary condition residual--see Eq. (8.3)
    X,Y = run_kut4.integrate(F,xStart,initCond(u),xStop,h)
    y = Y[len(Y) - 1]
    r = y[0] - 1.0 
    return r

def F(x,y):       # First-order differential equations
    F = numpy.zeros((2),dtype=float)
    F[0] = y[1]
    F[1] = -3.0*y[0]*y[1]
    return F

xStart = 0.0        # Start of integration
xStop = 2.0         # End of integration
u1 = 1.0            # 1st trial value of unknown init. cond.
u2 = 2.0            # 2nd trial value of unknown init. cond.
h = 0.1             # Step size
freq = 2            # Printout frequency
u = brent.brent(r,u1,u2)  # Compute the correct initial condition
X,Y = run_kut4.integrate(F,xStart,initCond(u),xStop,h)
printSoln.printSoln(X,Y,freq)
raw_input("\nPress return to exit")
