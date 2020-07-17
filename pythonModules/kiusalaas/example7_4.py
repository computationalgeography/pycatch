## example7_4
import numpy
import printSoln
import run_kut4

def F(x,y):
    F = numpy.zeros((2),dtype=float)
    F[0] = y[1]
    F[1] = -0.1*y[1] - x
    return F

x = 0.0                 # Start of integration
xStop = 2.0             # End of integration
y = numpy.array([0.0, 1.0])   # Initial values of {y}
h = 0.25                # Step size
freq = 1                # Printout frequency

X,Y = run_kut4.integrate(F,x,y,xStop,h)
printSoln.printSoln(X,Y,freq)
raw_input("Press return to exit")

