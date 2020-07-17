## example7_5
import numpy
import run_kut4
import printSoln
import math

def F(x,y):                     
    F = numpy.zeros((1),dtype=float) 
    F[0] = 3.0*y[0] - 4.0*math.exp(-x)
    return F
  
x = 0.0           # Start of integration
xStop = 10.0      # End of integration
y = numpy.array([1.0])  # Initial values of {y}
h = 0.1           # Step size
freq = 10         # Printout frequency

X,Y = run_kut4.integrate(F,x,y,xStop,h)
printSoln.printSoln(X,Y,freq)
raw_input("\nPress return to exit")
