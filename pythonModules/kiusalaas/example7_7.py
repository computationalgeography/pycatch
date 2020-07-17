## example7_7
import numpy
import run_kut4
import printSoln

def F(x,y):                     
    F = numpy.zeros((2),dtype=float) 
    F[0] = y[1]
    F[1] = -4.75*y[0] - 10.0*y[1]
    return F
  
x = 0.0
xStop = 10.0
y = numpy.array([-9.0, 0.0])
h = 0.1
freq = 0

X,Y = run_kut4.integrate(F,x,y,xStop,h)
printSoln.printSoln(X,Y,freq)
raw_input("\nPress return to exit")
