## example6_7
import math
import romberg

def f(x): return 2.0*(x**2)*math.cos(x**2)

I,n = romberg.romberg(f,0,math.sqrt(math.pi))
print "Integral =",I
print "numEvals =",n
raw_input("\nPress return to exit")
