## example 6_11
import math
import gaussQuad

def f(x): return (math.sin(x)/x)**2

a = 0.0; b = math.pi;
Iexact = 1.41815
for m in range(2,12):
    I = gaussQuad.gaussQuad(f,a,b,m)
    if abs(I - Iexact) < 0.00001:
        print "Number of nodes =",m
        print "Integral =", gaussQuad.gaussQuad(f,a,b,m)
        break
raw_input("\nPress return to exit")
