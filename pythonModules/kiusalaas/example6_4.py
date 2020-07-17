## example6_4
import math
import trapezoid

def f(x): return math.sqrt(x)*math.cos(x)

Iold = 0.0
for k in range(1,21):
    Inew = trapezoid.trapezoid(f,0.0,math.pi,Iold,k)
    if (k > 1) and (abs(Inew - Iold)) < 1.0e-6: break
    Iold = Inew
print "Integral =",Inew
print "nPanels =",2**(k-1)
raw_input("\nPress return to exit")




