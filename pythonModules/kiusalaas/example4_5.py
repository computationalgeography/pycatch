#!/usr/bin/python
## example4_5
import math
import brent

def f(x): return x*abs(math.cos(x)) - 1.0

print "root =",brent.brent(f,0.0,4.0)
raw_input("Press return to exit")
