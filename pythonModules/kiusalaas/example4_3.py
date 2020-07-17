#!/usr/bin/python
## example4_3

import math
import rootsearch
import  bisect

def f(x): return x - math.tan(x)

a,b,dx = (0.0, 20.0, 0.01)
print "The roots are:"
while 1:
    x1,x2 = rootsearch.rootsearch(f,a,b,dx)
    if x1 != None:
        a = x2
        root = bisect.bisect(f,x1,x2,1)
        if root != None: print root
    else:
        print "\nDone"
        break
raw_input("Press return to exit")
