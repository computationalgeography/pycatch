import string
import os, shutil, numpy
import matplotlib.pyplot as plt
from PCRaster import *
from PCRaster.Framework import *

ndof = 3
x = numpy.random.standard_t(ndof,size=200)
print x

x.sort()

ecdf=stats.rankdata(x) / len(x)

