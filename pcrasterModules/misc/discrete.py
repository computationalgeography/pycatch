import generalfunctions
from pcraster import *

resultA=generalfunctions.discreteProbabilityDistributionPerArea('discrete.tbl','discrete.map')
resultC=generalfunctions.discreteProbabilityDistributionPerCell('discrete.tbl','discrete.map')

report(resultA,'outA.map')
report(resultC,'outC.map')

