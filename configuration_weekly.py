import pathlib

# use for other runs
#numberOfTimeSteps=1500000   # long run (for hysteresis)
numberOfTimeSteps=5200   # test run

# option to fix both the regolith and the vegetation, not typically used
# in normal simulations
fixedStates=False

# option to call the methods that change the geomorphology
# this is typically on
changeGeomorphology=True

# number of realizations
nrOfSamples = 1

# calculation of early warning signals
intervalForStatsCalculated=100
#intervalForStatsCalculated=1

variances=False

# option to do data assimilation, always False, not implemented
filtering=False

# option to define some parameters as realizations from random
# functions, typically False
createRealizations=False

theDurationOfRainstorm=2.0

# option to define the variables to report to disc
# these are defined in the modules, it is typically either 'full' or 'filtering'
setOfVariablesToReport = 'filtering'
