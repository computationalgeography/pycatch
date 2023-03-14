import pathlib

# define number of hourly timesteps to run
numberOfTimeSteps = 360

# folder with input files (maps, timeseries)
inputFolder = "inputs"

# Define the number of Monte Carlo samples or particles
# first time users will use 1 and results for that realization are written to
# the folder '1'
nrOfSamples = 1

# when classes of components are initialized, we pass a list with the time steps
# that are reported. These are defined here. In principle for each component
# a different set of time steps can be reported, by just passing another list
# but this version uses three different ones

# definition for components were all timesteps should be reported
timeStepsToReportAll = list(range(1, numberOfTimeSteps + 1, 1))

# used for discharge only
timeStepsToReportRqs = list(range(20, numberOfTimeSteps + 1, 20))

# definition for components were a subset of timesteps should be reported
timeStepsToReportSome = list(range(100, numberOfTimeSteps + 1, 100))

################
# model inputs #
################

# general ########

# set clone
cloneString = str(pathlib.Path(inputFolder, "mergeClone.map"))

# dem
dem = str(pathlib.Path(inputFolder, "mergeDem.map"))

# report locations, i.e. outflow points, for instance, at the outlet
locations = str(pathlib.Path(inputFolder, "mergeOutFlowsNominal.map"))

# dem geometry ###########
lddMap = str(pathlib.Path(inputFolder, "mergeldd.map"))


# real time of first time step, duration of time step
# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING
print("# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING")
startTimeYearValue = 2005
startTimeMonthValue = 7
startTimeDayValue = 1
timeStepDurationHoursFloatingPointValue = 10.0/3600.0 
