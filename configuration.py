import pathlib

# define number of hourly timesteps to run
numberOfTimeSteps = 10968

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

# switch to report for locations as small numpy files
# mainly used for particle filtering
doReportComponentsDynamicAsNumpy = False

# when True, a particle filtering run is done
# first time users should have this False
filtering = False

# selects whether a single, given, value is used for a number of parameters
# or whether a realization for that parameters is drawn
# first time users will use a single, fixed value for these parameters, so
# use False and search on createRealizations in the script to see which
# parameters are defined like this
createRealizations = False

# switch to swap parameter values between two catchments
# first time users will need to set this to False
swapCatchments = False

# when True, one can read a set of parameters for all Monte Carlo realizations
# from disk (e.g. representing probability distributions from a calibration)
# first time users should have a False here
readDistributionOfParametersFromDisk = False

# switch to define which set of variables are reported
# either 'full' or 'filtering'. These are passed to the class of a component
# where it the variables that are reported can be defined, i.e. either full or filtering
# early users will always use full
setOfVariablesToReport = 'full'
# setOfVariablesToReport = 'filtering'


with_shading = True

if with_shading is False:
  fractionReceivedValue = 1.0
  fractionReceivedFlatSurfaceValue = 1.0

################
# model inputs #
################

# general ########

# set clone
cloneString = str(pathlib.Path(inputFolder, "mergeClone.map"))

dem = str(pathlib.Path(inputFolder, "mergeDem.map"))
# report locations, i.e. outflow points, for instance, at the outlet
locations = str(pathlib.Path(inputFolder, "mergeOutFlowsNominal.map"))
# map with forest or no forest, only used when swapCatchments is True
forestNoForest = str(pathlib.Path(inputFolder, "mergeForestNoForest.map"))
areas = str(pathlib.Path(inputFolder, "mergeForestNoForest.map"))


# meteorology #######

# observed precipitation
rainfallFluxDetermTimeSeries = str(pathlib.Path(inputFolder, "rainfallFluxTwoCatchsJulAugSep0506.tss"))
# areas linked to rainfallFluxDetermTimeSeries
rainfallFluxDetermTimeSeriesAreas = str(pathlib.Path(inputFolder, "mergeArnasSansaNominal.map"))

airTemperatureDetermString = str(pathlib.Path(inputFolder, "airTemperatureArnaJulAugSep0506.tss"))
relativeHumidityDetermString = str(pathlib.Path(inputFolder, "relativeHumidityArnasJulAugSep0506.tss"))
incomingShortwaveRadiationFlatSurfaceString = str(pathlib.Path(inputFolder, "incomingShortwaveRadiationArnasJulAugSep0506.tss"))
windVelocityDetermString = str(pathlib.Path(inputFolder, "windVelocityArnasJulAugSep0506.tss"))
elevationAboveSeaLevelOfMeteoStationValue = 900.0


# interception #######
maximumInterceptionCapacityValue = 0.0002
leafAreaIndexValue = str(pathlib.Path(inputFolder, "mergeVegLAIFS.map"))


# surface storage ######
maxSurfaceStoreValue = 0.001


# infiltration #######

# green and ampt
ksatValue = 0.0163
initialSoilMoistureFractionFromDiskValue = str(pathlib.Path(inputFolder, "mergeFieldCapacityFractionFS.map"))
soilPorosityFractionValue = str(pathlib.Path(inputFolder, "mergePorosityFractionFS.map"))


# regolith geometry ########
regolithThicknessHomogeneousValue = 0.5

# location of the stream, used to adjust regolith thickness there
streamValue = str(pathlib.Path(inputFolder, "mergeStream.map"))


# 'groundwater' (saturated flow) ##########
saturatedConductivityMetrePerDayValue = 37.0
limitingPointFractionValue = str(pathlib.Path(inputFolder, "mergeLimitingPointFractionFS.map"))
mergeWiltingPointFractionFSValue = str(pathlib.Path(inputFolder, "mergeWiltingPointFractionFS.map"))
fieldCapacityFractionValue = str(pathlib.Path(inputFolder, "mergeFieldCapacityFractionFS.map"))

# evapotranspiration ###########

# penman
multiplierMaxStomatalConductanceValue = 1.0
albedoValue = str(pathlib.Path(inputFolder, "mergeVegAlbedoFS.map"))
maxStomatalConductanceValue = str(pathlib.Path(inputFolder, "mergeVegStomatalFS.map"))
vegetationHeightValue = str(pathlib.Path(inputFolder, "mergeVegHeightFS.map"))


# dem geometry ###########
lddMap = str(pathlib.Path(inputFolder, "mergeldd.map"))


# real time of first time step, duration of time step
# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING
print("# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING")
startTimeYearValue = 2005
startTimeMonthValue = 7
startTimeDayValue = 1
timeStepDurationHoursFloatingPointValue = 1.0  # only tested for one hour!!!!


# lat long for shading (solar radiation)
latitudeOfCatchment = 52.12833333
longitudeOfCatchment = 5.19861111
timeZone = "Europe/Madrid"

#
# Reporting for the model components
#
if setOfVariablesToReport == 'full':
  interception_report_rasters = ["Vot", "Vo", "Vi", "Vgf", "Vms", "Vot"]
  evapotrans_report_rasters = ["Ep", "Epc"]
  infiltration_report_rasters = ["Ii", "Ij", "Is", "Iks"]
  runoff_report_rasters = ["Rq", "Rqs"]
  shading_report_rasters = ["Mfs", "Msc", "Msh"]
  subsurface_report_rasters = ["Gs", "Go", "Gxt", "Got"]
  surfacestore_report_rasters = ["Ss", "Sc"]
  randomparameters_report_rasters = ["RPic", "RPks", "RPrt", "RPsc", "RPmm"]
  exchange_report_rasters = ["Xrc"]
elif setOfVariablesToReport == 'filtering':
  interception_report_rasters = []
  evapotrans_report_rasters = []
  infiltration_report_rasters = []
  runoff_report_rasters = []
  shading_report_rasters = []
  subsurface_report_rasters = []
  surfacestore_report_rasters = []
  randomparameters_report_rasters = ["RPic", "RPks", "RPrt", "RPsc", "RPmm"]
  exchange_report_rasters = ["Xrc", "Xra"]
