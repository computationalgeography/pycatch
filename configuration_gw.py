import pathlib

# define number of hourly timesteps to run
numberOfTimeSteps = 10968

# folder with input files (maps, timeseries)
#inputFolder = "inputs"
inputFolder = "../switzerland/40m_small_area"

# select maps as input parameters and initial values or uniform values over the area
mapsAsInput = False

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
#timeStepsToReportAll = list(range(100, numberOfTimeSteps + 1, 100))

# used for discharge only
timeStepsToReportRqs = list(range(1, numberOfTimeSteps + 1, 1))
#timeStepsToReportRqs = list(range(100, numberOfTimeSteps + 1, 100))

# definition for components were a subset of timesteps should be reported
timeStepsToReportSome = list(range(1, numberOfTimeSteps + 1, 1))
#timeStepsToReportSome = list(range(100, numberOfTimeSteps + 1, 100))

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
#setOfVariablesToReport = 'full'
setOfVariablesToReport = 'filtering'

#with_shading = True
with_shading = False

if with_shading is False:
  print("TEMPORARY SHADING SETTING IN CONFIGURATION.py")
  #fractionReceivedValue = 1.0
  #fractionReceivedFlatSurfaceValue = 1.0
  fractionReceivedValue = 0.0
  fractionReceivedFlatSurfaceValue = 0.0 

################
# model inputs #
################

#########################
# always needed as maps #
#########################

# set clone
cloneString = str(pathlib.Path(inputFolder, "mergeClone.map"))

# digital elevation model (m)
dem = str(pathlib.Path(inputFolder, "mergeDem.map"))

# ldd map
lddMap = str(pathlib.Path(inputFolder, "mergeldd.map"))

# report locations, i.e. outflow points, for instance, at the outlet
locations = str(pathlib.Path(inputFolder, "mergeOutFlowsNominal.map"))


#####################
# timeseries inputs #
#####################

# meteorology 
# same value across area (timseries has two columns) except rainfall (two areas, three columns)
rainfallFluxDetermTimeSeries = str(pathlib.Path(inputFolder, "rainfallFluxTwoCatchsJulAugSep0506.tss"))
airTemperatureDetermString = str(pathlib.Path(inputFolder, "airTemperatureArnaJulAugSep0506.tss"))
relativeHumidityDetermString = str(pathlib.Path(inputFolder, "relativeHumidityArnasJulAugSep0506.tss"))
incomingShortwaveRadiationFlatSurfaceString = str(pathlib.Path(inputFolder, "incomingShortwaveRadiationArnasJulAugSep0506.tss"))
windVelocityDetermString = str(pathlib.Path(inputFolder, "windVelocityArnasJulAugSep0506.tss"))

######
# inputs as maps or as uniform (constant) value over the area
######

if mapsAsInput:

  # forest (0) or no forest (1), only used when swapCatchments is True
  forestNoForest = str(pathlib.Path(inputFolder, "mergeForestNoForest.map"))
  areas = str(pathlib.Path(inputFolder, "mergeForestNoForest.map"))

  # meteorology
  # areas linked to rainfallFluxDetermTimeSeries (area code 1 and 2)
  rainfallFluxDetermTimeSeriesAreas = str(pathlib.Path(inputFolder, "mergeArnasSansaNominal.map"))
  elevationAboveSeaLevelOfMeteoStationValue = 900.0

  # interception
  maximumInterceptionCapacityValue = 0.0002
  leafAreaIndexValue = str(pathlib.Path(inputFolder, "mergeVegLAIFS.map"))

  # surface storage
  maxSurfaceStoreValue = 0.001

  # infiltration
  ksatValue = 0.0163
  initialSoilMoistureFractionFromDiskValue = str(pathlib.Path(inputFolder, "mergeFieldCapacityFractionFS.map"))
  soilPorosityFractionValue = str(pathlib.Path(inputFolder, "mergePorosityFractionFS.map"))

  # regolith geometry
  regolithThicknessHomogeneousValue = 1.0

  # groundwater layer geometry
  groundWaterLayerThicknessHomogeneousValue = 5.0
  #groundWaterLayerThicknessHomogeneousValue = 0.01

  # location of the stream, used to adjust regolith thickness there
  streamValue = str(pathlib.Path(inputFolder, "mergeStream.map"))

  # subsurface water 
  #saturatedConductivityMetrePerDayValue = 37.0
  #limitingPointFractionValue = str(pathlib.Path(inputFolder, "mergeLimitingPointFractionFS.map"))
  #mergeWiltingPointFractionFSValue = str(pathlib.Path(inputFolder, "mergeWiltingPointFractionFS.map"))
  #fieldCapacityFractionValue = str(pathlib.Path(inputFolder, "mergeFieldCapacityFractionFS.map"))
  #saturatedConductivityMetrePerDayValue = 0.008    # sand 8, silty sand 1, silt 0.008, clay much lower
  #saturatedConductivityMetrePerDayValue = 40.0
  saturatedConductivityMetrePerDayValue = 1.0
  limitingPointFractionValue = 0.2
  mergeWiltingPointFractionFSValue = 0.1
  fieldCapacityFractionValue = 0.4

  # evapotranspiration
  # penman
  multiplierMaxStomatalConductanceValue = 1.0
  albedoValue = str(pathlib.Path(inputFolder, "mergeVegAlbedoFS.map"))
  maxStomatalConductanceValue = str(pathlib.Path(inputFolder, "mergeVegStomatalFS.map"))
  vegetationHeightValue = str(pathlib.Path(inputFolder, "mergeVegHeightFS.map"))

else:
  
  # forest (0) or no forest (1), only used when swapCatchments is True
  forestNoForest = 0
  areas = 0

  # meteorology
  # areas linked to rainfallFluxDetermTimeSeries (area code 1 and 2)
  rainfallFluxDetermTimeSeriesAreas = 1
  elevationAboveSeaLevelOfMeteoStationValue = 900.0

  # interception
  maximumInterceptionCapacityValue = 0.0002
  leafAreaIndexValue = 3.0

  # surface storage
  maxSurfaceStoreValue = 0.001

  # infiltration
  ksatValue = 0.0163
  initialSoilMoistureFractionFromDiskValue = 0.327306
  soilPorosityFractionValue = 0.5

  # regolith geometry
  regolithThicknessHomogeneousValue = 1.0

  # groundwater layer geometry
  groundWaterLayerThicknessHomogeneousValue = 5.0
  #groundWaterLayerThicknessHomogeneousValue = 0.01

  # subsurface water 
  #saturatedConductivityMetrePerDayValue = 37.0
  #limitingPointFractionValue = 0.276293
  #mergeWiltingPointFractionFSValue = 0.1
  #fieldCapacityFractionValue = 0.327306
  # simple values for debugging
  #saturatedConductivityMetrePerDayValue = 0.008
  #saturatedConductivityMetrePerDayValue = 40.0  
  saturatedConductivityMetrePerDayValue = 1.0    # gravel high, sand 8, silty sand 1, silt 0.008, clay much lower
  limitingPointFractionValue = 0.2
  mergeWiltingPointFractionFSValue = 0.1
  fieldCapacityFractionValue = 0.4

  # evapotranspiration
  # penman
  multiplierMaxStomatalConductanceValue = 1.0
  albedoValue = 0.27
  maxStomatalConductanceValue = 0.0067
  vegetationHeightValue = 1.3


# real time of first time step, duration of time step
# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING
print("# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING")
print("# IMPORTANT NOTE: ALSO ONE NEEDS TO ENTER GEOG COORDINATES")
startTimeYearValue = 2005
startTimeMonthValue = 7
startTimeDayValue = 1
timeStepDurationHoursFloatingPointValue = 1.0  # only tested for one hour!!!!


# lat long for shading (solar radiation)
latitudeOfCatchment = 52.12833333
longitudeOfCatchment = 5.19861111
timeZone = "Europe/Madrid"

# calculate upstream totals (with accuflux) in subsurfacewateronelayer module
# and interceptionuptomaxstore module
# at least needed for some reports and possibly for budget checks (if one
# needs these) and almost certainly for report as numpy in subsurfacewateronelayer and
# interceptionuptomaxstore modules (not tested). For normal functioning of the model False is OK
calculateUpstreamTotals = False

#
# Reporting for the model components
#
if setOfVariablesToReport == 'full':
  interception_report_rasters = ["Vo", "Vi", "Vgf", "Vms"] # reports of totals (Vot) only make sense if calculateUpstreamTotals is True
  evapotrans_report_rasters = ["Ep", "Epc"]
  infiltration_report_rasters = ["Ii", "Ij", "Is", "Iks"]
  runoff_report_rasters = ["Rq"]
  shading_report_rasters = ["Mfs", "Msc", "Msh"]
  subsurface_report_rasters = ["Gs", "Go", "Gppa", "Gkun"]   # reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
  subsurface_report_rasters_gw = ["GGs", "GGgwd", "GGpcr", "GGpcrf", "GGdos", "GGkun"]   # reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
  surfacestore_report_rasters = ["Ss", "Sc"]
  randomparameters_report_rasters = ["RPic", "RPks", "RPrt", "RPsc", "RPmm"]
  exchange_report_rasters = ["Xrc"]
elif setOfVariablesToReport == 'filtering':
  interception_report_rasters = []  # reports of totals (Vot) only make sense if calculateUpstreamTotals is True
  evapotrans_report_rasters = []
  infiltration_report_rasters = []
  runoff_report_rasters = ["Rq"]
  #shading_report_rasters = []
  shading_report_rasters = ["Mfs", "Msc", "Msh"]
  subsurface_report_rasters = ["Gs", "Gppa", "Gdos", "Gkun"]   # reports of totals (accuflux) only make sense if calculateUpstreamTotals is True
  subsurface_report_rasters_gw = ["GGs", "GGgwd", "GGpcr", "GGpcrf", "GGdos", "GGkun"]   # reports of totals (accuflux) only make sense if calculateUpstreamTotals is True
  surfacestore_report_rasters = []
  #randomparameters_report_rasters = ["RPic", "RPks", "RPrt", "RPsc", "RPmm"]
  randomparameters_report_rasters = []
  #exchange_report_rasters = ["Xrc", "Xra"]
  exchange_report_rasters = []
