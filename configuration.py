import pathlib

# define number of hourly timesteps to run
numberOfTimeSteps = 10968

# folder with input files (maps, timeseries)
inputFolder = "inputs"

with_shading = True

if with_shading == False:
  fractionReceivedValue = 1.0
  fractionReceivedFlatSurfaceValue = 1.0

################
# model inputs #
################

# general ########

# set clone
cloneString = str(pathlib.Path(inputFolder, "mergeClone.map"))

dem = str(pathlib.Path(inputFolder, 'mergeDem.map'))
# report locations, i.e. outflow points, for instance, at the outlet
locations = str(pathlib.Path(inputFolder, 'mergeOutFlowsNominal.map'))
# map with forest or no forest, only used when swapCatchments is True
forestNoForest = str(pathlib.Path(inputFolder, 'mergeForestNoForest.map'))
areas = str(pathlib.Path(inputFolder, 'mergeForestNoForest.map'))


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
streamValue = str(pathlib.Path(inputFolder, 'mergeStream.map'))


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
lddMap = str(pathlib.Path(inputFolder, 'mergeldd.map'))


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
