import pathlib
import pcraster as pcr

# define number of hourly timesteps to run
numberOfTimeSteps = 10986

# folder with input files (maps, timeseries)
inputFolder = "inputs"


################
# model inputs #
################

# general ########

# set clone
cloneString = str(pathlib.Path(inputFolder, "mergeClone.map"))
pcr.setclone(cloneString)

dem = pcr.scalar(str(pathlib.Path(inputFolder, 'mergeDem.map')))
# report locations, i.e. outflow points, for instance, at the outlet
locations = pcr.nominal(str(pathlib.Path(inputFolder, 'mergeOutFlowsNominal.map')))
# map with forest or no forest, only used when swapCatchments is True
forestNoForest = pcr.boolean(str(pathlib.Path(inputFolder, 'mergeForestNoForest.map')))
areas = str(pathlib.Path(inputFolder, 'mergeForestNoForest.map'))


# meteorology #######

# observed precipitation
rainfallFluxDetermTimeSeries = str(pathlib.Path(inputFolder, "rainfallFluxTwoCatchsJulAugSep0506.tss"))
# areas linked to rainfallFluxDetermTimeSeries
rainfallFluxDetermTimeSeriesAreas = pcr.nominal(str(pathlib.Path(inputFolder, "mergeArnasSansaNominal.map")))

airTemperatureDetermString = str(pathlib.Path(inputFolder, "airTemperatureArnaJulAugSep0506.tss"))
relativeHumidityDetermString = str(pathlib.Path(inputFolder, "relativeHumidityArnasJulAugSep0506.tss"))
incomingShortwaveRadiationFlatSurfaceString = str(pathlib.Path(inputFolder, "incomingShortwaveRadiationArnasJulAugSep0506.tss"))
windVelocityDetermString = str(pathlib.Path(inputFolder, "windVelocityArnasJulAugSep0506.tss"))
elevationAboveSeaLevelOfMeteoStationValue = 900.0


# interception #######
maximumInterceptionCapacityValue = pcr.scalar(0.0002)
leafAreaIndexValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeVegLAIFS.map")))


# surface storage ######
maxSurfaceStoreValue = pcr.scalar(0.001)


# infiltration #######

# green and ampt
ksatValue = pcr.scalar(0.0163)
initialSoilMoistureFractionFromDiskValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeFieldCapacityFractionFS.map")))
soilPorosityFractionValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergePorosityFractionFS.map")))


# regolith geometry ########
regolithThicknessHomogeneousValue = pcr.scalar(0.5)

# location of the stream, used to adjust regolith thickness there
streamValue = pcr.boolean(str(pathlib.Path(inputFolder, 'mergeStream.map')))


# 'groundwater' (saturated flow) ##########
saturatedConductivityMetrePerDayValue = pcr.scalar(37.0)
limitingPointFractionValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeLimitingPointFractionFS.map")))
mergeWiltingPointFractionFSValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeWiltingPointFractionFS.map")))
fieldCapacityFractionValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeFieldCapacityFractionFS.map")))

# evapotranspiration ###########

# penman
multiplierMaxStomatalConductanceValue = pcr.scalar(1.0)
albedoValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeVegAlbedoFS.map")))
maxStomatalConductanceValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeVegStomatalFS.map")))
vegetationHeightValue = pcr.scalar(str(pathlib.Path(inputFolder, "mergeVegHeightFS.map")))


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
