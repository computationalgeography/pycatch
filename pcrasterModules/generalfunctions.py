import pcraster as pcr
import pcraster.framework as pcrfw
import scipy.stats

import numpy, os, operator, glob, subprocess
##
#import rpy2.robjects.numpy2ri
# rpy2.robjects.numpy2ri.activate() # nodig bij nieuwere rpy2 versie
#import rpy2.robjects as robjects
##
from collections import deque
#from PCRaster.NumPy import *
import random
#from osgeo import gdal

# time in hours

###########################################
# retrieving data from maps, writing maps #
###########################################

# at first time step and sample, create self.a=numpy.empty(1000*100*2*1).reshape(1000,100,2,1)
# of als object self.rainAsNumpy=class(filenameToBeWritten,locs,nrsamples,nrtimesteps,rows,cols)
# en dan in de dynamic self.rainAsNumpy.report(Rain,currentSample,currentTimestep)
# of in premcloop
# c is a 2 value array
# >>> b[999,99,:,0]=c
# as soon as all timesteps and samples have been done, write (or at last time step and last sample)
# self.a=writeToNumpy(self.a,locs,variable,nameString,currentSample,currentTimestep,nrOfSamples,nrOfTimesteps)

def getCellValue(Map, Row, Column):
  Value, Valid = pcr.cellvalue(Map, Row, Column)
  if Valid:
    return Value
  else:
    print('missing value in input of getCellValue')

def printErrorMessageIfACellContainsTrue(booleanMap, errorMessage):
  scalarMap = pcr.cover(pcr.scalar(booleanMap), 0)
  cellContainsTrueMap = pcr.boolean(pcr.mapmaximum(scalarMap))
  cellContainsTrue = getCellValue(cellContainsTrueMap, 1, 1)
  if cellContainsTrue > 0.5:
    print(errorMessage)

def getCellValueAtBooleanLocation(location, map):
  # map can be any type, return value always float
  valueMap = pcr.mapmaximum(pcr.ifthen(location, pcr.scalar(map)))
  # to get rid of bug that gives very low value in case of missing value on location
  valueMap = pcr.ifthen(pcr.pcrgt(valueMap, -1e10), valueMap)
  value = getCellValue(valueMap, 1, 1)
  return value

def getCellValueAtBooleanLocationReturnMVOrNot(location, map):
  # map can be any type, return value always float
  valueMap = pcr.mapmaximum(pcr.ifthen(location, pcr.scalar(map)))
  # to get rid of bug that gives very low value in case of missing value on location
  valueMap = pcr.ifthen(pcr.pcrgt(valueMap, -1e10), valueMap)
  value, valid = pcr.cellvalue(valueMap, 1, 1)
  return value, valid

def printCellValue(self, mapVariable, variableNameToPrint, unit, row, column):
  cellValue = getCellValue(mapVariable, row, column)
  print(variableNameToPrint + ' (' + unit + ') at row ' + str(row) + ', column: ' + str(column) + ' is: ' + str(cellValue))

def onePeriod(self, startTime, endTime, timeStepDuration, currentTimeStep):
  # this could be separated in two functions, one converting hours to
  # time steps, one creating the period
  time = float(currentTimeStep) * float(timeStepDuration)
  period = (time > startTime) & (time < endTime)
  return period

def returnCellValuesAtLocationsAsNumpyArray(locations, values):
  # output gives values at locations as numpy array
  # ordered by id's of locations, so lowest location number
  # is first in out
  locsAsNumpy = pcr.pcr2numpy(locations, 9999)
  valuesAsNumpy = pcr.pcr2numpy(values, 9999)
  locsOnly = locsAsNumpy[locsAsNumpy > 0.5]
  valuesOnly = valuesAsNumpy[locsAsNumpy > 0.5]
  all = numpy.dstack((locsOnly, valuesOnly))
  allZero = all[0]
  final = allZero[allZero[:, 0].argsort(), ]
  out = final[:, 1]
  return out

def writeNumpyArrayAsScalarPCRasterMap(numpyArray, fileName):
  driver = None
  cols = len(numpyArray)
  rows = 2
  driver = None
  dst_ds = None
  dst = None
  driver = gdal.GetDriverByName('RST')
  #dst_ds = driver.Create('piet.rst', cols, rows, 1, gdal.GDT_Float32 )
  tmpFile = fileName + 'piet.rst'
  dst_ds = driver.Create(tmpFile, cols, rows, 1, gdal.GDT_Float32 )
  numpyArrayTwoRows = numpy.array((numpyArray, numpyArray))
  dst_ds.GetRasterBand(1).WriteArray(numpyArrayTwoRows)

  for filename in glob.glob('piet.*'):
    os.remove(filename)

  driver = gdal.GetDriverByName('PCRaster')
  print(fileName)
  dst_ds_pcr = driver.CreateCopy( fileName, dst_ds, 1, ['Type=Float32'] )
  print('create copy gedaan ' + fileName)
  del dst_ds
  del dst_ds_pcr

def writeNumpyArrayAsScalarPCRasterMapAsc2Map(numpyArray, fileName):
  numpy.savetxt('tmp.asc', numpyArray)
  length = len(numpyArray)
  command = 'mapattr -S -R ' + str(length) + ' -s -C 1 tmp.clone'
  p = subprocess.call(command, shell=True)
  command = 'asc2map tmp.asc ' + fileName + ' -S --clone tmp.clone'
  p = subprocess.call(command, shell=True)
  os.remove('tmp.clone')
  os.remove('tmp.asc')

def reportLocations(locations, values, basename, sampleNumber, timeStep):
  # let op locations should have no mv's
  fileName = pcrfw.generateNameST(basename, sampleNumber, timeStep)
  numpyArray = returnCellValuesAtLocationsAsNumpyArray(locations, values)
  writeNumpyArrayAsScalarPCRasterMapAsc2Map(numpyArray, fileName)

def reportLocationsAsNumpyArray(locations, values, basename, sampleNumber, timeStep):
  # reports one file per realization and per time step
  fileName = pcrfw.generateNameST(basename, sampleNumber, timeStep)
  numpyArray = returnCellValuesAtLocationsAsNumpyArray(locations, pcr.spatial(values))
  numpyArrayAsMapWithOneRow = numpyArray.reshape(1, len(numpyArray))
  numpy.savetxt(fileName + '.numpy.txt', numpyArrayAsMapWithOneRow)


def reportLocationsAsNumpyArrayOneFilePerRealization(locations, values,
                       basename, sampleNumber, timeStep, endTimeStep):

  fileName = pcrfw.generateNameS(basename, sampleNumber) + '.numpy.txt'

  if timeStep == 1:
    theFile = file(fileName, 'w')
  else:
    theFile = file(fileName, 'a')

  numpyArray = returnCellValuesAtLocationsAsNumpyArray(locations, pcr.spatial(values))
  numpyArrayAsMapWithOneRow = numpyArray.reshape(1, len(numpyArray))
  numpy.savetxt(theFile, numpyArrayAsMapWithOneRow)

  theFile.close()


def reportAsNumpyArray(values, basename, sampleNumber, timeStep):
  fileName = pcrfw.generateNameST(basename, sampleNumber, timeStep)
  valuesAsNumpy = pcr.pcr2numpy(values, 9999)
  numpy.savetxt(fileName + '.numpy.txt', valuesAsNumpy)

def openFileAsNumpyArray(name):
  src_ds = gdal.Open(name)
  cols = src_ds.RasterXSize
  rows = src_ds.RasterYSize
  rasterBand = src_ds.GetRasterBand(1)
  numpyArray = numpy.array(rasterBand.ReadAsArray())
  return numpyArray

def openAsNumpyArray(basename, sampleNumber, timeStep):
  fileName = pcrfw.generateNameST(basename, sampleNumber, timeStep)
  numpyArray = openFileAsNumpyArray(fileName)
  return numpyArray

def openSamplesAndTimestepsAsNumpyArray(basename, samples, timesteps):
  t = 1
  output = []
  for timestep in timesteps:
    print('timestep ' + str(timestep) + ' done,',)
    allSamples = []
    for sample in samples:
      array = openAsNumpyArray(basename, sample, timestep)
      allSamples.append(array)
    output.append(allSamples)
  outputAsArray = numpy.array(output)
  return outputAsArray

def openSamplesAndTimestepsAsNumpyArraysAsNumpyArrayTimeThenSamples(basename, samples, timesteps):
  # this is the same (older) as openSamplesAndTimestepsAsNumpyArraysAsNumpyArray
  # but it loops for each time step over all samples, which appears to be slower
  t = 1
  output = []
  print('doing basename ', basename)
  for timestep in timesteps:
    allSamples = []
    for sample in samples:
      #fileName = pcrfw.generateNameST(basename,sample,timestep) + '.npy'
      # array=numpy.load(fileName)
      fileName = pcrfw.generateNameST(basename, sample, timestep) + '.numpy.txt'
      array = numpy.atleast_2d(numpy.loadtxt(fileName))
      allSamples.append(array)
    print(timestep)
    output.append(allSamples)
  outputAsArray = numpy.array(output)
  return outputAsArray

def convertTimeseriesOfMapFilesToNumpyArray(basename, samples, timesteps):
  t = 1
  output = []
  for timestep in timesteps:
    allSamples = []
    for sample in samples:
      fileName = pcrfw.generateNameST(basename, sample, timestep)
      valuesAsNumpy = pcr.pcr2numpy(pcr.spatial(fileName), 9999)
      # array=numpy.atleast_2d(numpy.loadtxt(fileName))
      allSamples.append(valuesAsNumpy)
    output.append(allSamples)
  outputAsArray = numpy.array(output)
  return outputAsArray

def openSamplesAndTimestepsAsNumpyArraysAsNumpyArray(basename, samples, timesteps):
  print('doing basename ', basename)
  done = 0
  firstFileName = pcrfw.generateNameST(basename, samples[0], timesteps[0]) + '.numpy.txt'
  array = numpy.atleast_2d(numpy.loadtxt(firstFileName))
  a = numpy.shape(array)
  b = (len(timesteps), len(samples))
  outputAsArray = numpy.ones(b + a)
  sampleIndex = 0
  for sample in samples:
    timestepIndex = 0
    for timestep in timesteps:
      fileName = pcrfw.generateNameST(basename, sample, timestep) + '.numpy.txt'
      array = numpy.atleast_2d(numpy.loadtxt(fileName))
      outputAsArray[timestepIndex, sampleIndex, ] = array
      done = done + 1
      timestepIndex += 1
    print(sample, done)
    sampleIndex += 1
  print('old', outputAsArray)
  return outputAsArray

def openSamplesAsNumpyArrays(basename, samples, timesteps):
  # opens for each realization a timeseries (stored with generalfunctions.reportLocationsAsNumpyArrayOneFilePerRealization)
  # and stores it in a multi dimensional numpy array and writes to disk
  print('doing basename ', basename)
  done = 0
  #firstFileName = pcrfw.generateNameST(basename,samples[0],timesteps[0]) + '.numpy.txt'
  firstFileName = pcrfw.generateNameS(basename, samples[0]) + '.numpy.txt'
  array = numpy.atleast_2d(numpy.loadtxt(firstFileName)[0])
  a = numpy.shape(array)
  b = (len(timesteps), len(samples))
  outputAsArray = numpy.ones(b + a)
  sampleIndex = 0
  for sample in samples:
    fileName = pcrfw.generateNameS(basename, sample) + '.numpy.txt'
    timeSeries = numpy.atleast_2d(numpy.loadtxt(fileName))
    timestepIndex = 0
    for timestep in timesteps:
      array = timeSeries[timestepIndex]
      outputAsArray[timestepIndex, sampleIndex, ] = array
      done = done + 1
      timestepIndex += 1
    print(sample, done)
    sampleIndex += 1
  print('new', outputAsArray)
  return outputAsArray


def createList(samples, timesteps):
  wholeList = []
  allSamples = samples[:]
  for timestep in timesteps:
    wholeList.append(allSamples)
  return wholeList

def test(basename, samples, timesteps):
  t = 1
  output = []
  print('doing basename ', basename)
  for timestep in timesteps:
    allSamples = []
    for sample in samples:
      print(timestep, sample)
      fileName = pcrfw.generateNameST(basename, sample, timestep) + '.numpy.txt'
      array = numpy.atleast_2d(numpy.loadtxt(fileName))
      allSamples.append(array)
    output.append(allSamples)
  outputAsArray = numpy.array(output)
  return outputAsArray

def testTwo(basename, samples, timesteps):
  print('doing basename ', basename)
  result = createList(samples, timesteps)
  t = 0
  for timestep in timesteps:
    s = 0
    for sample in samples:
      print(timestep, sample)
      fileName = pcrfw.generateNameST(basename, sample, timestep) + '.numpy.txt'
      theArray = numpy.atleast_2d(numpy.loadtxt(fileName))
      result[t][s] = theArray
      s = s + 1
    t = t + 1
  outputAsArray = numpy.array(result)
  return outputAsArray




###########################
#      MAP ALGEBRA        #
###########################

def selectACell(Map, XInNrCells, YInNrCells):
  xMap = pcr.nominal(pcr.xcoordinate(pcr.boolean(Map)) / pcr.celllength())
  yMap = pcr.nominal(pcr.ycoordinate(pcr.boolean(Map)) / pcr.celllength())
  location = pcr.pcrand((xMap == XInNrCells), (yMap == YInNrCells))
  return location

def mapeq(mapOne, mapTwo):
  mapOneScalar = pcr.scalar(mapOne)
  mapTwoScalar = pcr.scalar(mapTwo)
  difference = mapOneScalar - mapTwoScalar
  cellEqual = pcr.pcreq(difference, pcr.scalar(0))
  mapEqual = pcr.pcrgt(pcr.mapminimum(pcr.scalar(cellEqual)), pcr.scalar(0.5))
  return getCellValue(mapEqual, 1, 1)

def slopeToDownstreamNeighbour(dem, ldd):
  slopeToDownstreamNeighbour = (dem - pcr.downstream(ldd, dem)) / pcr.downstreamdist(ldd)
  return slopeToDownstreamNeighbour

def slopeToDownstreamNeighbourNotFlat(dem, ldd, minSlope):
  slopeToDownstreamNeighbourMap = slopeToDownstreamNeighbour(dem, ldd)
  lddArea = pcr.defined(ldd)
  minSlopeCover = pcr.ifthen(lddArea, pcr.scalar(minSlope))
  slopeToDownstreamNeighbourNotFlat = pcr.cover(pcr.max(minSlopeCover, slopeToDownstreamNeighbourMap), minSlopeCover)
  return slopeToDownstreamNeighbourNotFlat

def distancetodownstreamcell(Ldd):
  distanceToDownstreamCell = pcr.max(pcr.downstreamdist(Ldd), pcr.celllength())
  return distanceToDownstreamCell

def normalcorrelated(normalX, normalY, correlation):
  # returns realizations of two normal variables with
  # mean zero and var 1 having correlation of correlation
  # based on:
  # x=normal()
  # y=ax+b*normal()
  # correlation = a /  pcr.sqrt( pcr.sqr(a) + pcr.sqr(b) )
  x = pcr.scalar(normalX)
  y = (x + pcr.sqrt((1 / pcr.sqr(correlation)) - 1) * pcr.scalar(normalY)) * pcr.scalar(correlation)
  return x, y

def swapValuesOfTwoRegions(regions, values, doIt):
  # assigns the highest value found in region False to all cells
  # in region True, and vice versa
  # regions, a boolean map with two regions
  # values, a map of scalar data type
  if doIt:
    valueInRegionFalse = pcr.mapmaximum(pcr.ifthen(pcr.pcrnot(regions), values))
    valueInRegionTrue = pcr.mapmaximum(pcr.ifthen(regions, values))
    swapped = pcr.ifthenelse(regions, valueInRegionFalse, valueInRegionTrue)
    return swapped
  else:
    return values


##############################
# converting to numpy stuff  #
##############################

def createTimeSeriesList(timeSeriesFile):
  file = open(timeSeriesFile, 'r')
  piet = file.readlines()
  newList = []
  for line in piet:
    lineList = string.split(line)
    newList.append(lineList)
  file.close()
  return newList

def timeInputSparse(fileName):
  return os.path.exists(fileName)

def mapToColAsArray(name):
  """Selects values at row, col from raster name in Monte Carlo samples.

  name -- Name of raster.
  row -- Row index of cell to read.
  col -- Col index of cell to read.
  The returned array does not contain missing values so the size is maximimal
  the number of cells. It contains three columns, x, y, name
  x,y are given as xcoordinate and ycoordinate values

  Returned array has elements of type numpy.float32"""
  nrRows = clone().nrRows()
  nrCols = clone().nrCols()
  nrCells = nrRows * nrCols

  mask = numpy.zeros(nrCells).astype(numpy.bool_)
  arrayX = numpy.zeros(nrCells).astype(numpy.float32)
  arrayY = numpy.zeros(nrCells).astype(numpy.float32)
  arrayName = numpy.zeros(nrCells).astype(numpy.float32)

  xMap = pcr.xcoordinate(pcr.defined(name))
  yMap = pcr.ycoordinate(pcr.defined(name))

  # For each cell.
  c = 0
  while c < nrCells:
    arrayName[c], mask[c] = pcr.cellvalue(name, c + 1)
    arrayX[c], dummy = pcr.cellvalue(xMap, c + 1)
    arrayY[c], dummy = pcr.cellvalue(yMap, c + 1)
    c += 1

  arrayName = numpy.compress(mask, arrayName)
  arrayX = numpy.compress(mask, arrayX)
  arrayY = numpy.compress(mask, arrayY)
  mapAsColArray = numpy.column_stack((arrayX, arrayY, arrayName))
  return mapAsColArray

def addTimeColumnToMapAsColArray(mapAsColArray, time):
  b = numpy.insert(mapAsColArray, 2, float(time), axis=1)
  return b

def stackOfMapsToColAsArray(stackOfMapsAsList, currentTime):
  timeOfFirstMap = currentTime - len(stackOfMapsAsList) + 1
  t = timeOfFirstMap
  stackOfMapsAsColsList = []
  for map in stackOfMapsAsList:
    mapArray = mapToColAsArray(map)
    mapAsColArrayWithTime = addTimeColumnToMapAsColArray(mapArray, t)
    stackOfMapsAsColsList.append(mapAsColArrayWithTime)
    t =  t + 1
  array = numpy.concatenate(stackOfMapsAsColsList, axis=0)
  return array

def stackOfMapsToRDataFrame(stackOfMapsAsList, currentTime):
  colAsArray = stackOfMapsToColAsArray(stackOfMapsAsList, currentTime)
  dataFrame = convertStackOfMapsToRDataFrame(colAsArray)
  return dataFrame

def loadTimeseries(timeSeriesFileName):
  """Reads a PCRaster timeseries that should not contain a header

  timeSeriesFileName -- Name of timeseries.

  The returned numpy array has two dimensions
  First dimension: timesteps, second dimension values, where first item is timestep
  second item is first value column, third item is second value column, et.  """

  a = numpy.loadtxt(timeSeriesFileName)
  return a


##################################
#            links to R          #
##################################

def convertStackOfMapsToRDataFrame(mapAsColArray):
  # note z is time and v is variable
  robjects.r('''
  convertToDataFrame <- function(x) {
        frame=as.data.frame(x)
        colnames(frame)[1] <- "x"
        colnames(frame)[2] <- "y"
        colnames(frame)[3] <- "z"
        colnames(frame)[4] <- "v"
        frame
        }
        ''')
  convertToDataFrame = robjects.r['convertToDataFrame']
  a = convertToDataFrame(mapAsColArray)
  return a

def experimentalVariogramValues(stackOfMapsAsList, boundariesVector, space, savePlot, fileName, maxVarPlot):
  # returns distances and semivariances for steps defined by boundaries vector
  # note that length of returned vector equals number of intervals that is
  # available, thus, len(returnedList) can be smaller than len(boundariesVector)!
  # space (TRUE) -> spatial correlation
  # space (FALSE) -> temporal correlation
  stackOfMapsAsRDataFrame = stackOfMapsToRDataFrame(stackOfMapsAsList, 10)
  if space:
    robjects.r('''
    experimentalVariogramValues <- function(dataFrame,boundariesVector) {
        require("gstat")
        require("automap")
        #gstatDataFrame <- gstat(id = "v",formula = v~1, locations = ~x+y+z, data = dataFrame)
        # spatial only
        a = variogram(v ~ 1, ~ x + y + z, dataFrame, beta=0,tol.ver=0.1, boundaries=boundariesVector)
        #plotting seems not to work
        #data(dataFrame)
        #coordinates(dataFrame)=~x+y+z
        #variogram=autofitVariogram(v~1, dataFrame,model="Exp")
        #b=fit.variogram(a, vgm(1,"Exp",3))
        #pdf("test.pdf")
        #plot(a$dist,a$gamma)
        #dev.off()
        a
        }
        ''')
  else:
    robjects.r('''
    experimentalVariogramValues <- function(dataFrame,boundariesVector) {
        require("gstat")
        require("automap")
        #gstatDataFrame <- gstat(id = "v",formula = v~1, locations = ~x+y+z, data = dataFrame)
        # temporal only
        colnames(dataFrame)[1] <- "z"
        colnames(dataFrame)[3] <- "x"
        #gstatDataFrame <- gstat(id = "v",formula = v~1, locations = ~x+y+z, data = dataFrame)
        # normal
        #a = variogram(v ~ 1, ~ x + y + z, dataFrame, beta=0,tol.ver=0.1, alpha=90,tol.hor=0.1, boundaries=boundariesVector)
        # remove trend (universal kriging variogram)
        a = variogram(v ~ x, ~ x + y + z, dataFrame, beta=0,tol.ver=0.1, alpha=90,tol.hor=0.1, boundaries=boundariesVector)
        a
        }
        ''')
  experimentalVariogramValues = robjects.r['experimentalVariogramValues']
  boundariesVectorR = robjects.FloatVector(boundariesVector)
  expVariogram = experimentalVariogramValues(stackOfMapsAsRDataFrame, boundariesVectorR)
  if savePlot:
    robjects.r('''
    saveExperimentalVariogram <- function(experimentalVariogram,fileName,maxVarPlot) {
        require("gstat")
        png(fileName)
        plot(experimentalVariogram$dist,experimentalVariogram$gamma,ylim=c(0,maxVarPlot))
        dev.off()
        }
        ''')
    saveExperimentalVariogram = robjects.r['saveExperimentalVariogram']
    saveExperimentalVariogram(expVariogram, fileName, maxVarPlot)
  # return expVariogram.r['dist'][0], expVariogram.r['gamma'][0]
  return expVariogram[1], expVariogram[2]

def semvar(firstMap, secondMap):
  nrPairs = getCellValue( pcr.cover(pcr.maptotal(pcr.scalar(pcr.pcrand(pcr.defined(firstMap), pcr.defined(secondMap)))), 0), 1, 1)
  sumOfSquaredDiff = getCellValue( pcr.cover(pcr.maptotal(pcr.sqr(firstMap - secondMap) / 2.0), 0), 1, 1)
  return nrPairs, sumOfSquaredDiff

def experimentalVariogramValuesInTime(stackOfMapsAsList, bounds):
  nrPairsOfLags = [0.0] * len(bounds)
  print(nrPairsOfLags)
  sumOfSquaredDiffOfLags = [0.0] * len(bounds)
  sumOfDists = [0.0] * len(bounds)
  nMaps = len(stackOfMapsAsList)
  for i in range(0, nMaps):
    for j in range(i + 1, nMaps):
      dist = math.fabs(i - j)
      nrPairs, sumOfSquaredDiff = semvar(stackOfMapsAsList[i], stackOfMapsAsList[j])
      k = 0
      used = 0
      for bound in bounds:
        if (dist < bound) and (used == 0):
          nrPairsOfLags[k] = nrPairsOfLags[k] + nrPairs
          sumOfSquaredDiffOfLags[k] = sumOfSquaredDiffOfLags[k] + sumOfSquaredDiff
          sumOfDists[k] = sumOfDists[k] + nrPairs * float(dist)
          used = 1
        k = k + 1
  semvarList = map(operator.div, sumOfSquaredDiffOfLags, nrPairsOfLags)
  distList = map(operator.div, sumOfDists, nrPairsOfLags)
  return distList, semvarList

# jan=pcr.ifthen(pcr.pcrle(pcr.uniqueid(pcr.defined('jet00000.001')),100),pcr.scalar('jet00000.001'))
#test=[pcr.scalar('jet00000.001'),jan, pcr.scalar('jet00000.003'),pcr.scalar('jet00000.004')]
# a=experimentalVariogramValuesInTime(test,[1.5,7.0])

def semvarOfStackOfMapsInSpace(stackOfMapsAsList, lagX, lagY):
  # lagX is shift of cells to right (positive)
  # lagY is shift of cells up (positive)
  nrPairsTot = 0.0
  sumOfSquaredDiffTot = 0.0
  for theMap in stackOfMapsAsList:
    shiftedMap = pcr.shift(theMap, lagY, 0 - lagX)
    nrPairs, sumOfSquaredDiff = semvar(theMap, shiftedMap)
    nrPairsTot = nrPairsTot + nrPairs
    sumOfSquaredDiffTot = sumOfSquaredDiffTot + sumOfSquaredDiff
  return nrPairsTot, sumOfSquaredDiffTot

def createLagXAndLagYForBounds(bounds):
  lagXlagYDist = []
  maxLag = int(round(bounds[-1]))
  possibleLags = range(0, maxLag, 1)
  for i in possibleLags:
    for j in possibleLags:
      dist = math.pcr.sqrt((float(i)**2) + (float(j)**2))
      if (dist < bounds[-1]) and not ((i == 0) and (j == 0)):
        a = [ i, j, dist]
        lagXlagYDist.append(a)
  return lagXlagYDist

# b=createLagXAndLagYForBounds([0.2,4.0])
# print b


def experimentalVariogramValuesInSpace(stackOfMapsAsList, bounds):
  nrPairsOfLags = [0.0] * len(bounds)
  sumOfSquaredDiffOfLags = [0.0] * len(bounds)
  sumOfDists = [0.0] * len(bounds)
  nMaps = len(stackOfMapsAsList)
  lagXAndLagY = createLagXAndLagYForBounds(bounds)
  for i in lagXAndLagY:
    dist = i[2]
    nrPairs, sumOfSquaredDiff = semvarOfStackOfMapsInSpace(stackOfMapsAsList, i[0], i[1])
    k = 0
    used = 0
    for bound in bounds:
      if (dist < bound) and (used == 0):
        nrPairsOfLags[k] = nrPairsOfLags[k] + nrPairs
        sumOfSquaredDiffOfLags[k] = sumOfSquaredDiffOfLags[k] + sumOfSquaredDiff
        sumOfDists[k] = sumOfDists[k] + nrPairs * float(dist)
        used = 1
      k = k + 1
  semvarList = map(operator.div, sumOfSquaredDiffOfLags, nrPairsOfLags)
  distList = map(operator.div, sumOfDists, nrPairsOfLags)
  return distList, semvarList

# jan=pcr.ifthen(pcr.pcrle(pcr.uniqueid(pcr.defined('jet00000.001')),100),pcr.scalar('jet00000.001'))
##test=[pcr.scalar('jet00000.001'),jan, pcr.scalar('jet00000.003'),pcr.scalar('jet00000.004')]
#test=[pcr.scalar('jet00000.001'),pcr.scalar('jet00000.002'), pcr.scalar('jet00000.003'),pcr.scalar('jet00000.004')]
# a,b=experimentalVariogramValuesInSpace(test,[2.3,15,17.8])
# print a, b

def descriptiveStatistics(stackOfMapsAsRDataFrame):
  robjects.r('''
  descriptiveStatistics <- function(dataFrame) {
      mean <- mean(dataFrame$v)
      variance <- var(dataFrame$v)
      c(mean,variance)
      }
      ''')
  descriptiveStatistics = robjects.r['descriptiveStatistics']
  var = descriptiveStatistics(stackOfMapsAsRDataFrame)
  return var

###########################
#   some data management  #
###########################

def cornerMap(cloneMap):
  x = pcr.xcoordinate(cloneMap)
  y = pcr.ycoordinate(cloneMap)
  corner = pcr.pcrand(pcr.pcreq(x, pcr.mapminimum(x)), pcr.pcreq(y, pcr.mapmaximum(y)))
  return corner

def convertListOfValuesToListOfNonSpatialMaps(listOfValues, cloneMap):
  # non spatial, i.e. value is put in the corner with the largest x and y
  # return values is always scalar
  corner = cornerMap(cloneMap)
  listOfMaps = []
  for value in listOfValues:
    map = pcr.ifthen(corner, pcr.scalar(value))
    listOfMaps.append(map)
  return listOfMaps

def convertListOfNonSpatialMapsToListOfValues(listOfMaps):
  listOfValues = []
  clone = pcr.ifthenelse(pcr.defined(listOfMaps[0]), pcr.boolean(1), 1)
  x = pcr.xcoordinate(clone)
  y = pcr.ycoordinate(clone)
  corner = pcr.pcrand(pcr.pcreq(x, pcr.mapminimum(x)), pcr.pcreq(y, pcr.mapmaximum(y)))
  for map in listOfMaps:
    value = getCellValueAtBooleanLocation(corner, map)
    listOfValues.append(value)
  return listOfValues

# tests
# setclone("cloneSmall.map")
#a = pcr.scalar("norSmall.map")
#b = pcr.scalar("norSmall.map")*2
#c = pcr.scalar("norSmall.map")*3
#d = pcr.scalar("norSmall.map")*4
#e = pcr.scalar("norSmall.map")*5
# stackOfMapsAsList=[a,b,c,d,e]
# pcr.report(d,"testje")
#test = convertListOfNonSpatialMapsToListOfValues(stackOfMapsAsList)
#c = stackOfMapsToColAsArray(stackOfMapsAsList, 10)
#d = convertStackOfMapsToRDataFrame(c)
# boundVector=(1.5,2.5,30.5)
# dist,gamma=experimentalVariogramValues(d,boundVector,1,1,'pietje.pdf')
# print dist
# print gamma
# aListOfMaps=convertListOfValuesToListOfNonSpatialMaps(gamma,"cloneSmall.map")
# pcr.report(aListOfMaps[0],"testje")
#
#descrStats = descriptiveStatistics(d)
# bListOfMaps=convertListOfValuesToListOfNonSpatialMaps(descrStats,"cloneSmall.map")
# pcr.report(bListOfMaps[0],"mean")
#
#import time
# time.sleep(100)

def keepHistoryOfMaps(currentHistoryOfMaps, mapOfCurrentTimeStep, numberOfTimeStepsToKeep):
  # uses deque objects instead of lists, thus, conversion is required to get a list simply
  # by list(currentHistoryOfMaps)
  currentHistoryOfMaps.append(mapOfCurrentTimeStep)
  if len(currentHistoryOfMaps) > numberOfTimeStepsToKeep:
    currentHistoryOfMaps.popleft()
  if len(currentHistoryOfMaps) > numberOfTimeStepsToKeep:
    print('warning: length of keepHistoryOfMaps is greater than number of timesteps to keep')
  return currentHistoryOfMaps


##############################
#     map algebra            #
##############################

def nrCols(map):
  x = pcr.xcoordinate(pcr.boolean(map))
  xMax = pcr.mapmaximum(x)
  xMin = pcr.mapminimum(x)
  nrCols = ((xMax - xMin) / pcr.celllength()) + 1
  return nrCols

def nrRows(map):
  y = pcr.ycoordinate(pcr.boolean(map))
  yMax = pcr.mapmaximum(y)
  yMin = pcr.mapminimum(y)
  nrCols = ((yMax - yMin) / pcr.celllength()) + 1
  return nrCols

def nrCells(map):
  return nrCols(map) * nrRows(map)

def corners(map):
  left = edge(map, 4, 0)
  right = edge(map, 6, 0)
  top = edge(map, 8, 0)
  bottom = edge(map, 2, 0)
  return pcr.pcrgt(pcr.scalar(left) + pcr.scalar(right) + pcr.scalar(top) + pcr.scalar(bottom), 1.5)

def edges(map):
  left = edge(map, 4, 0)
  right = edge(map, 6, 0)
  top = edge(map, 8, 0)
  bottom = edge(map, 2, 0)
  return pcr.pcror(pcr.pcror(left, right), pcr.pcror(top, bottom))

def edgeZone(map, nrCells):
  # nrCells can be (should be) floating point
  edgeMap = edges(map)
  distToEdge = pcr.spread(edgeMap, 0, 1) / pcr.celllength()
  edgeZoneMap = distToEdge < (nrCells - 1.0)
  return edgeZoneMap


def booleanTrue(map):
  # returns a map that is everywhere true
  # removes mvs
  noMVs = pcr.cover(map, 1)
  return pcr.defined(noMVs)

def edge(map, side, distance):
  # map should have y incr. bot to top
  # returns boolean map with edge
  # distance defines distance to map boundary of edge,
  # e.g. distance 0 returns the real edge, distance
  # distance is an integer or nominal or ordinal
  # 1 returns one to the left/right/top/bottom
  # side is ldd dirs, e.g. 4 returns left side
  # side is an integer
  # works only with whole coordinates (as cells are selected using equals on coors as floating points)
  realDist = pcr.celllength() * pcr.scalar(distance)
  if ((side == 4) or (side == 6)):
    x = pcr.xcoordinate(booleanTrue(map))
    if side == 4:
      sideMap = pcr.pcreq(x, pcr.mapminimum(x) + realDist)
    if side == 6:
      sideMap = pcr.pcreq(x, pcr.mapmaximum(x) - realDist)
  if ((side == 8) or (side == 2)):
    y = pcr.ycoordinate(booleanTrue(map))
    if side == 2:
      sideMap = pcr.pcreq(y, pcr.mapminimum(y) + realDist)
    if side == 8:
      sideMap = pcr.pcreq(y, pcr.mapmaximum(y) - realDist)
  return sideMap

def bottom(map):
  '''
  returns the bottom line of cells
  works only with y increases bottom to top
  any cell size (unlike edge)
  '''
  yCoordinate = pcr.ycoordinate(pcr.defined(map))
  bottom = (yCoordinate == pcr.mapminimum(yCoordinate))
  return bottom

def neighbourIsMissingComponent(nbShift, noMVOnMap):
  NBIsMissingAll = pcr.ifthenelse(pcr.defined(nbShift), pcr.boolean(0), pcr.boolean(1))
  NBIsMissing = pcr.pcrand(NBIsMissingAll, noMVOnMap)
  return NBIsMissing

def neighbourIsMissingValueOrEdgeAndCellItselfIsDefined(map):
  noMVOnMap = pcr.defined(map)
  rightNBIsMissing = neighbourIsMissingComponent(pcr.shift(map, pcr.scalar(0), pcr.scalar(1)), noMVOnMap)
  leftNBIsMissing = neighbourIsMissingComponent(pcr.shift(map, pcr.scalar(0), pcr.scalar(-1)), noMVOnMap)
  upperNBIsMissing = neighbourIsMissingComponent(pcr.shift(map, pcr.scalar(-1), pcr.scalar(0)), noMVOnMap)
  lowerNBIsMissing = neighbourIsMissingComponent(pcr.shift(map, pcr.scalar(1), pcr.scalar(0)), noMVOnMap)
  return upperNBIsMissing, rightNBIsMissing, lowerNBIsMissing, leftNBIsMissing


# aMap=pcr.scalar("idOth.map")
#aMap=pcr.ifthen(pcr.uniform(1) < 0.9,pcr.scalar(2))
# one,two,three,four=neighbourIsMissingValueOrEdgeAndCellItselfIsDefined(aMap)
# pcr.report(one,'one.map')
# pcr.report(two,'two.map')
# pcr.report(three,'three.map')
# pcr.report(four,'four.map')
# pcr.report(aMap,'amap.map')


def moveRowsOrColumnsForPeriodicBoundaryCondition(map, direction):
  # moves row/column next to edge cells to other edge
  # direction is ldd dirs (e.g. 8 is bottom to top)
  # direction is an integer
  if (direction == 6):
    cellsToShift = edge(map, 4, 1)
    distanceToShift = pcr.scalar(nrCols(map) - 2)
    shiftedMap = pcr.ifthen(edge(map, 6, 0), pcr.shift(map, 0, 0 - distanceToShift))
  if (direction == 4):
    cellsToShift = edge(map, 6, 1)
    distanceToShift = pcr.scalar(nrCols(map) - 2)
    shiftedMap = pcr.ifthen(edge(map, 4, 0), pcr.shift(map, 0, distanceToShift))
  if (direction == 2):
    cellsToShift = edge(map, 8, 1)
    distanceToShift = pcr.scalar(nrRows(map) - 2)
    shiftedMap = pcr.ifthen(edge(map, 2, 0), pcr.shift(map, 0 - distanceToShift, 0))
  if (direction == 8):
    cellsToShift = edge(map, 2, 1)
    distanceToShift = pcr.scalar(nrRows(map) - 2)
    shiftedMap = pcr.ifthen(edge(map, 8, 0), pcr.shift(map, distanceToShift, 0))
  return shiftedMap

def periodicBoundaryCondition(map):
  # note positive shifts are to left or to top..
  # first value is vertical shift, second value is horizontal shift
  left = moveRowsOrColumnsForPeriodicBoundaryCondition(map, 4)
  right = moveRowsOrColumnsForPeriodicBoundaryCondition(map, 6)
  top = moveRowsOrColumnsForPeriodicBoundaryCondition(map, 8)
  bottom = moveRowsOrColumnsForPeriodicBoundaryCondition(map, 2)
  tmp = pcr.cover(left, right, top, bottom, map)
  newMap = pcr.ifthen(pcr.pcrnot(corners(map)), tmp)
  return newMap

# test=pcr.scalar("idOth.map")
# nrCols=nrCols(test)
# pcr.report(nrCols,"test")
# testje=edge(test,8,1)
# pcr.report(testje,"testje")
#
# test2=moveRowsOrColumnsForPeriodicBoundaryCondition(test,8)
# pcr.report(test2,"test2")
# test3=periodicBoundaryCondition(test)
# pcr.report(test3,"test3")
# test4=edges(test)
# pcr.report(test4,"test4")

def periodicBoundaryConditionNumpy(map):
  a = pcr.pcr2numpy(map, 1)
  # second left edge
  secondLeft = a[:, 1]
  # second right edge
  secondRight = a[:, a.shape[1] - 2]
  # remove left and right edge
  new = numpy.delete(a, (0, a.shape[1] - 1), 1)
  # add new left edge
  newLeft = numpy.insert(new, 0, secondRight, 1)
  # add new right edge
  b = numpy.insert(newLeft, newLeft.shape[1], secondLeft, 1)
  # second upper edge
  secondUpper = b[1, :]
  # second right edge
  secondLower = b[a.shape[0] - 2, :]
  # remove upper and lower edge
  bNew = numpy.delete(b, (0, a.shape[0] - 1), 0)
  # add upper edge
  newTop = numpy.insert(bNew, 0, secondLower, 0)
  # add new lower edge
  c = numpy.insert(newTop, newTop.shape[0], secondUpper, 0)
  outMap = pcr.numpy2pcr(Scalar, c, 0)
  return outMap

# test=pcr.scalar("idOth.map")
# test2=periodicBoundaryConditionNumpy(test)
# pcr.report(test2,"test2")
# test3=periodicBoundaryCondition(test)
# pcr.report(test3,"test3")

def createToCellsPeriodicBoundaryCondition(clone):
  list = []
  list.append(edge(clone, 6, 0))
  list.append(edge(clone, 8, 0))
  list.append(edge(clone, 4, 0))
  list.append(edge(clone, 2, 0))
  return list

def createFromCellsPeriodicBoundaryCondition(clone):
  list = []
  list.append(edge(clone, 4, 1))
  list.append(edge(clone, 2, 1))
  list.append(edge(clone, 6, 1))
  list.append(edge(clone, 8, 1))
  return list

def colNumber(clone):
  return pcr.ordinal((pcr.roundoff((pcr.xcoordinate(clone) / pcr.celllength()))))

def rowNumber(clone):
  return pcr.ordinal((pcr.roundoff((pcr.ycoordinate(clone) / pcr.celllength()))))

def periodicBoundaryConditionAreatotal(map, fromCells, toCells, cols, rows):
  right = pcr.ifthen(toCells[0], pcr.areatotal(pcr.ifthen(fromCells[0], map), rows))
  bottom = pcr.ifthen(toCells[1], pcr.areatotal(pcr.ifthen(fromCells[1], map), cols))
  left = pcr.ifthen(toCells[2], pcr.areatotal(pcr.ifthen(fromCells[2], map), rows))
  top = pcr.ifthen(toCells[3], pcr.areatotal(pcr.ifthen(fromCells[3], map), cols))
  return pcr.cover(right, bottom, left, top, map)


# test=pcr.scalar("idOth.map")
# clone=pcr.defined(test)
#
# fromCells=createFromCellsPeriodicBoundaryCondition(clone)
# toCells=createToCellsPeriodicBoundaryCondition(clone)
# cols=colNumber(clone)
# rows=rowNumber(clone)
# test2=periodicBoundaryConditionAreatotal(test,fromCells,toCells,cols,rows)
# pcr.report(test,"test")
# pcr.report(test2,"test2")
# test3=periodicBoundaryCondition(test)
# pcr.report(test2-test3,"diff")

##########################
#     sampling schemes   #
##########################

def samplingScheme(clone, nrSamples, fractionShortDistance, separationDistance, nrCellsToRight, nrCellsToTop):
  # get numb. of samples
  nrSamplesGrid = pcr.rounddown((1.0 - fractionShortDistance) * nrSamples)
  nrSamplesShortDistance = nrSamples - nrSamplesGrid
  # get possible locs
  colnumber = pcr.roundoff((pcr.xcoordinate(clone) / pcr.celllength()))
  tmp = pcr.roundoff((pcr.ycoordinate(clone) / pcr.celllength()))
  rownumber = tmp - pcr.mapminimum(tmp) + 1
  possibleCols = pcr.pcreq(pcr.pcrmod(colnumber, separationDistance), 0)
  possibleRows = pcr.pcreq(pcr.pcrmod(rownumber, separationDistance), 0)
  possibleLocations = pcr.pcrand(possibleCols, possibleRows)
  # get grid area
  nrRowsCols = pcr.roundup(pcr.sqrt(nrSamplesGrid))
  locColSel = pcr.ifthenelse(pcr.pcrle(colnumber, (separationDistance * nrRowsCols)), possibleLocations, 0)
  locRowSel = pcr.ifthenelse(pcr.pcrle(rownumber, (separationDistance * nrRowsCols)), possibleLocations, 0)
  locsArea = pcr.pcrand(locColSel, locRowSel)
  # get samples
  samplesAtGrid = pcr.cover(pcr.pcrle(pcr.order(pcr.ifthen(locsArea, pcr.uniqueid(clone))), nrSamplesGrid), 0)
  # get samples with extra samples
  randomSampleNumbers = pcr.ordinal(pcr.ifthen(samplesAtGrid, pcr.order(pcr.uniform(pcr.ifthen(samplesAtGrid, pcr.boolean(1))))))
  alreadyCovered = pcr.defined(randomSampleNumbers)
  for sample in range(1, int(getCellValue(pcr.mapmaximum(randomSampleNumbers), 1, 1)) + 1):
    theSample = pcr.cover(pcr.pcreq(sample, randomSampleNumbers), 0)
    nbSample = pcr.pcrand(pcr.pcrgt(pcr.window4total(pcr.scalar(theSample)), 0.5), pcr.pcrnot( alreadyCovered))
    extraSample = pcr.cover(pcr.pcrlt(pcr.order(pcr.ifthen(nbSample, pcr.uniform(1))), 1.5), 0)
    alreadyCovered = pcr.pcror(alreadyCovered, extraSample)
    totalSamples = pcr.maptotal(pcr.scalar(alreadyCovered))
    # print getCellValue(totalSamples,1,1)
    if getCellValue(totalSamples, 1, 1) > (nrSamples - 1.0 + 0.1):
      break
  sampleNumbers = pcr.ordinal(pcr.ifthen(alreadyCovered, pcr.order(pcr.uniqueid(pcr.ifthen(alreadyCovered, pcr.boolean(1))))))
  # shift the cells to centre
  centreCol = pcr.maptotal(colnumber) / nrCells(clone)
  centreRow = pcr.maptotal(rownumber) / nrCells(clone)
  centreColSamples = pcr.maptotal(pcr.ifthen(alreadyCovered, colnumber)) / pcr.maptotal(pcr.ifthen(alreadyCovered, pcr.scalar(1)))
  centreRowSamples = pcr.maptotal(pcr.ifthen(alreadyCovered, rownumber)) / pcr.maptotal(pcr.ifthen(alreadyCovered, pcr.scalar(1)))
  sampleNumbersShifted = pcr.shift(sampleNumbers, centreRow - centreRowSamples - nrCellsToTop, centreColSamples - centreCol - nrCellsToRight)
  return pcr.cover(sampleNumbersShifted, 0)

def samplingSchemeRandomShift(clone, nrSamples, fractionShortDistance, separationDistance, maxNrCellsToRight, maxNrCellsToTop):
  uniformRight = pcr.mapuniform()
  shiftRight = pcr.roundoff((uniformRight - 0.5) * maxNrCellsToRight)
  uniformTop = pcr.mapuniform()
  shiftTop = pcr.roundoff((uniformTop - 0.5) * maxNrCellsToTop)
  sampleNumbersShifted = samplingScheme(clone, nrSamples, fractionShortDistance, separationDistance, shiftRight, shiftTop)
  return sampleNumbersShifted

def samplingSchemeSubset(clone, nrSamples, fractionShortDistance, separationDistance, nrCellsToRight, nrCellsToTop, nrSamplesRemove, uniformMap,
                         realNrSamples ):
  # nrSamples, nr of samples as if it were a normal scheme (regular, all grid positions used)
  # nrSamplesRemove, nr of samples removed again from normal scheme, random positions
  # realNrSamples real number of samples (to represent removal)

  # get numb. of samples
  nrSamplesGrid = pcr.rounddown((1.0 - fractionShortDistance) * nrSamples)
  nrSamplesShortDistance = nrSamples - nrSamplesGrid
  # get possible locs
  colnumber = pcr.roundoff((pcr.xcoordinate(clone) / pcr.celllength()))
  tmp = pcr.roundoff((pcr.ycoordinate(clone) / pcr.celllength()))
  rownumber = tmp - pcr.mapminimum(tmp) + 1
  possibleCols = pcr.pcreq(pcr.pcrmod(colnumber, separationDistance), 0)
  possibleRows = pcr.pcreq(pcr.pcrmod(rownumber, separationDistance), 0)
  possibleLocations = pcr.pcrand(possibleCols, possibleRows)
  # get grid area
  nrRowsCols = pcr.roundup(pcr.sqrt(nrSamplesGrid))
  locColSel = pcr.ifthenelse(pcr.pcrle(colnumber, (separationDistance * nrRowsCols)), possibleLocations, 0)
  locRowSel = pcr.ifthenelse(pcr.pcrle(rownumber, (separationDistance * nrRowsCols)), possibleLocations, 0)
  locsArea = pcr.pcrand(locColSel, locRowSel)
  # get samples
  samplesAtGrid = pcr.cover(pcr.pcrle(pcr.order(pcr.ifthen(locsArea, pcr.uniqueid(clone))), nrSamplesGrid), 0)
  # randomly remove samples
  samplesAtGridRandomSampleNumbers = pcr.order(pcr.ifthen(pcr.pcrne(samplesAtGrid, 0), uniformMap) )
  pcr.report(samplesAtGridRandomSampleNumbers, 'rem.map')
  remove = pcr.cover(samplesAtGridRandomSampleNumbers < nrSamplesRemove, 0)
  newSamples = pcr.cover(pcr.ifthenelse(pcr.pcrnot(remove), samplesAtGrid, 0), 0)
  samplesAtGrid = newSamples
  pcr.report(samplesAtGrid, 'testje.map')
  # get samples with extra samples
  randomSampleNumbers = pcr.ordinal(pcr.ifthen(samplesAtGrid, pcr.order(pcr.uniform(pcr.ifthen(samplesAtGrid, pcr.boolean(1))))))
  alreadyCovered = pcr.defined(randomSampleNumbers)
  for sample in range(1, int(getCellValue(pcr.mapmaximum(randomSampleNumbers), 1, 1)) + 1):
    theSample = pcr.cover(pcr.pcreq(sample, randomSampleNumbers), 0)
    nbSample = pcr.pcrand(pcr.pcrgt(pcr.window4total(pcr.scalar(theSample)), 0.5), pcr.pcrnot( alreadyCovered))
    extraSample = pcr.cover(pcr.pcrlt(pcr.order(pcr.ifthen(nbSample, pcr.uniform(1))), 1.5), 0)
    alreadyCovered = pcr.pcror(alreadyCovered, extraSample)
    totalSamples = pcr.maptotal(pcr.scalar(alreadyCovered))
    # print getCellValue(totalSamples,1,1)
    if getCellValue(totalSamples, 1, 1) > (realNrSamples - 1.0 + 0.1):
      break
  sampleNumbers = pcr.ordinal(pcr.ifthen(alreadyCovered, pcr.order(pcr.uniqueid(pcr.ifthen(alreadyCovered, pcr.boolean(1))))))
  # shift the cells to centre
  centreCol = pcr.maptotal(colnumber) / nrCells(clone)
  centreRow = pcr.maptotal(rownumber) / nrCells(clone)
  centreColSamples = pcr.maptotal(pcr.ifthen(alreadyCovered, colnumber)) / pcr.maptotal(pcr.ifthen(alreadyCovered, pcr.scalar(1)))
  centreRowSamples = pcr.maptotal(pcr.ifthen(alreadyCovered, rownumber)) / pcr.maptotal(pcr.ifthen(alreadyCovered, pcr.scalar(1)))
  sampleNumbersShifted = pcr.shift(sampleNumbers, centreRow - centreRowSamples - nrCellsToTop, centreColSamples - centreCol - nrCellsToRight)
  return pcr.cover(sampleNumbersShifted, 0)

def samplingSchemeRandomShiftSubset(clone, nrSamples, fractionShortDistance, separationDistance, maxNrCellsToRight, maxNrCellsToTop, nrSamplesRemove,
                                    uniformMap, realNrSamples):
  uniformRight = pcr.mapuniform()
  shiftRight = pcr.roundoff((uniformRight - 0.5) * maxNrCellsToRight)
  uniformTop = pcr.mapuniform()
  shiftTop = pcr.roundoff((uniformTop - 0.5) * maxNrCellsToTop)
  sampleNumbersShifted = samplingSchemeSubset(clone, nrSamples, fractionShortDistance, separationDistance, shiftRight, shiftTop, nrSamplesRemove, uniformMap,
                       realNrSamples)
  return sampleNumbersShifted



# clone=pcr.defined("clone.map")
# test=samplingSchemeRandomShift(clone,100,0.3,3,4,40)
# test2=samplingScheme(clone,100,0.3,3,0,0)
# pcr.report(test,"pietje")
# pcr.report(test2,"pietje2")

def createGifAnimation(name, samples):
  # postmcloop function to convert individual png or gif
  # files to animated gif
  # use gifview (from gifsicle package) to animate, man gifview for
  # options
  # or gimp, filters -> animation -> playback
  for sample in samples:
    baseName = pcrfw.generateNameS(name, sample)
    command = "convert " + baseName + "* " + baseName + "_ani.gif"
    os.system(command)
    # alternatief: for i in ete*; do convert $i ${i%%png}gif; done
    # en dan gifsicle to make the animation (maybe faster)

def mapaverage(aMap):
  total = pcr.maptotal(aMap)
  nrCellsNoMV = pcr.maptotal(pcr.scalar(pcr.defined(aMap)))
  return total / nrCellsNoMV

########################################
### FUNCTIONS FOR PARTICLE FILTERING ###
########################################

def createFileNameToReportAVariableForSuspend(classInstance, variableName, currentTimeStep, currentSampleNumber):
  className = classInstance.__class__.__name__
  classNameDir = str(currentSampleNumber) + '/stateVar/' + className
  if not os.path.exists(classNameDir):
    os.mkdir(classNameDir)
  return pcrfw.generateNameST('stateVar/' + className + '/' + variableName, currentSampleNumber, currentTimeStep)
  # print className + "/" variableName

def letterSequence():
  alphabeth = string.ascii_letters[0:26]
  letters = []
  for firstLetter in alphabeth:
    for secondLetter in alphabeth:
      letters.append(firstLetter + secondLetter)
  return letters

def reportMemberVariablesOfAClassForSuspend(classInstance, currentTimeStep, currentSampleNumber):
  b = vars(classInstance)
  names = letterSequence()
  i = 0
  for name in sorted(b):
    if type(b[name]) == PCRaster._PCRaster.Field:
      fileName = createFileNameToReportAVariableForSuspend(classInstance, names[i], currentTimeStep, currentSampleNumber)
      pcr.report(b[name], fileName)
    # switch on for testing
    #fileName = createFileNameToReportAVariableForSuspend(classInstance,names[i],currentTimeStep,currentSampleNumber)
    # print 'report', classInstance, name, fileName
    i = i + 1
  print('In reporting, the number of member vars in ', classInstance, ' is ', i)

def readMemberVariablesOfAClassForResume(classInstance, currentTimeStep, currentSampleNumber):
  b = vars(classInstance)
  names = letterSequence()
  i = 0
  for name in sorted(b):
    if type(b[name]) == PCRaster._PCRaster.Field:
      fileName = createFileNameToReportAVariableForSuspend(classInstance, names[i], currentTimeStep, currentSampleNumber)
      mapToResume = pcr.readmap(fileName)
      vars(classInstance)[name]  = mapToResume
    # switch on for testing
    #fileName = createFileNameToReportAVariableForSuspend(classInstance,names[i],currentTimeStep,currentSampleNumber)
    # print 'read', classInstance, name, fileName
    i = i + 1
  print('In reading, the number of member vars in ', classInstance, ' is ', i)

def printMemberVariables(classInstance):
  a = sorted(vars(classInstance))
  print(classInstance)
  print('number of member variables ', len(a))
  print(a)
  print

def removePeriodsFromAListOfTimesteps(filterTimesteps, periodsToExclude):
  newFilterTimesteps = filterTimesteps[:]
  for period in periodsToExclude:
    lower = period[0]
    upper = period[1]
    oldFilterTimesteps = newFilterTimesteps[:]
    newFilterTimesteps = []
    for timestep in oldFilterTimesteps:
      if (timestep < lower) or (timestep > upper):
        newFilterTimesteps.append(timestep)
  return newFilterTimesteps


##########################
### OBJECTIVE FUNCTIONS ##
##########################


def hydrographOF(observedDischarge, modelledDischarge):
  """Returns nash sutcliffe coefficient and mean square error

  observedDischarge -- numpy array
  modelledDischarge -- numpy array
  """
  squares = numpy.square(observedDischarge - modelledDischarge)
  squaresO = numpy.square(observedDischarge - numpy.mean(observedDischarge))
  NSE = 1.0 - (numpy.sum(squares) / numpy.sum(squaresO))
  MSE = numpy.sum(squares) / len(squares)
  a = numpy.corrcoef(observedDischarge, modelledDischarge)
  sumObs = numpy.sum(observedDischarge)
  sumMod = numpy.sum(modelledDischarge)
  slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(observedDischarge, modelledDischarge)
  print('sumObs is', sumObs)
  print('sumMod is', sumMod)
  print('lin. regression; slope, intercept, r2, p value are', slope, intercept, r_value**2.0, p_value)
  return NSE, MSE

def hydrographOFProb(observedDischarge, modelledDischarge):
  """Returns nash sutcliffe coefficient and mean square error, probabilistic, i.e. modelled is a MC sample

  observedDischarge -- numpy array
  modelledDischarge -- numpy array
  """
  squares = numpy.square(observedDischarge - modelledDischarge.T)
  squaresO = numpy.square(observedDischarge - numpy.mean(observedDischarge))
  NSE = 1.0 - (numpy.sum(squares, axis=1) / numpy.sum(squaresO))
  # MSE=numpy.sum(squares)/len(squares)
  # sumObs=numpy.sum(observedDischarge)
  # sumMod=numpy.sum(modelledDischarge)
  # print 'sumObs is', sumObs
  # print 'sumMod is', sumMod
  return NSE



##############################################
### FUNCTIONS FOR GENERATING RANDOM FIELDS ###
##############################################

def mapuniformBounds(min, max, fixedValue='aString..', useRealization=True):
  """Assigns a random value taken from a uniform distribution

  min -- lower bound, python floating point or pcraster type (no default)
  max -- upper bound of interval, python or pcraster type (no default)
  fixedValue -- value assigned when useRealization is False (python floating point or pcraster type)
  useRealization -- switch to assign fixedValue (default is True)

  fixedValue should be a PCRaster type when output is used in PCRaster functions
  """
  if useRealization:
    a = pcr.mapuniform()
    range = pcr.scalar(max) - pcr.scalar(min)
    field = (a * range) + min
    return field
  else:
    return fixedValue

def areauniformBounds(min, max, areaMap, fixedValue='aString..', useRealization=True):
  """Assigns a random value taken from a uniform distribution

  min -- lower bound, python or pcraster type
  max -- upper bound of interval, python or pcraster type
  areamap -- classified map
  fixedValue -- value assigned when useRealization is False
  useRealization -- switch to assign fixedValue

  fixedValue should be a PCRaster type when output is used in PCRaster functions"""
  if useRealization:
    a = pcr.areauniform(pcr.spatial(areaMap))
    range = pcr.scalar(max) - pcr.scalar(min)
    field = (a * range) + min
    return field
  else:
    return fixedValue

def mapNormalRelativeError(input, standardDeviationForInputOfOne):
  a = pcr.mapnormal()
  error = a * input * standardDeviationForInputOfOne
  realization = input + error
  return realization

def mapgamma(shapeParameter):
  '''Returns a realization from the gamma distribution with a mean of one
  shapeParameter is a Python floating point
  return value is a Python floating point
  '''

  scaleParameter = 1.0 / shapeParameter
  realization = random.gammavariate(shapeParameter, scaleParameter)
  return realization

def booleanMapWithOnlyMissingValues(map):
  nonMV = pcr.defined(map)
  result = pcr.ifthen(pcr.pcrnot(nonMV), pcr.boolean(map))
  return result

def lookupColumnsInTableScalar(table, classifiedMap):
  '''Reads 2nd to nth column in table and returns for each of these
  columns a map. Maps are returned as a list of maps. First column is
  key column linked to the classified input map classifiedMap. Returns
  a scalar maps.
  table         -- input table giving for each class on classifiedMap a value
  classifiedMap -- classes
  Example of a table (ascii file):
  1  0.2 0.5  0.3
  2  0.9 0.05 0.05
  3  0.0 0.01 0.99
  First column should give all codes on classifiedMap, second up to n columns
  give values assigned to output, second column value of code 1, third
  value of code 2, etc.
  '''
  # read the table
  tableNP = numpy.loadtxt(table)

  # create the output, empty
  emptyMap = pcr.scalar(booleanMapWithOnlyMissingValues(classifiedMap))
  nrOfOutputMaps = len(tableNP[0]) - 1
  resultMaps = []
  i = 0
  while i < nrOfOutputMaps:
    resultMaps.append(emptyMap)
    i += 1

  for row in tableNP:
    key = row[0]
    outputs = row[1:]
    areaToBeAssigned = pcr.pcreq(pcr.nominal(key), classifiedMap)
    outputNumber = 0
    for value in outputs:
      resultMap = pcr.ifthenelse(areaToBeAssigned, value, resultMaps[outputNumber])
      resultMaps[outputNumber] = resultMap
      outputNumber += 1
  return resultMaps

def discreteProbabilityDistributionPerArea(table, classifiedMap):
  '''Draws a realization from a discrete probability distribution given
  for each class. All cells in a class get the same realization.
  table         -- input table giving for each class on classifiedMap the distribution
  classifiedMap -- classes
  Returns a nominal map with values between 1 and the number of columns in table - 1.
  Example of a table (ascii file):
  1  0.2 0.5  0.3
  2  0.9 0.05 0.05
  3  0.0 0.01 0.99
  First column should give all codes on classifiedMap, second up to n columns
  give probability distribution, second column probability of code 1, third
  probability of code 2, etc. Result may contain values between 1 and n-1.
  '''
  randomValue = pcr.areauniform(classifiedMap)
  result = discreteProbabilityDistribution(table, classifiedMap, randomValue)
  return result

def discreteProbabilityDistributionPerCell(table, classifiedMap):
  '''Same as discreteProbabilityDistributionPerArea, but each cell in output
  is a independent realization
  '''
  randomValue = pcr.uniform(pcr.boolean(1))
  result = discreteProbabilityDistribution(table, classifiedMap, randomValue)
  return result

def discreteProbabilityDistribution(table, classifiedMap, uniformMap):
  '''Draws a realization from a discrete probability distribution given
  for each class. All cells in a class get the same realization
  table         -- input table giving for each class on classifiedMap the distribution
  classifiedMap -- classes
  uniformMap    -- random values between zero and one used to create the realizations
  Returns a nominal map with values between 1 and the number of columns in table - 1.
  Example of a table (ascii file):
  1  0.2 0.5  0.3
  2  0.9 0.05 0.05
  3  0.0 0.01 0.99
  First column should give all codes on classifiedMap, second up to n columns
  give probability distribution, second column probability of code 1, third
  probability of code 2, etc. Result may contain values between 1 and n-1.
  '''
  # create series of maps with the probabilities
  probabilities = lookupColumnsInTableScalar(table, classifiedMap)

  # get number of output classes
  tableNP = numpy.loadtxt(table)
  nrOfOutputClasses = len(tableNP[0]) - 1

  # calculate cumulative probabilities
  cumulativeProbabilities = []
  cumulativeProbability = pcr.scalar(0)
  for probability in probabilities:
    cumulativeProbability = probability + cumulativeProbability
    cumulativeProbabilities.append(cumulativeProbability)

  randomValue = uniformMap

  # create empty result map
  result = pcr.nominal(booleanMapWithOnlyMissingValues(classifiedMap))
  value = nrOfOutputClasses
  reversedCumProbs = reversed(cumulativeProbabilities)
  for cumProb in reversedCumProbs:
    print(value)
    result = pcr.ifthenelse(pcr.pcrlt(randomValue, cumProb), value, result)
    value = value - 1

  return result



# a=0.0
# for i in range(1,100000):
#  test = mapgamma(6.0)
#  a=a+test
# print 'gamma real :', test
# print 'gamma mean :', a/100000.0



###############################
### FUNCTIONS FOR REPORTING ###
###############################

def reportAListOfMaps(listOfMaps, baseName, timeStep, sample):
  alphabeth = string.ascii_letters[0:26]
  i = 0
  for map in listOfMaps:
    totBaseName = baseName + str(alphabeth[i])
    reportName = pcrfw.generateNameST(totBaseName, sample, timeStep)
    pcr.report(map, reportName)
    i = i + 1

def reportAVariableAtTheLastTimeStep(currentTimeStep, currentSampleNumber, nrTimeSteps, variable, name):
  if currentTimeStep == nrTimeSteps:
    pcr.report(variable, pcrfw.generateNameS(name + '.map', currentSampleNumber))

#####################################
### FUNCTIONS FOR FRAGSTATS STUFF ###
#####################################

def proportionOfClassInAreas(areas, classMap, selectedClass):
  '''
  Calculates the fractional area of a class on classMap
  for each area on areas map. The class on classMap
  is given by selectedClass.
  '''
  areaOfAreas = pcr.areatotal(pcr.spatial(pcr.scalar(1)), areas)
  selectedClassMap = pcr.scalar((classMap == selectedClass))
  proportionOfClassInAreasMap = pcr.areatotal(selectedClassMap, areas) / areaOfAreas
  return proportionOfClassInAreasMap

def convertValuesInAreasToSeparateMaps(areas, numberOfAreas, valuesInAreas):
  '''
  Takes the maximum value of valuesInAreas in each area on areas and
  stores it in a map. These maps are collected and returned as a list
  of maps.
  '''
  separateMaps = []
  for i in range(1, numberOfAreas + 1):
    separateMap = pcr.mapmaximum(pcr.ifthen(areas == i, valuesInAreas))
    separateMaps.append(separateMap)
  return separateMaps

def selectFromEachAreaOneCell(areas):
  '''
  Selects from each class on areas one cell and returns
  the class value from area at that cell. Remaining cells
  have zero value at result.
  '''
  uniqueId = pcr.uniqueid(pcr.defined(areas))
  oneCellPerArea = pcr.ifthenelse(pcr.areamaximum(uniqueId, areas) == uniqueId, areas, 0)
  return oneCellPerArea

def patchSize(areas):
  areasClumped = pcr.clump(areas)
  nrCellsPerPatch = pcr.areatotal(pcr.spatial(pcr.scalar(1)), areasClumped)
  singleCellPerArea = selectFromEachAreaOneCell(areasClumped) != 0
  nrCellsPerPatchStoredOnOneCellPerPatch = pcr.ifthen(singleCellPerArea, nrCellsPerPatch)
  meanPatchSize = pcr.areaaverage(nrCellsPerPatchStoredOnOneCellPerPatch, areas)
  return meanPatchSize


# areas=(pcr.nominal('cl000000.014'))
# pcr.report(areas,'areas.map')
# test=selectFromEachAreaOneCell(areas)
# pcr.report(test,'test.map')
# patchsize=patchSize(areas)
# pcr.report(patchsize,'patchsize.map')

