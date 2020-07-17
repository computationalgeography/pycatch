import string
import os, shutil, numpy
import generalfunctions
import matplotlib.pyplot as plt
from PCRaster import *
from PCRaster.Framework import *

#def varmean(names,sampleNumbers, timeSteps):
#  nrSamples=scalar(len(sampleNumbers))
#  for name in names:
#    for step in timeSteps:
#      sumSquared=scalar(0.0)
#      sum=scalar(0.0)
#      for sample in sampleNumbers:
#        realization=scalar(generateNameST(name,sample,step))
#        sumSquared=sumSquared+realization*realization
#        sum=sum+realization
#      std=(scalar(1.0)/nrSamples)*sqrt(nrSamples*sumSquared-sum*sum)
#      var=std*std
#      mean=sum/nrSamples
#      report(var, generateNameT(name + '-var', step))
#      report(mean, generateNameT(name + '-ave', step))


def correlation(location, independentName, dependentName, locationName, sampleNumbers, timeSteps):
  location=boolean(location)
  name=independentName + '_' + dependentName + '_' + locationName
  tssFileIntercept = file("%s%s.tss" % (name,'_int'), "w")
  tssFileSlope = file("%s%s.tss" % (name,'_slope'), "w")
  tssFileRSquared = file("%s%s.tss" % (name,'_rSq'), "w")
  for step in timeSteps:
    print step, '#############'
    values=[]
    for sample in sampleNumbers:
      smallValue=0.0000000000000000001
      fileNameOne=generateNameST(dependentName,sample,step)
      valueOne=generalfunctions.getCellValueAtBooleanLocation(location,scalar(fileNameOne))
      pairList=[valueOne + smallValue]
      fileNameTwo=generateNameST(independentName,sample,step)
      valueTwo=generalfunctions.getCellValueAtBooleanLocation(location,scalar(fileNameTwo))
      pairList.append(valueTwo + smallValue)
      values.append(pairList)
      print valueOne + smallValue, valueTwo + smallValue
    reg=regression.linearRegression(values, 1)
    rSq=regression.linearRSquared(values,reg)
    print step, reg, rSq
    tssFileIntercept.write("%d %g\n" % (step, reg[0]))
    tssFileSlope.write("%d %g\n" % (step, reg[1]))
    tssFileRSquared.write("%d %g\n" % (step, rSq))
  tssFileIntercept.close()
  tssFileSlope.close()
  tssFileRSquared.close()

def staticInput(timeSteps):
  return len(timeSteps) == 1 and timeSteps[0] == 0

def deterministicInput(sampleNumbers):
  return len(sampleNumbers) == 1 and sampleNumbers[0] == 0

def minOfSamples(name, sampleNumbers):
  """Calculates the minimum value of each cell.

  name             Name of the scalar raster for which each sample has a
                   realization.
  sampleNumbers    List of numbers of samples to aggregate.
  Returns a raster with minimum values.
  """
  minimum = scalar(1e31)
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = scalar(readmap(filename))
    minimum      = ifthenelse(pcrlt(raster,minimum),raster,minimum)
  return minimum

def maxOfSamples(name, sampleNumbers):
  """Calculates the maximum value of each cell.

  name             Name of the scalar raster for which each sample has a
                   realization.
  sampleNumbers    List of numbers of samples to aggregate.
  Returns a raster with maximum values.
  """
  maximum = scalar(-1e31)
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = scalar(readmap(filename))
    maximum      = ifthenelse(pcrgt(raster,maximum),raster,maximum)
  return maximum

def uniquesamples(name, sampleNumbers):
  """Retrieves the unique samples.

  name             Name of the raster for which each sample has a
                   realization.
  sampleNumbers    List of numbers of samples to aggregate.
  Returns a list with sets of corresponding loops.
  """
  uniqueSets=[]
  for sample in sampleNumbers:
    filename = generateNameS(name, sample)
    raster   = readmap(filename)
    setNumber=0
    sampleAddedToExistingSet=False
    for uniqueSet in uniqueSets:
      if generalfunctions.mapeq(uniqueSet[0],raster):
        uniqueSets[setNumber].append(sample)
        sampleAddedToExistingSet=True
        break
      setNumber += 1
    if sampleAddedToExistingSet==False:
      uniqueSets.append([raster,sample])
  firstLoopOfEachUniqueSet=[]
  for uniqueSet in uniqueSets:
    firstLoopOfEachUniqueSet.append(uniqueSet[1])
  return len(firstLoopOfEachUniqueSet),firstLoopOfEachUniqueSet

def mcaveragevariance(names,sampleNumbers, timeSteps):
  if staticInput(timeSteps):
    for name in names:
      mean=average(name + '.map', sampleNumbers)
      var=variance(name + '.map', sampleNumbers)
      minimum=minOfSamples(name + '.map', sampleNumbers)
      maximum=maxOfSamples(name + '.map', sampleNumbers)
      #std=stddev(name + '.map', sampleNumbers)
      report(mean, name + '-ave.map')
      report(var, name + '-var.map')
      report(minimum, name + '-min.map')
      report(maximum, name + '-max.map')
      report(sqrt(var)/mean, name + '-err.map')
  else:
    nrSamples=scalar(len(sampleNumbers))
    for name in names:
      for step in timeSteps:
        var=variance(generateNameT(name,step),sampleNumbers)
        mean=average(generateNameT(name,step), sampleNumbers)
        report(mean, generateNameT(name + '-ave', step))
        report(var, generateNameT(name + '-var', step))
        report(sqrt(var)/mean, generateNameT(name + '-err', step))

def mcpercentiles(names,percentiles,sampleNumbers,timeSteps):
  if staticInput(timeSteps):
    for name in names:
      results = percentile(name + '.map', sampleNumbers, percentiles)
      for i in range(len(percentiles)):
        report(results[i], (name + "_%s.map" % (percentiles[i])))
  else:
    for name in names:
      for step in timeSteps:
        results = percentile(generateNameT(name, step), sampleNumbers, percentiles)
        assert len(results) == len(percentiles)
        for i in range(len(percentiles)):
          report(results[i], (name + "_%d_%s.map" % (step, percentiles[i])))
          ## Add a record for this percentile to the timeseries file for a cell.
          #timeseriesValue=mapmaximum(ifthen('outflow.map',results[i]))
          ##value, valid = cellvalue(results[i], 4, 1); assert valid
          #value, valid = cellvalue(timeseriesValue, 1, 1); assert valid
          #timeseriesFiles[name][i].write("%d %g\n" % (step, value))
      #timeseriesFiles[name][i].close()  # dj this is required!

def createtimeseries(names, nameExtension, locations,sampleNumbers,timeSteps):
  if deterministicInput(sampleNumbers):
    for name in names:
      tssFile = file(name + nameExtension + '.tss', "w")
      tssFile.write("timeseries scalar\n")
      tssFile.write("2\n")
      tssFile.write("timestep\n")
      tssFile.write("%s\n" % (name))
      for step in timeSteps:
        timeseriesValue=mapmaximum(ifthen(locations,generateNameT(name,step)))
        value, valid = cellvalue(timeseriesValue, 1, 1); assert valid
        tssFile.write("%d %g\n" % (step, value))
      tssFile.close()
  else:
    for name in names:
      for sample in sampleNumbers:
        tssFile = file(generateNameS("%s%s.tss" % (name,nameExtension), sample), "w")
        tssFile.write("timeseries scalar\n")
        tssFile.write("2\n")
        tssFile.write("timestep\n")
        tssFile.write("%s\n" % (name))
        for step in timeSteps:
          timeseriesValue=mapmaximum(ifthen(locations,generateNameST(name,sample,step)))
          value, valid = cellvalue(timeseriesValue, 1, 1); assert valid
          tssFile.write("%d %g\n" % (step, value))
        tssFile.close()

def createtimeseriesnewfileformat(names,locations,sampleNumbers,timeSteps, quantiles):
  if deterministicInput(sampleNumbers):
    for quantile in quantiles:
      for name in names:
        tssFile = file(name + "_" + str(quantile) + '.tss', "w")
        tssFile.write("timeseries scalar\n")
        tssFile.write("2\n")
        tssFile.write("timestep\n")
        tssFile.write("%s\n" % (name))
        for step in timeSteps:
          filename = name + ("_%d_" % (step)) + str(quantile) + ".map"
          timeseriesValue=mapmaximum(ifthen(locations,filename))
          value, valid = cellvalue(timeseriesValue, 1, 1); assert valid
          tssFile.write("%d %g\n" % (step, value))
        tssFile.close()
  else:
    print 'timeseries for monte carlo loops not yet available'

def createGstatRealizations(setOfRealizations, nameCommandFile, nameOutMapList):
  # number of realizations required
  nSim=len(setOfRealizations)
  # open template gstat script and replace
  print nameCommandFile
  gstatTemplate=file(nameCommandFile + '.gst','r')
  gstatTemplateString=gstatTemplate.read()
  gstatTemplate.close()
  gstatString=string.replace(gstatTemplateString,'NSIM',str(nSim))
  gstatFile=file('tmpGstat.gst','w')
  gstatFile.write(gstatString)
  gstatFile.close()
  # run gstat
  os.system('gstat tmpGstat.gst')
  os.remove('tmpGstat.gst')
  # rename files
  i=1
  for realization in setOfRealizations:
    print realization
    item=0
    for name in nameOutMapList:
      gstatOutputFileName=generateNameT('g_' + name,i)
      #print gstatOutputFileName, realization[item]
      shutil.move(gstatOutputFileName,realization[item])
      item=item+1
    i = i + 1

def createAllGstatRealizations(nameCommandFile,nameOutMapList,nrRealPerGstatCall,sampleNumbers,timeSteps):
  # create list names with filenames required
  names=[]
  names.append([])
  i = 0
  j = 0
  for sample in sampleNumbers:
    for step in timeSteps:
      if j == nrRealPerGstatCall:
        j = 1
        i = i+1
        names.append([])
      else:
        j = j+1
      namesOneSampleOneTimeStep=[]
      for nameOutMap in nameOutMapList:
        if staticInput(timeSteps):
          fileName=generateNameS(nameOutMap,sample) + '.map'
        else:
          fileName=generateNameST(nameOutMap,sample,step)
        namesOneSampleOneTimeStep.append(fileName)
      names[i].append(namesOneSampleOneTimeStep)
  for setOfRealizations in names:
    createGstatRealizations(setOfRealizations,nameCommandFile, nameOutMapList)
    print setOfRealizations

def mySelectSArray(name, sampleNumbers, row, col):
  """Selects values at row, col from raster name in Monte Carlo samples.

  name -- Name of raster.
  sampleNumber -- Numbers of MC samples to use.
  row -- Row index of cell to read.
  col -- Col index of cell to read.
  The returned array does not contain missing values so the size is maximimal
  sampleNumbers but possibly smaller.

  Returned array has elements of type numpy.float32"""
  mask = numpy.zeros(len(sampleNumbers)).astype(numpy.bool_)
  array = numpy.zeros(len(sampleNumbers)).astype(numpy.float32)
  i = 0
  while i < len(sampleNumbers):
    filename = generateNameS(name, sampleNumbers[i])
    array[i], mask[i] = readFieldCell(filename, row, col)
    i += 1
  array = numpy.compress(mask, array)
  return array

def selectSArrayMultipleRasters(names,sampleNumbers,row,col):
  """Selects at row, col from each raster name
  Returned array is 'nested', i.e. each element contains
  an array with the values of a raster name"""
  a = []
  for name in names:
    arrayOfRaster = mySelectSArray(name,sampleNumbers,row,col)
    a.append(arrayOfRaster)
  c = numpy.vstack(a)
  return c

def covarMatrix(names,sampleNumbers,row,col,covarMatrixName, corrMatrixName):
  dataMatrix = selectSArrayMultipleRasters(names,sampleNumbers,row,col)
  covarMatrix = numpy.cov(dataMatrix)
  corrMatrix = numpy.corrcoef(dataMatrix)
  if len(names) == 1:
    covarMatrixDifferentType=numpy.array([float(covarMatrix)])
    corrMatrixDifferentType=numpy.array([float(corrMatrix)])
    numpy.savetxt(covarMatrixName,covarMatrixDifferentType)
    numpy.savetxt(corrMatrixName,corrMatrixDifferentType,fmt="%3.2f")
  else:
    numpy.savetxt(covarMatrixName,covarMatrix)
    numpy.savetxt(corrMatrixName,corrMatrix,fmt="%3.2f")

def mcCovarMatrix(names,sampleNumbers,timeSteps,row,col,covarMatrixBaseName,corrMatrixBaseName):
  for step in timeSteps:
    namesForTimestep=[]
    for name in names:
      namesForTimestep.append(generateNameT(name,step)) 
    covarMatrixName=generateNameT(covarMatrixBaseName,step)
    corrMatrixName=generateNameT(corrMatrixBaseName,step)
    covarMatrix(namesForTimestep, sampleNumbers, row, col, covarMatrixName,corrMatrixName)

def mySelectTArray(name, timeSteps, row, col):
  # kan weg zodra alles is overgezet op
  # def mySelectTArrayWithTimesteps(name, timeSteps, row, col):
  """Selects values at row, col from raster name in time steps.

  name -- Name of raster.
  timeSteps -- Numbers of time steps to use.
  row -- Row index of cell to read, starts with 1!.
  col -- Col index of cell to read, starts with 1!.

  Returned array has elements of type numpy.float32"""
  mask = numpy.zeros(len(timeSteps)).astype(numpy.bool_)
  array = numpy.zeros(len(timeSteps)).astype(numpy.float32)
  i = 0
  while i < len(timeSteps):
    filename = generateNameT(name, timeSteps[i])
    array[i], mask[i] = readFieldCell(filename, row, col)
    i += 1
  #array = numpy.compress(mask, array)
  return array

def mySelectTArrayWithTimesteps(name, timeSteps, row, col):
  """Selects values at row, col from raster name in time steps.

  name -- Name of raster.
  timeSteps -- Numbers of time steps to use.
  row -- Row index of cell to read, starts with 1!.
  col -- Col index of cell to read, starts with 1!.

  Returned array has elements of type numpy.float32"""
  mask = numpy.zeros(len(timeSteps)).astype(numpy.bool_)
  array = numpy.zeros(len(timeSteps)).astype(numpy.float32)
  i = 0
  while i < len(timeSteps):
    filename = generateNameT(name, timeSteps[i])
    array[i], mask[i] = readFieldCell(filename, row, col)
    i += 1
  #array = numpy.compress(mask, array)
  timesteps = numpy.array(timeSteps)
  outArray=numpy.column_stack((timesteps,array))
  return outArray

def mySelectTArrayForSamples(name, sampleNumbers, timeSteps, row, col):
  """Selects values at row, col from raster name in time steps
  and sample numbers. Returns a numpy array with first column
  time steps and remaining columns map values for each sample."""

  array = numpy.array(timeSteps)
  i = 0
  while i < len(sampleNumbers):
    filename=generateNameS(name,sampleNumbers[i])
    newArray=mySelectTArray(filename, timeSteps, row, col)
    i = i + 1
    array=numpy.column_stack((array,newArray))
  return array

def mySelectTArrayNewNameFormatPercentiles(name, timeSteps, row, col,percentileString):
  # deze kan weg als alles met
  # def mySelectTArrayNewNameFormatPercentilesWithTimesteps(name, timeSteps, row, col,percenti
  # werkt
  """Selects values at row, col from raster name in time steps.

  name -- Name of raster.
  timeSteps -- Numbers of time steps to use.
  row -- Row index of cell to read, starts with 1!.
  col -- Col index of cell to read, starts with 1!.

  #Returned array has elements of type numpy.float32"""
  mask = numpy.zeros(len(timeSteps)).astype(numpy.bool_)
  array = numpy.zeros(len(timeSteps)).astype(numpy.float32)
  i = 0
  while i < len(timeSteps):
    filename = name + "_" + str(timeSteps[i]) + "_" + percentileString + ".map"
    array[i], mask[i] = readFieldCell(filename, row, col)
    i += 1
  return array

def mySelectTArrayNewNameFormatPercentilesWithTimesteps(name, timeSteps, row, col,percentileFloat):
  """Selects values at row, col from raster name in time steps.

  name -- Name of raster.
  timeSteps -- Numbers of time steps to use.
  row -- Row index of cell to read, starts with 1!.
  col -- Col index of cell to read, starts with 1!.

  #Returned array has elements of type numpy.float32"""
  mask = numpy.zeros(len(timeSteps)).astype(numpy.bool_)
  array = numpy.zeros(len(timeSteps)).astype(numpy.float32)
  i = 0
  while i < len(timeSteps):
    filename = name + "_" + str(timeSteps[i]) + "_" + str(percentileFloat) + ".map"
    array[i], mask[i] = readFieldCell(filename, row, col)
    i += 1
  timesteps = numpy.array(timeSteps)
  outArray=numpy.column_stack((timesteps,array))
  return outArray

###########################################
# converting written maps to ascii output #
###########################################

def writeTimeseriesPercentiles(names, percentiles, timeSteps, row, col):
  for name in names:
    for percentile in percentiles:
      array=mySelectTArrayNewNameFormatPercentilesWithTimesteps(name, timeSteps, row, col,percentile)
      filename = name + "_" + str(percentile) + "_r" + str(row) + "_c" + str(col) + ".tss.txt"
      numpy.savetxt(filename,array)

def writeTimeseriesForSamples(names, sampleNumbers, timeSteps, row, col):
  for name in names:
    array=mySelectTArrayForSamples(name, sampleNumbers, timeSteps, row, col)
    filename = name + "_r" + str(row) + "_c" + str(col) + ".tss.txt"
    numpy.savetxt(filename,array)

def writeStaticForSamples(names,sampleNumbers, row,col):
  for name in names:
    array=mySelectSArray(name, sampleNumbers, row, col)
    filename = name + "_r" + str(row) + "_c" + str(col) + ".static.txt"
    numpy.savetxt(filename,array)

def writeTimeseries(names,timeSteps,row,col):
  for name in names:
    array=mySelectTArrayWithTimesteps(name, timeSteps, row, col)
    filename = name + "_r" + str(row) + "_c" + str(col) + ".tss.txt"
    print filename
    numpy.savetxt(filename,array)


#####################
# functions on maps #
#####################

def shiftInTimeSamplesAndTimesteps(inputName,outputName,shiftTimeName,addTime,sampleNumbers,timeSteps, roundingTimestep):
  shiftArray=mySelectSArray(shiftTimeName, sampleNumbers, 1, 1)
  i = 0
  for sample in sampleNumbers:
    timeShiftR=int(round(shiftArray[i]/float(roundingTimestep))*float(roundingTimestep))
    i = i + 1
    for timeStep in timeSteps:
      inputNameString=generateNameST(inputName,sample,timeStep)
      if (timeStep-timeShiftR+addTime) > 0:
        outputNameString=generateNameST(outputName,sample,timeStep-timeShiftR+addTime)
        print inputNameString, outputNameString
        shutil.copy(inputNameString, outputNameString)

#shiftInTimeSamplesAndTimesteps("jan","klaas","piet0000.001",2000,range(1,6,1),range(5,105,5),5)

#################
# renaming maps #
#################

def renamePercentileStackOfMapsToNormal(name,timesteps,percentile):
  for t in timesteps:
    tString=str(t)
    percentileString=str(percentile)
    oldName=name + "_" + tString + "_" + percentileString + ".map"
    newName=generateNameT(name,t)
    shutil.copy(oldName, newName)

##############
## plotting  #
##############

def createTimeseries(name, sampleNumbers, timeSteps, row, col, filename):
  array=mySelectTArrayForSamples(name,sampleNumbers,timeSteps, row, col)
  i = 0
  while i < len(sampleNumbers):
    plt.plot(array[:,0],array[:,i+1])
    i = i + 1
  #plt.axis([2.5,3.5,2.2,18.3])
  plt.savefig(filename)
  plt.gcf().clear()

def createTimeseriesShiftTime(name, shiftTimeName, sampleNumbers, timeSteps, row, col, filename):
  shiftArray=mySelectSArray(shiftTimeName, sampleNumbers, row, col)
  array=mySelectTArrayForSamples(name,sampleNumbers,timeSteps, row, col)
  i = 0
  while i < len(sampleNumbers):
    print i, shiftArray[i]
    plt.plot(array[:,0]-shiftArray[i],array[:,i+1])
    i = i + 1
  plt.savefig(filename)
  plt.gcf().clear()

def createTimeseriesTwoInputsShiftTime(nameOne, nameTwo, shiftTimeName, sampleNumbers, timeSteps, row, col, filename):
  shiftArray=mySelectSArray(shiftTimeName, sampleNumbers, row, col)
  print shiftArray
  arrayOne=mySelectTArrayForSamples(nameOne,sampleNumbers,timeSteps, row, col)
  arrayTwo=mySelectTArrayForSamples(nameTwo,sampleNumbers,timeSteps, row, col)
  i = 0
  while i < len(sampleNumbers):
    print i, shiftArray[i]
    plt.plot(arrayOne[:,0]-shiftArray[i],arrayOne[:,i+1]/arrayTwo[:,i+1])
    i = i + 1
  plt.savefig(filename)
  plt.gcf().clear()

def createTimeseriesConfInt(name, timeSteps, row, col,filename,lowerPercentile,higherPercentile):
  arrayLower=mySelectTArrayNewNameFormatPercentiles(name, timeSteps, row, col,lowerPercentile)
  arrayHigher=mySelectTArrayNewNameFormatPercentiles(name, timeSteps, row, col,higherPercentile)
  arrayMedian=mySelectTArrayNewNameFormatPercentiles(name, timeSteps, row, col,"0.5")
  plt.fill_between(timeSteps,arrayLower, arrayHigher,facecolor="gray",linewidth=0)
  plt.plot(timeSteps,arrayMedian,linewidth=1,color="black")
  ##plt.axis([2.5,3.5,2.2,18.3])
  plt.savefig(filename)
  plt.gcf().clear()



# test
#createTimeseries("jan", [1,2,3],[1,2], 1, 1,"test.pdf")
#createTimeseriesConfInt("1/jan", "1/piet","2/jan", [1,2,3,4,5,6,7,8,9], 1, 1,"test.pdf")
