from PCRaster import *
from PCRaster.Framework import *
import numpy, os, operator
##
import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
##
from collections import deque
from PCRaster.NumPy import *
import random

# time in hours

def getCellValue(Map, Row, Column):
  Value, Valid=cellvalue(Map, Row, Column)
  if Valid:
    return Value
  else:
    print 'missing value in input of getCellValue'

def selectACell(Map, XInNrCells, YInNrCells):
  xMap=nominal(xcoordinate(boolean(Map))/celllength())
  yMap=nominal(ycoordinate(boolean(Map))/celllength())
  location=pcrand((xMap == XInNrCells), (yMap == YInNrCells))
  return location
  
def printErrorMessageIfACellContainsTrue(booleanMap,errorMessage):
  scalarMap=cover(scalar(booleanMap),0)
  cellContainsTrueMap=boolean(mapmaximum(scalarMap))
  cellContainsTrue=getCellValue(cellContainsTrueMap,1,1)
  if cellContainsTrue > 0.5:
    print errorMessage

def getCellValueAtBooleanLocation(location,map):
  # map can be any type, return value always float
  valueMap=mapmaximum(ifthen(location,scalar(map)))
  value=getCellValue(valueMap,1,1)
  return value

def printCellValue(self, mapVariable, variableNameToPrint, unit, row, column):
  cellValue=getCellValue(mapVariable,row,column)
  print variableNameToPrint + ' (' + unit + ') at row ' + str(row) + ', column: ' + str(column) + ' is: ' + str(cellValue)

def onePeriod(self, startTime, endTime, timeStepDuration, currentTimeStep):
  # this could be separated in two functions, one converting hours to
  # time steps, one creating the period
  time = float(currentTimeStep) * float(timeStepDuration)
  period = (time > startTime) & (time < endTime)
  return period

def mapeq(mapOne, mapTwo):
  mapOneScalar=scalar(mapOne)
  mapTwoScalar=scalar(mapTwo)
  difference=mapOneScalar-mapTwoScalar
  cellEqual=pcreq(difference,scalar(0))
  mapEqual=pcrgt(mapminimum(scalar(cellEqual)),scalar(0.5))
  return getCellValue(mapEqual,1,1)

def slopeToDownstreamNeighbour(dem, ldd):
  slopeToDownstreamNeighbour=(dem-downstream(ldd,dem))/downstreamdist(ldd)
  return slopeToDownstreamNeighbour

def slopeToDownstreamNeighbourNotFlat(dem,ldd,minSlope):
  slopeToDownstreamNeighbourMap=slopeToDownstreamNeighbour(dem,ldd)
  lddArea=defined(ldd)
  minSlopeCover=ifthen(lddArea,scalar(minSlope))
  slopeToDownstreamNeighbourNotFlat=cover(max(minSlopeCover,slopeToDownstreamNeighbourMap), minSlopeCover)
  return slopeToDownstreamNeighbourNotFlat

def distancetodownstreamcell(Ldd):
  distanceToDownstreamCell=max(downstreamdist(Ldd),celllength())
  return distanceToDownstreamCell

def createTimeSeriesList(timeSeriesFile):
  file = open(timeSeriesFile,'r')
  piet = file.readlines()
  newList=[]
  for line in piet:
    lineList=string.split(line)
    newList.append(lineList)
  file.close()
  return newList

def timeInputSparse(fileName):
  return os.path.exists(fileName)

def normalcorrelated(normalX, normalY, correlation):
  # returns realizations of two normal variables with
  # mean zero and var 1 having correlation of correlation
  # based on:
  # x=normal()
  # y=ax+b*normal()
  # correlation = a /  sqrt( sqr(a) + sqr(b) ) 
  x=scalar(normalX)
  y=(x+sqrt((1/sqr(correlation))-1)*scalar(normalY)) * scalar(correlation)
  return x,y


def mapToColAsArray(name):
  """Selects values at row, col from raster name in Monte Carlo samples.

  name -- Name of raster.
  row -- Row index of cell to read.
  col -- Col index of cell to read.
  The returned array does not contain missing values so the size is maximimal
  the number of cells. It contains three columns, x, y, name
  x,y are given as xcoordinate and ycoordinate values

  Returned array has elements of type numpy.float32"""
  nrRows=clone().nrRows()
  nrCols=clone().nrCols()
  nrCells=nrRows*nrCols

  mask = numpy.zeros(nrCells).astype(numpy.bool_)
  arrayX = numpy.zeros(nrCells).astype(numpy.float32)
  arrayY = numpy.zeros(nrCells).astype(numpy.float32)
  arrayName = numpy.zeros(nrCells).astype(numpy.float32)

  xMap=xcoordinate(defined(name))
  yMap=ycoordinate(defined(name))

  # For each cell.
  c = 0
  while c < nrCells:
    arrayName[c], mask[c] = cellvalue(name, c + 1)
    arrayX[c], dummy = cellvalue(xMap, c + 1)
    arrayY[c], dummy = cellvalue(yMap, c + 1)
    c += 1

  arrayName = numpy.compress(mask, arrayName)
  arrayX = numpy.compress(mask, arrayX)
  arrayY = numpy.compress(mask, arrayY)
  mapAsColArray = numpy.column_stack((arrayX,arrayY,arrayName))
  return mapAsColArray 

def addTimeColumnToMapAsColArray(mapAsColArray,time):
  b = numpy.insert(mapAsColArray,2,float(time), axis=1)
  return b

def stackOfMapsToColAsArray(stackOfMapsAsList, currentTime):
  timeOfFirstMap=currentTime-len(stackOfMapsAsList)+1
  t = timeOfFirstMap
  stackOfMapsAsColsList=[] 
  for map in stackOfMapsAsList:
    mapArray = mapToColAsArray(map)
    mapAsColArrayWithTime = addTimeColumnToMapAsColArray(mapArray,t)
    stackOfMapsAsColsList.append(mapAsColArrayWithTime)
    t =  t + 1
  array = numpy.concatenate(stackOfMapsAsColsList, axis=0)
  return array

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
  convertToDataFrame=robjects.r['convertToDataFrame']
  a = convertToDataFrame(mapAsColArray)
  return a

def experimentalVariogramValues(stackOfMapsAsRDataFrame,boundariesVector,space,savePlot,fileName,maxVarPlot):
  # returns distances and semivariances for steps defined by boundaries vector
  # note that length of returned vector equals number of intervals that is
  # available, thus, len(returnedList) can be smaller than len(boundariesVector)!
  # space (TRUE) -> spatial correlation
  # space (FALSE) -> temporal correlation
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
  experimentalVariogramValues=robjects.r['experimentalVariogramValues']
  boundariesVectorR=robjects.FloatVector(boundariesVector)
  expVariogram = experimentalVariogramValues(stackOfMapsAsRDataFrame,boundariesVectorR)
  if savePlot:
    robjects.r('''
    saveExperimentalVariogram <- function(experimentalVariogram,fileName,maxVarPlot) {
        require("gstat")
        png(fileName)
        plot(experimentalVariogram$dist,experimentalVariogram$gamma,ylim=c(0,maxVarPlot))
        dev.off()
        }
        ''')
    saveExperimentalVariogram=robjects.r['saveExperimentalVariogram']
    saveExperimentalVariogram(expVariogram,fileName,maxVarPlot)
  #return expVariogram.r['dist'][0], expVariogram.r['gamma'][0]
  return expVariogram[1], expVariogram[2]

def semvar(firstMap,secondMap):
  nrPairs=getCellValue( cover(maptotal(scalar(pcrand(defined(firstMap),defined(secondMap)))),0)  ,1,1)
  sumOfSquaredDiff=getCellValue( cover(maptotal(sqr(firstMap-secondMap)/2.0),0)   ,1,1)
  return nrPairs,sumOfSquaredDiff

def experimentalVariogramValuesInTime(stackOfMapsAsList,bounds):
  nrPairsOfLags=[0.0]*len(bounds)
  print nrPairsOfLags
  sumOfSquaredDiffOfLags=[0.0]*len(bounds)
  sumOfDists=[0.0]*len(bounds)
  nMaps=len(stackOfMapsAsList)
  for i in range(0,nMaps):
    for j in range(i+1,nMaps):
      dist=math.fabs(i-j)
      nrPairs,sumOfSquaredDiff=semvar(stackOfMapsAsList[i],stackOfMapsAsList[j])
      k=0
      used=0
      for bound in bounds:
        if (dist < bound) and (used==0):
          nrPairsOfLags[k]=nrPairsOfLags[k]+nrPairs
          sumOfSquaredDiffOfLags[k]=sumOfSquaredDiffOfLags[k]+sumOfSquaredDiff
          sumOfDists[k]=sumOfDists[k]+nrPairs*float(dist)
          used=1
        k=k+1
  semvarList=map(operator.div,sumOfSquaredDiffOfLags,nrPairsOfLags)
  distList=map(operator.div,sumOfDists,nrPairsOfLags)
  return distList,semvarList

#jan=ifthen(pcrle(uniqueid(defined('jet00000.001')),100),scalar('jet00000.001'))
#test=[scalar('jet00000.001'),jan, scalar('jet00000.003'),scalar('jet00000.004')]
#a=experimentalVariogramValuesInTime(test,[1.5,7.0])

def semvarOfStackOfMapsInSpace(stackOfMapsAsList,lagX,lagY):
  # lagX is shift of cells to right (positive)
  # lagY is shift of cells up (positive)
  nrPairsTot=0.0
  sumOfSquaredDiffTot=0.0
  for theMap in stackOfMapsAsList:
    shiftedMap=shift(theMap,lagY,0-lagX)
    nrPairs,sumOfSquaredDiff=semvar(theMap,shiftedMap)
    nrPairsTot=nrPairsTot+nrPairs
    sumOfSquaredDiffTot=sumOfSquaredDiffTot+sumOfSquaredDiff
  return nrPairsTot, sumOfSquaredDiffTot

def createLagXAndLagYForBounds(bounds):
  lagXlagYDist=[]
  maxLag=int(round(bounds[-1]))
  possibleLags=range(0,maxLag,1)
  for i in possibleLags:
    for j in possibleLags:
      dist= math.sqrt((float(i)**2)+(float(j)**2))
      if (dist < bounds[-1]) and not ((i==0) and (j==0)):
        a=[ i, j,dist]
        lagXlagYDist.append(a)
  return lagXlagYDist

#b=createLagXAndLagYForBounds([0.2,4.0])
#print b
  

def experimentalVariogramValuesInSpace(stackOfMapsAsList,bounds):
  nrPairsOfLags=[0.0]*len(bounds)
  sumOfSquaredDiffOfLags=[0.0]*len(bounds)
  sumOfDists=[0.0]*len(bounds)
  nMaps=len(stackOfMapsAsList)
  lagXAndLagY=createLagXAndLagYForBounds(bounds)
  for i in lagXAndLagY:
    dist=i[2]
    nrPairs,sumOfSquaredDiff=semvarOfStackOfMapsInSpace(stackOfMapsAsList,i[0],i[1])
    k=0
    used=0
    for bound in bounds:
      if (dist < bound) and (used==0):
        nrPairsOfLags[k]=nrPairsOfLags[k]+nrPairs
        sumOfSquaredDiffOfLags[k]=sumOfSquaredDiffOfLags[k]+sumOfSquaredDiff
        sumOfDists[k]=sumOfDists[k]+nrPairs*float(dist)
        used=1
      k=k+1
  semvarList=map(operator.div,sumOfSquaredDiffOfLags,nrPairsOfLags)
  distList=map(operator.div,sumOfDists,nrPairsOfLags)
  return distList,semvarList

#jan=ifthen(pcrle(uniqueid(defined('jet00000.001')),100),scalar('jet00000.001'))
##test=[scalar('jet00000.001'),jan, scalar('jet00000.003'),scalar('jet00000.004')]
#test=[scalar('jet00000.001'),scalar('jet00000.002'), scalar('jet00000.003'),scalar('jet00000.004')]
#a,b=experimentalVariogramValuesInSpace(test,[2.3,15,17.8])
#print a, b

def descriptiveStatistics(stackOfMapsAsRDataFrame):
  robjects.r('''
  descriptiveStatistics <- function(dataFrame) {
      mean <- mean(dataFrame$v)
      variance <- var(dataFrame$v)
      c(mean,variance)
      }
      ''')
  descriptiveStatistics=robjects.r['descriptiveStatistics']
  var = descriptiveStatistics(stackOfMapsAsRDataFrame)
  return var

def cornerMap(cloneMap):
  x=xcoordinate(cloneMap)
  y=ycoordinate(cloneMap)
  corner= pcrand(pcreq(x, mapminimum(x)),pcreq(y,mapmaximum(y)))
  return corner

def convertListOfValuesToListOfNonSpatialMaps(listOfValues,cloneMap):
  # non spatial, i.e. value is put in the corner with the largest x and y
  # return values is always scalar
  corner=cornerMap(cloneMap)
  listOfMaps=[]
  for value in listOfValues:
    map=ifthen(corner,scalar(value))
    listOfMaps.append(map)
  return listOfMaps

def convertListOfNonSpatialMapsToListOfValues(listOfMaps):
  listOfValues=[]
  clone=ifthenelse(defined(listOfMaps[0]),boolean(1),1)
  x=xcoordinate(clone)
  y=ycoordinate(clone)
  corner= pcrand(pcreq(x, mapminimum(x)),pcreq(y,mapmaximum(y)))
  for map in listOfMaps:
    value=getCellValueAtBooleanLocation(corner,map)
    listOfValues.append(value)
  return listOfValues
    
## tests
#setclone("cloneSmall.map")
#a = scalar("norSmall.map")
#b = scalar("norSmall.map")*2
#c = scalar("norSmall.map")*3
#d = scalar("norSmall.map")*4
#e = scalar("norSmall.map")*5
#stackOfMapsAsList=[a,b,c,d,e]
#report(d,"testje")
#test = convertListOfNonSpatialMapsToListOfValues(stackOfMapsAsList)
#c = stackOfMapsToColAsArray(stackOfMapsAsList, 10)
#d = convertStackOfMapsToRDataFrame(c)
#boundVector=(1.5,2.5,30.5)
#dist,gamma=experimentalVariogramValues(d,boundVector,1,1,'pietje.pdf')
#print dist
#print gamma
#aListOfMaps=convertListOfValuesToListOfNonSpatialMaps(gamma,"cloneSmall.map")
#report(aListOfMaps[0],"testje")
#
#descrStats = descriptiveStatistics(d)
#bListOfMaps=convertListOfValuesToListOfNonSpatialMaps(descrStats,"cloneSmall.map")
#report(bListOfMaps[0],"mean")
#
#import time
#time.sleep(100)

def keepHistoryOfMaps(currentHistoryOfMaps,mapOfCurrentTimeStep,numberOfTimeStepsToKeep):
  # uses deque objects instead of lists, thus, conversion is required to get a list simply
  # by list(currentHistoryOfMaps)
  currentHistoryOfMaps.append(mapOfCurrentTimeStep)
  if len(currentHistoryOfMaps) > numberOfTimeStepsToKeep:
    currentHistoryOfMaps.popleft()
  if len(currentHistoryOfMaps) > numberOfTimeStepsToKeep:
    print 'warning: length of keepHistoryOfMaps is greater than number of timesteps to keep'
  return currentHistoryOfMaps

def nrCols(map):
  x=xcoordinate(boolean(map))
  xMax=mapmaximum(x)
  xMin=mapminimum(x)
  nrCols=((xMax-xMin)/celllength())+1
  return nrCols

def nrRows(map):
  y=ycoordinate(boolean(map))
  yMax=mapmaximum(y)
  yMin=mapminimum(y)
  nrCols=((yMax-yMin)/celllength())+1
  return nrCols

def nrCells(map):
  return nrCols(map)*nrRows(map)

def corners(map):
  left=edge(map,4,0)
  right=edge(map,6,0)
  top=edge(map,8,0)
  bottom=edge(map,2,0)
  return pcrgt(scalar(left)+scalar(right)+scalar(top)+scalar(bottom),1.5)

def edges(map):
  left=edge(map,4,0)
  right=edge(map,6,0)
  top=edge(map,8,0)
  bottom=edge(map,2,0)
  return pcror(pcror(left,right), pcror(top,bottom))
  
def booleanTrue(map):
  # returns a map that is everywhere true
  # removes mvs
  noMVs=cover(map,1)
  return defined(noMVs)

def edge(map,side,distance):
  # map should have y incr. bot to top
  # returns boolean map with edge
  # distance defines distance to map boundary of edge,
  # e.g. distance 0 returns the real edge, distance
  # distance is an integer or nominal or ordinal
  # 1 returns one to the left/right/top/bottom
  # side is ldd dirs, e.g. 4 returns left side
  # side is an integer
  # works only with whole coordinates (as cells are selected using equals on coors as floating points)
  realDist=celllength()*scalar(distance)
  if ((side == 4) or (side == 6)):
    x=xcoordinate(booleanTrue(map))
    if side == 4:
      sideMap=pcreq(x,mapminimum(x)+realDist)
    if side == 6:
      sideMap=pcreq(x,mapmaximum(x)-realDist)
  if ((side == 8) or (side == 2)):
    y=ycoordinate(booleanTrue(map))
    if side == 2:
      sideMap=pcreq(y,mapminimum(y)+realDist)
    if side == 8:
      sideMap=pcreq(y,mapmaximum(y)-realDist)
  return sideMap

def bottom(map):
  '''
  returns the bottom line of cells
  works only with y increases bottom to top
  any cell size (unlike edge)
  '''
  yCoordinate=ycoordinate(defined(map))
  bottom = (yCoordinate == mapminimum(yCoordinate))
  return bottom

def neighbourIsMissingComponent(nbShift,noMVOnMap):
  NBIsMissingAll=ifthenelse(defined(nbShift),boolean(0),boolean(1))
  NBIsMissing=pcrand(NBIsMissingAll,noMVOnMap)
  return NBIsMissing

def neighbourIsMissingValueOrEdgeAndCellItselfIsDefined(map):
  noMVOnMap=defined(map)
  rightNBIsMissing=neighbourIsMissingComponent(shift(map,scalar(0),scalar(1)),noMVOnMap)
  leftNBIsMissing=neighbourIsMissingComponent(shift(map,scalar(0),scalar(-1)),noMVOnMap)
  upperNBIsMissing=neighbourIsMissingComponent(shift(map,scalar(-1),scalar(0)),noMVOnMap)
  lowerNBIsMissing=neighbourIsMissingComponent(shift(map,scalar(1),scalar(0)),noMVOnMap)
  return upperNBIsMissing,rightNBIsMissing,lowerNBIsMissing,leftNBIsMissing


#aMap=scalar("idOth.map")
#aMap=ifthen(uniform(1) < 0.9,scalar(2))
#one,two,three,four=neighbourIsMissingValueOrEdgeAndCellItselfIsDefined(aMap)
#report(one,'one.map')
#report(two,'two.map')
#report(three,'three.map')
#report(four,'four.map')
#report(aMap,'amap.map')


def moveRowsOrColumnsForPeriodicBoundaryCondition(map,direction):
  # moves row/column next to edge cells to other edge
  # direction is ldd dirs (e.g. 8 is bottom to top)
  # direction is an integer
  if (direction == 6):
    cellsToShift=edge(map,4,1)
    distanceToShift=scalar(nrCols(map)-2)
    shiftedMap=ifthen(edge(map,6,0),shift(map,0,0-distanceToShift))
  if (direction == 4):
    cellsToShift=edge(map,6,1)
    distanceToShift=scalar(nrCols(map)-2)
    shiftedMap=ifthen(edge(map,4,0),shift(map,0,distanceToShift))
  if (direction == 2):
    cellsToShift=edge(map,8,1)
    distanceToShift=scalar(nrRows(map)-2)
    shiftedMap=ifthen(edge(map,2,0),shift(map,0-distanceToShift,0))
  if (direction == 8):
    cellsToShift=edge(map,2,1)
    distanceToShift=scalar(nrRows(map)-2)
    shiftedMap=ifthen(edge(map,8,0),shift(map,distanceToShift,0))
  return shiftedMap

def periodicBoundaryCondition(map):
  # note positive shifts are to left or to top..
  # first value is vertical shift, second value is horizontal shift
  left=moveRowsOrColumnsForPeriodicBoundaryCondition(map,4)
  right=moveRowsOrColumnsForPeriodicBoundaryCondition(map,6)
  top=moveRowsOrColumnsForPeriodicBoundaryCondition(map,8)
  bottom=moveRowsOrColumnsForPeriodicBoundaryCondition(map,2)
  tmp=cover(left,right,top,bottom,map)
  newMap=ifthen(pcrnot(corners(map)),tmp)
  return newMap

#test=scalar("idOth.map")
##nrCols=nrCols(test)
##report(nrCols,"test")
#testje=edge(test,8,1)
#report(testje,"testje")
#
#test2=moveRowsOrColumnsForPeriodicBoundaryCondition(test,8)
#report(test2,"test2")
#test3=periodicBoundaryCondition(test)
#report(test3,"test3")
#test4=edges(test)
#report(test4,"test4")

def periodicBoundaryConditionNumpy(map):
  a=pcr2numpy(map,1)
  #second left edge
  secondLeft = a[:,1]
  # second right edge
  secondRight = a[:,a.shape[1]-2]
  # remove left and right edge
  new = numpy.delete(a,(0,a.shape[1]-1),1)
  # add new left edge
  newLeft = numpy.insert(new,0,secondRight,1)
  # add new right edge
  b = numpy.insert(newLeft,newLeft.shape[1],secondLeft,1)
  #second upper edge
  secondUpper = b[1,:]
  # second right edge
  secondLower = b[a.shape[0]-2,:]
  # remove upper and lower edge
  bNew = numpy.delete(b,(0,a.shape[0]-1),0)
  # add upper edge
  newTop = numpy.insert(bNew,0,secondLower,0)
  # add new lower edge
  c = numpy.insert(newTop,newTop.shape[0],secondUpper,0)
  outMap=numpy2pcr(Scalar,c,0)
  return outMap

#test=scalar("idOth.map")
#test2=periodicBoundaryConditionNumpy(test)
#report(test2,"test2")
#test3=periodicBoundaryCondition(test)
#report(test3,"test3")

def createToCellsPeriodicBoundaryCondition(clone):
  list=[]
  list.append(edge(clone,6,0))
  list.append(edge(clone,8,0))
  list.append(edge(clone,4,0))
  list.append(edge(clone,2,0))
  return list

def createFromCellsPeriodicBoundaryCondition(clone):
  list=[]
  list.append(edge(clone,4,1))
  list.append(edge(clone,2,1))
  list.append(edge(clone,6,1))
  list.append(edge(clone,8,1))
  return list

def colNumber(clone):
  return ordinal((roundoff((xcoordinate(clone)/celllength()))))

def rowNumber(clone):
  return ordinal((roundoff((ycoordinate(clone)/celllength()))))

def periodicBoundaryConditionAreatotal(map,fromCells,toCells,cols,rows):
  right=ifthen(toCells[0],areatotal(ifthen(fromCells[0],map),rows))
  bottom=ifthen(toCells[1],areatotal(ifthen(fromCells[1],map),cols))
  left=ifthen(toCells[2],areatotal(ifthen(fromCells[2],map),rows))
  top=ifthen(toCells[3],areatotal(ifthen(fromCells[3],map),cols))
  return cover(right,bottom,left,top,map)


#test=scalar("idOth.map")
#clone=defined(test)
#
#fromCells=createFromCellsPeriodicBoundaryCondition(clone)
#toCells=createToCellsPeriodicBoundaryCondition(clone)
#cols=colNumber(clone)
#rows=rowNumber(clone)
#test2=periodicBoundaryConditionAreatotal(test,fromCells,toCells,cols,rows)
#report(test,"test")
#report(test2,"test2")
#test3=periodicBoundaryCondition(test)
#report(test2-test3,"diff")
  

def samplingScheme(clone,nrSamples,fractionShortDistance,separationDistance,nrCellsToRight,nrCellsToTop):
  # get numb. of samples
  nrSamplesGrid=rounddown((1.0-fractionShortDistance)*nrSamples)
  nrSamplesShortDistance=nrSamples-nrSamplesGrid
  # get possible locs
  colnumber=roundoff((xcoordinate(clone)/celllength()))
  tmp=roundoff((ycoordinate(clone)/celllength()))
  rownumber=tmp-mapminimum(tmp)+1
  possibleCols=pcreq(pcrmod(colnumber,separationDistance),0)
  possibleRows=pcreq(pcrmod(rownumber,separationDistance),0)
  possibleLocations=pcrand(possibleCols, possibleRows)
  # get grid area
  nrRowsCols=roundup(sqrt(nrSamplesGrid))
  locColSel=ifthenelse(pcrle(colnumber,(separationDistance*nrRowsCols)),possibleLocations,0)
  locRowSel=ifthenelse(pcrle(rownumber,(separationDistance*nrRowsCols)),possibleLocations,0)
  locsArea=pcrand(locColSel,locRowSel)
  # get samples
  samplesAtGrid=cover(pcrle(order(ifthen(locsArea,uniqueid(clone))),nrSamplesGrid),0)
  # get samples with extra samples
  randomSampleNumbers=ordinal(ifthen(samplesAtGrid,order(uniform(ifthen(samplesAtGrid,boolean(1))))))
  alreadyCovered=defined(randomSampleNumbers)
  for sample in range(1,int(getCellValue(mapmaximum(randomSampleNumbers),1,1))+1):
    theSample=cover(pcreq(sample,randomSampleNumbers),0)
    nbSample=pcrand(pcrgt(window4total(scalar(theSample)),0.5), pcrnot( alreadyCovered))
    extraSample=cover(pcrlt(order(ifthen(nbSample,uniform(1))),1.5),0)
    alreadyCovered=pcror(alreadyCovered,extraSample)
    totalSamples=maptotal(scalar(alreadyCovered))
    #print getCellValue(totalSamples,1,1)
    if getCellValue(totalSamples,1,1) > (nrSamples-1.0+0.1):
      break
  sampleNumbers=ordinal(ifthen(alreadyCovered,order(uniqueid(ifthen(alreadyCovered,boolean(1))))))
  # shift the cells to centre
  centreCol=maptotal(colnumber)/nrCells(clone)
  centreRow=maptotal(rownumber)/nrCells(clone)
  centreColSamples=maptotal(ifthen(alreadyCovered,colnumber))/maptotal(ifthen(alreadyCovered,scalar(1)))
  centreRowSamples=maptotal(ifthen(alreadyCovered,rownumber))/maptotal(ifthen(alreadyCovered,scalar(1)))
  sampleNumbersShifted=shift(sampleNumbers,centreRow-centreRowSamples-nrCellsToTop,centreColSamples-centreCol-nrCellsToRight)
  return cover(sampleNumbersShifted,0)

def samplingSchemeRandomShift(clone,nrSamples,fractionShortDistance,separationDistance,maxNrCellsToRight,maxNrCellsToTop):
  uniformRight=mapuniform()
  shiftRight=roundoff((uniformRight-0.5)*maxNrCellsToRight)
  uniformTop=mapuniform()
  shiftTop=roundoff((uniformTop-0.5)*maxNrCellsToTop)
  sampleNumbersShifted=samplingScheme(clone,nrSamples,fractionShortDistance,separationDistance,shiftRight,shiftTop)
  return sampleNumbersShifted
  
#clone=defined("clone.map")
#test=samplingSchemeRandomShift(clone,100,0.3,3,4,40)
#test2=samplingScheme(clone,100,0.3,3,0,0)
#report(test,"pietje")
#report(test2,"pietje2")

def createGifAnimation(name,samples):
  # postmcloop function to convert individual png or gif
  # files to animated gif
  # use gifview (from gifsicle package) to animate, man gifview for
  # options
  # or gimp, filters -> animation -> playback
  for sample in samples:
    baseName = generateNameS(name,sample)
    command = "convert " + baseName + "* " + baseName + "_ani.gif"
    os.system(command)
    # alternatief: for i in ete*; do convert $i ${i%%png}gif; done
    # en dan gifsicle to make the animation (maybe faster)

def mapaverage(aMap):
  total=maptotal(aMap)
  nrCellsNoMV=maptotal(scalar(defined(aMap)))
  return total/nrCellsNoMV

########################################
### FUNCTIONS FOR PARTICLE FILTERING ###
########################################

def createFileNameToReportAVariableForSuspend(classInstance,variableName,currentTimeStep,currentSampleNumber):
  className = classInstance.__class__.__name__
  classNameDir = str(currentSampleNumber) + '/stateVar/' + className
  if not os.path.exists(classNameDir):
    os.mkdir(classNameDir)
  return generateNameST('stateVar/' + className + '/' + variableName,currentSampleNumber,currentTimeStep)
  #print className + "/" variableName

def letterSequence():
  alphabeth=string.ascii_letters[0:26]
  letters=[]
  for firstLetter in alphabeth:
    for secondLetter in alphabeth:
      letters.append(firstLetter+secondLetter)
  return letters

def reportMemberVariablesOfAClassForSuspend(classInstance,currentTimeStep,currentSampleNumber):
  b=vars(classInstance)
  names = letterSequence()
  i=0
  for name in sorted(b):
    if type(b[name]) == PCRaster._PCRaster.Field:
      fileName = createFileNameToReportAVariableForSuspend(classInstance,names[i],currentTimeStep,currentSampleNumber)
      report(b[name],fileName)
    i = i+1

def readMemberVariablesOfAClassForResume(classInstance,currentTimeStep,currentSampleNumber):
  b=vars(classInstance)
  names = letterSequence()
  i=0
  for name in sorted(b):
    if type(b[name]) == PCRaster._PCRaster.Field:
      fileName = createFileNameToReportAVariableForSuspend(classInstance,names[i],currentTimeStep,currentSampleNumber)
      mapToResume = readmap(fileName)
      vars(classInstance)[name]  = mapToResume 
    i = i+1

def printMemberVariables(classInstance):
  a = sorted(vars(classInstance))
  print classInstance
  print 'number of member variables ', len(a)
  print a
  print



##############################################
### FUNCTIONS FOR GENERATING RANDOM FIELDS ###
##############################################

def mapuniformBounds(min,max,fixedValue='aString..',useRealization=True):
  """Assigns a random value taken from a uniform distribution

  min -- lower bound, python or pcraster type
  max -- upper bound of interval, python or pcraster type
  fixedValue -- value assigned when useRealization is False
  useRealization -- switch to assign fixedValue

  fixedValue should be a PCRaster type when output is used in PCRaster functions"""
  if useRealization:
    a = mapuniform()
    range = scalar(max) - scalar(min)
    field = (a * range) + min
    return field
  else:
    return fixedValue

def areauniformBounds(min,max,areaMap,fixedValue='aString..',useRealization=True):
  """Assigns a random value taken from a uniform distribution

  min -- lower bound, python or pcraster type
  max -- upper bound of interval, python or pcraster type
  areamap -- classified map
  fixedValue -- value assigned when useRealization is False
  useRealization -- switch to assign fixedValue

  fixedValue should be a PCRaster type when output is used in PCRaster functions"""
  if useRealization:
    a = areauniform(spatial(areaMap))
    range = scalar(max) - scalar(min)
    field = (a * range) + min
    return field
  else:
    return fixedValue

def mapNormalRelativeError(input,standardDeviationForInputOfOne):
  a = mapnormal() 
  error = a * input * standardDeviationForInputOfOne
  realization = input + error
  return realization

def mapgamma(shapeParameter):
  '''Returns a realization from the gamma distribution with a mean of one
  shapeParameter is a Python floating point
  return value is a Python floating point
  '''

  scaleParameter=1.0/shapeParameter
  realization=random.gammavariate(shapeParameter,scaleParameter)
  return realization 

#a=0.0
#for i in range(1,100000):
#  test = mapgamma(6.0)
#  a=a+test
#print 'gamma real :', test
#print 'gamma mean :', a/100000.0



###############################
### FUNCTIONS FOR REPORTING ###
###############################

def reportAListOfMaps(listOfMaps,baseName,timeStep,sample):
  alphabeth=string.ascii_letters[0:26]
  i=0
  for map in listOfMaps:
    totBaseName=baseName+str(alphabeth[i])
    reportName=generateNameST(totBaseName,sample,timeStep)
    report(map,reportName)
    i=i+1
  
#####################################
### FUNCTIONS FOR FRAGSTATS STUFF ###
#####################################

def proportionOfClassInAreas(areas,classMap,selectedClass):
  '''
  Calculates the fractional area of a class on classMap
  for each area on areas map. The class on classMap
  is given by selectedClass.
  '''
  areaOfAreas = areatotal(spatial(scalar(1)),areas)
  selectedClassMap = scalar((classMap == selectedClass))
  proportionOfClassInAreasMap = areatotal(selectedClassMap,areas)/areaOfAreas
  return proportionOfClassInAreasMap

def convertValuesInAreasToSeparateMaps(areas,numberOfAreas,valuesInAreas):
  '''
  Takes the maximum value of valuesInAreas in each area on areas and
  stores it in a map. These maps are collected and returned as a list
  of maps.
  '''
  separateMaps=[]
  for i in range(1,numberOfAreas+1):
    separateMap=mapmaximum(ifthen(areas == i,valuesInAreas))
    separateMaps.append(separateMap)
  return separateMaps

def selectFromEachAreaOneCell(areas):
  '''
  Selects from each class on areas one cell and returns
  the class value from area at that cell. Remaining cells
  have zero value at result.
  '''
  uniqueId=uniqueid(defined(areas))
  oneCellPerArea=ifthenelse(areamaximum(uniqueId,areas) == uniqueId,areas,0)
  return oneCellPerArea

def patchSize(areas):
  areasClumped=clump(areas)
  nrCellsPerPatch=areatotal(spatial(scalar(1)),areasClumped)
  singleCellPerArea=selectFromEachAreaOneCell(areasClumped) != 0
  nrCellsPerPatchStoredOnOneCellPerPatch=ifthen(singleCellPerArea,nrCellsPerPatch)
  meanPatchSize=areaaverage(nrCellsPerPatchStoredOnOneCellPerPatch,areas)
  return meanPatchSize


#areas=(nominal('cl000000.014'))
#report(areas,'areas.map')
#test=selectFromEachAreaOneCell(areas)
#report(test,'test.map')
#patchsize=patchSize(areas)
#report(patchsize,'patchsize.map')
  
