import sys
import numpy
import datetime
import matplotlib
import matplotlib.pyplot as plt
import generalfunctions
from pcraster import *
#from PCRaster.NumPy import *
from osgeo import gdal
import itertools
import scipy
#import scipy.stats
import scipy.interpolate
from itertools import chain

triu_indices = lambda x: zip(*list(chain(*[[(i, j) for j in range(i, x)] for i in range(x)])))

# 'conversions'

#def timeseriesAsDateTime(numberOfTimeSteps,startTime,timeStepDurationDays):
#  # starTime is datetime 
#  # timeStepDurationDays - floating point
#  steps=numpy.arange(0,numberOfTimeSteps)
#  stepsAsDays=steps*timeStepDurationDays
#  startTimeAsDays=matplotlib.dates.date2num(startTime)
#  realTimeAsDays=stepsAsDays+startTimeAsDays
#  stepsAsDateTime=matplotlib.dates.num2date(realTimeAsDays)
#  return stepsAsDateTime

def timeseriesAsDateTime(numberOfTimeSteps,startTime,timeStepDurationDays):
  # starTime is datetime 
  # timeStepDurationDays - floating point
  steps=numpy.arange(0,numberOfTimeSteps)
  stepsAsDateTime=timeStepsAsDateTime(steps,startTime,timeStepDurationDays)
  return stepsAsDateTime

def timeStepsAsDateTime(steps,startTime,timeStepDurationDays):
  stepsAsDays=steps*timeStepDurationDays
  startTimeAsDays=matplotlib.dates.date2num(startTime)
  realTimeAsDays=stepsAsDays+startTimeAsDays
  stepsAsDateTime=matplotlib.dates.num2date(realTimeAsDays)
  return stepsAsDateTime

def swapXandYInArray(a):
  b=numpy.reshape(a,a.size,order='F').reshape((a.shape[1],a.shape[0])) 
  return b

def timeAverage(data,listWithPeriods,row,col):
  output=[]
  for period in listWithPeriods:
    #dataForPeriod=data[period[0]:period[1],:,0,0] 
    dataForPeriod=data[period[0]:period[1],:,row,col] 
    averageForPeriod=numpy.average(dataForPeriod,axis=0) 
    output.append(averageForPeriod)
  outputAsArray=numpy.array(output)
  return outputAsArray

def scoreAtPercentileOfFlowDurationCurve(timeseries,percentile):
  # NOT USEFUL AT ALL
  # does not work in y direction, you need for each % all samples and
  # these are not available as there is only one bins, i.e. we have
  # realizations on the xaxis, not on the yaxis
  if len(numpy.shape(timeseries)) == 1:  # one sample
    print 'you supplied only one sample'
  fig=plt.figure()
  left=fig.add_subplot(211)
  n,bins,patches=left.hist(timeseries, bins=100,normed=True, cumulative=-1)
  print 'ns are', n
  #score=scipy.stats.scoreatpercentile(n, 50, limit=())
  score=numpy.percentile(n, 50, limit=())
  print score

def valuesInSelectedAreaOfVariablesInStackedList(listOfVariablesSelection):
    """Selects from each variable in listOfVariables the value at the locations
    defined by index. Stacks the results together where each variable will be in
    a row, number of rows in output is number of variables, number of columns
    in output is number of locations in index""" 
    oneRowMatrices=[]
    for selection in listOfVariablesSelection:
      oneRowMatrix=numpy.ravel(selection)
      oneRowMatrices.append(oneRowMatrix)
    stacked=numpy.vstack(oneRowMatrices) 
    return stacked
  

# Axes methods

def plotTimeSeries(self,timeseries,startTime,timeStepDurationDays,timeLoc,timeForm,**kwargs):
  numberOfTimeSteps = numpy.shape(timeseries)[0]
  stepsAsDateTime=timeseriesAsDateTime(numberOfTimeSteps,startTime,timeStepDurationDays)
  self.plot_date(stepsAsDateTime,timeseries, **kwargs)
  self.xaxis.set_major_locator(timeLoc)
  self.xaxis.set_major_formatter(timeForm)

def plotTimeSeriesBars(self,timeseries,startTime,timeStepDurationDays,timeLoc,timeForm,**kwargs):
  numberOfTimeSteps = numpy.shape(timeseries)[0]
  stepsAsDateTime=timeseriesAsDateTime(numberOfTimeSteps,startTime,timeStepDurationDays)
  self.bar(stepsAsDateTime,timeseries)
  self.xaxis.set_major_locator(timeLoc)
  self.xaxis.set_major_formatter(timeForm)

def plotTimeSeriesOfConfidenceInterval(self,timeseriesLower,timeseriesUpper,startTime, \
                                       timeStepDurationDays,timeLoc,timeForm,**kwargs):
  numberOfTimeSteps = numpy.shape(timeseriesLower)[0]
  stepsAsDateTime=timeseriesAsDateTime(numberOfTimeSteps,startTime,timeStepDurationDays)
  self.fill_between(stepsAsDateTime, timeseriesLower, timeseriesUpper,**kwargs)
  self.xaxis.set_major_locator(timeLoc)
  self.xaxis.set_major_formatter(timeForm)

def plotVerticalLinesInTimeSeries(self,timesteps,startTime,timeStepDurationDays):
  stepsAsDateTime=timeStepsAsDateTime(timesteps,startTime,timeStepDurationDays)
  for timestep in stepsAsDateTime:
    plt.axvline(timestep,linestyle=":") 

def interpolateFlowDurationCurve(timeseries,panel):
  n,bins,patches=panel.hist(timeseries, bins=100,normed=True, cumulative=-1)
  xVals=numpy.linspace(0,100,200)
  x=numpy.array([])
  i=0
  for realization in range(0,len(n)):
    nOfRealization = n[realization]
    yVals=numpy.interp(xVals,100.0*nOfRealization[::-1],bins[1:][::-1])
    if i == 0:
      x=yVals
    else:
      x=numpy.vstack((x,yVals))
    discharges=numpy.transpose(x)
    i=i+1
  return n, bins, patches, xVals, discharges

def plotFlowDurationCurve(self,timeseries,**kwargs):
  fig=plt.figure()
  left=fig.add_subplot(211)
  n,bins,patches,xVals,discharges=interpolateFlowDurationCurve(timeseries,left)
  if len(numpy.shape(timeseries)) == 1:  # one sample
    self.plot(n*100.0,bins[1:],**kwargs)
  else: # more than one sample
    self.plot(xVals,discharges,**kwargs)
  self.set_xlim(0,40)
  self.set_xlabel('% time above discharge')
  self.set_ylabel('discharge')

def getQInFlowDuration(percentiel,xVals,median):
  p=(len(xVals)/100.0)*percentiel
  position=int(round(p))
  print 'Q value ', percentiel, 
  print ' discharge ', median[position]

def plotConfidenceIntervalOfFlowDurationCurve(self,timeseries,percentileLower,percentileUpper,**kwargs):
  fig=plt.figure()
  left=fig.add_subplot(211)
  n,bins,patches,xVals,discharges=interpolateFlowDurationCurve(timeseries,left)
  median=numpy.percentile(discharges,50,axis=1)
  lower=numpy.percentile(discharges,percentileLower,axis=1)
  upper=numpy.percentile(discharges,percentileUpper,axis=1)
  self.fill_between(xVals, lower, upper,**kwargs)
  self.plot(xVals,median,color=kwargs['color'])
  self.set_xlim(0,100)
  self.set_xlabel('% time above discharge')
  self.set_ylabel('discharge')
  # print Q5
  getQInFlowDuration(5.0,xVals,median)
  getQInFlowDuration(50.0,xVals,median)
  getQInFlowDuration(95.0,xVals,median)

# Figures

def scatterPlotMatrix(listOfVariablesSelection,names,**kwargs):
  data = valuesInSelectedAreaOfVariablesInStackedList(listOfVariablesSelection)
  fig, correlationMatrix = scatterPlotMatrixOfDataFrame(data, names, **kwargs)
  return fig, correlationMatrix


def scatterPlotMatrixOfDataFrame(data, names, **kwargs):
    """Plots a scatterplot matrix of subplots.  Each row of "data" is plotted
    against other rows, resulting in a nrows by nrows grid of subplots with the
    diagonal subplots labeled with "names".  Additional keyword arguments are
    passed on to matplotlib's "plot" command. Returns the matplotlib figure
    object containg the subplot grid."""
    numvars, numdata = data.shape
    fig, axes = matplotlib.pyplot.subplots(nrows=numvars, ncols=numvars, figsize=(8,8))
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
        if ax.is_last_col():
            ax.yaxis.set_ticks_position('right')
        if ax.is_first_row():
            ax.xaxis.set_ticks_position('top')
        if ax.is_last_row():
            ax.xaxis.set_ticks_position('bottom')

    # Plot the data.

    ## Plot the data. only available in numpy 1.4
    #for i, j in zip(*numpy.triu_indices_from(axes, k=1)):
    #    for x, y in [(i,j), (j,i)]:
    #        axes[x,y].plot(data[x], data[y], **kwargs)

    # work around
    panels = zip(*triu_indices(numvars))
    for i, j in panels:
        for x, y in [(i,j), (j,i)]:
            axes[x,y].plot(data[x], data[y], **kwargs)

    # Label the diagonal subplots...
    for i, label in enumerate(names):
        axes[i,i].annotate(label, (0.5, 0.5), xycoords='axes fraction',
                ha='center', va='center')

    # Turn on the proper x or y axes ticks.
    for i, j in zip(range(numvars), itertools.cycle((-1, 0))):
        axes[j,i].xaxis.set_visible(True)
        axes[i,j].yaxis.set_visible(True)

    correlationMatrix=numpy.corrcoef(data)
    return fig, correlationMatrix


def mapsOfMapTimeSeries(mapsAsNumpy,timesteps,samples,labels):
  '''
  Plots on each row the maps for timesteps, on each column the samples (latter not tested)
  mapsAsNumpy -- time series of maps as numpy
  timesteps -- list
  samples -- list
  labels -- titles of panels, 2D list
  example: theFigure = mapsOfMapTimeSeries(c,[0,1,2],[0],[['1000','2000','3000']])
  Voor plotten van ldd inspireer je door http://matplotlib.org/examples/specialty_plots/hinton_demo.html
  '''
  numberOfCols=len(timesteps) # add one for colorbar
  fig, axes = matplotlib.pyplot.subplots(nrows=len(samples), ncols=numberOfCols,squeeze=False)

  minVal = numpy.min(mapsAsNumpy[timesteps,samples,:,:])
  maxVal = numpy.max(mapsAsNumpy[timesteps,samples,:,:])
  print minVal,maxVal
  a=matplotlib.colors.Normalize(vmin=minVal,vmax=maxVal)

  y=0
  for sample in samples:
    x=0
    for timestep in timesteps:
      data=mapsAsNumpy[timestep,sample,:,:]
      print data
      jan=axes[y,x].imshow(data,interpolation="nearest",norm=a)
      axes[y,x].axes.get_xaxis().set_ticks([])
      axes[y,x].axes.get_yaxis().set_ticks([])
      axes[y,x].set_title(labels[y][x])
      x=x+1
    y=y+1
 
  fig.subplots_adjust(right=0.80) 
  cax = fig.add_axes([0.85, 0.235, 0.045, 0.5])
  fig.colorbar(jan,cax=cax)
  #fig.colorbar(jan,cax=axes[0,numberOfCols-1],fraction=0.1)

  return fig



# helper functions

# moving average

def createBinBoundPairs(slices):
  binBoundPairs=[]
  i = 1
  while i < len(slices):
    pair=[slices[i-1],slices[i]]
    binBoundPairs.append(pair)
    i = i+1
  return binBoundPairs

def maskValuesNotInBin(bin,x,y):
  aboveLowerBound=x > bin[0]
  belowUpperBound=x < bin[1] 
  fallsInBin = aboveLowerBound & belowUpperBound
  xSelected=numpy.where(fallsInBin,x,numpy.zeros(numpy.shape(x))-9999)
  ySelected=numpy.where(fallsInBin,y,numpy.zeros(numpy.shape(x))-9999)
  xMasked=numpy.ma.masked_equal(xSelected,-9999)
  yMasked=numpy.ma.masked_equal(ySelected,-9999)
  return xMasked,yMasked
  
def griddataMean(x,y,slices):
  # x and y are arrays of equal length
  # slices is a range, e.g. slices=numpy.arange(0,1000,100)
  # x and y are binned according to slices and for each bin
  # the mean in x and y is calculated
  # accumulated means are returned
  # 'moving average along the x axis'
  # bio = numpy.arange(0,1000,1) + numpy.random.rand(1000)*100.0
  # growth = numpy.arange(0,1000,1)/100.0+numpy.random.rand(1000)*100.0
  # compressed statement is required because median does not work
  # on compressed arrays..
  binBoundPairs=createBinBoundPairs(slices)
  xOut=[]
  yOut=[]
  xOutAllPercentiles=[]
  yOutAllPercentiles=[]
  #percentiles=[10,20,30,40,50,60,70,80,90]
  percentiles=[20,30,40,50,60,70,80]
  print 'goes wrong when x or y has -9999!!!'
  print 'due to maskValuesNotInBin'
  for bin in binBoundPairs:
    xValues,yValues=maskValuesNotInBin(bin,x,y)
    xValuesCompressed=numpy.ma.compressed(xValues) 
    yValuesCompressed=numpy.ma.compressed(yValues) 
    meanX=numpy.percentile(xValuesCompressed,50)
    meanY=numpy.percentile(yValuesCompressed,50)
    percentilesX=numpy.percentile(xValuesCompressed,percentiles)
    percentilesY=numpy.percentile(yValuesCompressed,percentiles)
    xOut.append(meanX)
    yOut.append(meanY)
    xOutAllPercentiles.append(percentilesX)
    yOutAllPercentiles.append(percentilesY)
  return xOut,yOut,xOutAllPercentiles,yOutAllPercentiles
