# -*- coding: utf-8 -*-
from PCRaster import * 
import sys
from PCRaster.Framework import *
import generalfunctions,numpy

# class to write a variable as numpy array in PCRaster Python script
# works only for a set of locations
# reportLocations should be called each sample and timestep
# reports for all samples and timesteps
# output is a numpy array of PCRaster maps, where number of rows is 1 and number of columns
# equals the number of non zero values on the locations map

##########################################
# dit werkt alleen bij niet forken!!!!!!
##########################################

class VariableAsNumpy:
  def __init__(self,baseFileName,nrOfSamples,nrOfTimeSteps,nrOfRowsInArray,nrOfColsInArray):
    # if reportLocations is used, nrOfRowsInArray should be 1 and nrOfColsInArray should be equal to nr of cells on locations ne 0
    self.baseFileName=baseFileName
    self.nrOfSamples=nrOfSamples
    self.nrOfTimeSteps=nrOfTimeSteps
    self.nrOfRowsInArray=nrOfRowsInArray
    self.nrOfColsInArray=nrOfColsInArray
    # first index time steps # second index samples # third and fourth coordinates
    self.variableAsArray=numpy.empty(self.nrOfTimeSteps*self.nrOfSamples*self.nrOfRowsInArray*self.nrOfColsInArray) \
                             .reshape(self.nrOfTimeSteps,self.nrOfSamples,self.nrOfRowsInArray,self.nrOfColsInArray) \

  def reportLocations(self,variable,locations,currentSample,currentTimeStep):
    anArray=generalfunctions.returnCellValuesAtLocationsAsNumpyArray(locations,variable)
    self.variableAsArray[currentTimeStep-1,currentSample-1,0,:]=anArray
    if (currentSample == self.nrOfSamples) and (currentTimeStep == self.nrOfTimeSteps):
      numpy.save(self.baseFileName,self.variableAsArray)
