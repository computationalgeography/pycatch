# general
import datetime
import glob
import sys
import os

# add path to modules required
sys.path.append("./pcrasterModules/")

import configuration_seconds as cfg

# PCRaster itself
import pcraster as pcr
import pcraster.framework as pcrfw

# from pcrasterModules
import generalfunctions
import rainfallperiodconstant
import runoffkinematic

# from this folder
import exchangevariables

class CatchmentModel(pcrfw.DynamicModel, pcrfw.MonteCarloModel):
  def __init__(self):
    pcrfw.DynamicModel.__init__(self)
    pcrfw.MonteCarloModel.__init__(self)
    pcr.setclone(cfg.cloneString)

  def premcloop(self):
    self.clone = pcr.boolean(cfg.cloneString)
    self.dem = pcr.scalar(cfg.dem)

  def initial(self):
    # time step duration in hours
    self.timeStepDuration = cfg.timeStepDurationHoursFloatingPointValue
    self.createInstancesInitial()

  def dynamic(self):
    rainfallFlux = self.d_rainfall.update(self.currentTimeStep())
    self.report(rainfallFlux, 'p')

    # discharge calculated with kinematic
    actualAbstractionFromSurfaceWaterFlux, Q = self.d_runoffKinematic.update(rainfallFlux,0.0)
    # discharge in cubic metres per second
    self.report(Q,'q')

    # discharge calculated with accuthreshold 
    QAccu = ((pcr.accuthresholdflux(cfg.lddMap, rainfallFlux,0.0) * pcr.cellarea()) * self.timeStepDuration ) / (self.timeStepDuration * 3600.0)
    # discharge in cubic metres per second
    self.report(QAccu,'qa')

  def postmcloop(self):
    print("no postmcloop functons")

  def createInstancesInitial(self):
    self.startRainStorm=0.1
    self.endRainStorm=0.5
    self.d_rainfall=rainfallperiodconstant.RainfallPeriodConstant(0.04, self.startRainStorm,
                    self.endRainStorm,self.timeStepDuration)

    self.d_runoffKinematic=runoffkinematic.RunoffKinematic(
                           cfg.lddMap, cfg.dem, pcr.celllength(),
                           pcr.scalar(0.02), pcr.scalar(0.6), pcr.scalar(0.00000001),
                           pcr.scalar(0.000000001), self.timeStepDuration)

myModel = CatchmentModel()
dynamicModel = pcrfw.DynamicFramework(myModel, cfg.numberOfTimeSteps)
mcModel = pcrfw.MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
mcModel.setForkSamples(True, 10)
mcModel.run()
