from PCRaster import *
from PCRaster.Framework import *
from generalfunctions import *
from aggregationfunctionsDK import *

setclone("clone.map")

class MayModel(object):
  def __init__(self):
    pass

  def premcloop(self):
    pass

  def initial(self):
    self.clone="clone.map"
    henk=normal(1)*10+100
    self.report(henk,"henk")

  def dynamic(self):
    piet=normal(1)*40
    self.report(piet,"piet")
    jan=normal(1)
    self.report(jan,"jan")
    report(jan,generateNameT("jet",self.currentTimeStep()))

  def postmcloop(self):
    names=['piet','jan']
    percentiles=[0.1,0.5,0.9]
    sampleNumbers=self.sampleNumbers()
    timeSteps=self.timeSteps()
    row=1
    col=1
#    mcCovarMatrix(names,sampleNumbers,timeSteps,row,col,'test')
    mcpercentiles(names,percentiles,sampleNumbers,timeSteps)
    writeTimeseriesPercentiles(names, percentiles, timeSteps, row, col)
    writeTimeseriesForSamples(names, sampleNumbers, timeSteps, row, col)
    namen=['jet']
    writeTimeseries(namen,timeSteps,row,col)
    writeStaticForSamples(["henk"],sampleNumbers, 1,1)

myModel = MayModel()
dynamicModel = DynamicFramework(myModel, 10)
mcModel = MonteCarloFramework(dynamicModel, 25)
mcModel.run()
