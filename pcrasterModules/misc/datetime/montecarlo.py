from pcraster import *
from datetime import *
import datetimePCRasterPython

class SnowModel(object):
  def __init__(self):
    pass

  def premcloop(self):
    self.d_dateTimePCRasterPython=datetimePCRasterPython.DatetimePCRasterPython(datetime(1980,2,3),timedelta(0,10))

  def initial(self):
    pass

  def dynamic(self):
    self.d_dateTimePCRasterPython.update()
    self.d_dateTimePCRasterPython.printit()

  def postmcloop(self):
    pass

myModel = SnowModel()
dynamicModel = DynamicFramework(myModel, 1000)
mcModel = MonteCarloFramework(dynamicModel, 1)
mcModel.run()
