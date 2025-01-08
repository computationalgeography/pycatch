# notes
# time step duration in h

class DatetimePCRasterPython:
  """Maintains time in a PCRaster Python script"""

  def __init__(self, startTimeDatetimeFormat, timeStepDurationDatetimeFormat):
    """Create an instance of the datetime class

    Arguments:
    startTimeDatetimeFormat -- real time at timestep 0, given in the format of the Python datetime module
                               e.g., 26 february, 2001 is datetime(year=2001, month=2, day=26) 
    timeStepDurationDatetimeFormat -- real duration of a timestep, e.g. one hour is timedelta(hours=1)
    
    """
    self.startTimeDatetimeFormat=startTimeDatetimeFormat
    self.currentTimeDatetimeFormat=self.startTimeDatetimeFormat
    self.timeStepDurationDatetimeFormat= timeStepDurationDatetimeFormat
    self.initialTime=self.startTimeDatetimeFormat

  def update(self):
    """Update the time by adding a time duration to the time."""
    self.currentTimeDatetimeFormat=self.currentTimeDatetimeFormat + self.timeStepDurationDatetimeFormat


  def resume(self,timeStep):
    timeSinceStart=self.timeStepDurationDatetimeFormat*timeStep
    self.currentTimeDatetimeFormat=self.initialTime + timeSinceStart

  def getTimeDatetimeFormat(self):
    """Return the current time."""
    return self.currentTimeDatetimeFormat

  # KDJ
  # def __repr__(self):
  #   return "%s %s %s ..." % (self.rainFlux, 'rainfall flux', 'm/h', row, column)

  def printit(self):
    print(self.currentTimeDatetimeFormat)
