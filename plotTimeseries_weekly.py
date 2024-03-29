import matplotlib.pyplot as plt
import numpy

# script to create plots for timeseries output of main_weekly.py

# first index time steps, second index samples, third row number, fourth column number 
#gA = numpy.load('gA.npy')
bioA = numpy.load('bioA.npy')
regA = numpy.load('regA.npy')
#sfA = numpy.load('sfA.npy')       # moisture content
#qA = numpy.load('qA.npy')
#gpA = numpy.load('gpA.npy')      # growth part
#grA = numpy.load('grA.npy')      # grazing
#grNA = numpy.load('grNA.npy')      # net growth
#depA = numpy.load('depA.npy')    # net deposition
#weaA = numpy.load('weaA.npy')    # net weathering
#creA = numpy.load('creA.npy')    # net deposition due to creep

g = numpy.load('1/grazing.npy')


fig=plt.figure()

sample = 0
row = 0      # has to be zero here always

location = 0   # location zero is the one with code zero on the map

bioAMA=bioA[:,[0],row,location]
one=fig.add_subplot(211)
one.plot(bioAMA)

regAMA=regA[:,[0],row,location]
one=fig.add_subplot(212)
one.plot(regAMA)

fig.savefig("plotTimeseries_weekly.pdf",format="pdf")
