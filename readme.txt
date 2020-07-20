User manual
=============

Copyright of this software: Derek Karssenberg & Noemi Lana-Renault Monreal
Contact: d.karssenberg@uu.nl or noemi-solange.lana-renault@unirioja.es 


Installation
================

Create a conda environment by running pcraster_pycatch.yaml

To run the model, run main.py


Configuring the model
=============================

To remove all output, run clean.sh, this does NOT remove inputs, so there is no major risk, but check what is in clean.sh first

Reports are defined in the modules, in pcrasterPythonModules
As PCRaster maps:
  - switch on or off reporting of whole modules, comment/uncomment in reportComponentsDynamic (in main.py)
  - select maps reported by a module, edit in the module component def report
  - note, in module component there is 'full' and 'filtering' by selecting one of these when instantiating
    the module in main.py, you can also make selections in the reporting
As numpy arrays (stored as .txt files and then combined in a numpy array in the postmcloop)
  - switch on or off reporting of whole modules, comment/uncomment in reportComponentsDynamic (in main.py)
  - select maps reported as numpy by a module, edit in the module component, reportAsNumpy
    en reportAsNumpyComponentsPostmcloop


Output file name conventions
----------------------------------

# first letter:

exchange vars starting with                  X
precipitation files starting with            P
interception files starting with             V (from vegetation)
surface store files starting with            S
infiltration files staring with              I
evapotranspiration files starting with       E
subsurface store files starting with         G (from groundwater)
runoff files starting with                   R
shading files starting with                  M (from shading due to Mountains )
budgets files starting with                  B
soil wash files starting with                W
regolith files starting with                 A
bedrock weathering files starting with       C
base level files starting with               L
creep files starting with                    D

# second letter:

store (s, unit m)

actual flux in (i, unit m/h, for geomorphology m/year)
potential flux in (j, unit m/h, for geomorphology m/year)

actual flux out (o, unit m/h, for geomorphology m/year)
potential flux out (p, unit m/h, for geomorphology m/year)

actual flux in (in (positive) or out(negative)) (c, of change, unit m/h, for geomorphology m/year)
another flux (x, unit m/h, for geomorphology m/year)

lateral flux (q, cubic metre per hour)

two extra letters: other values, eg Ecl, cloud factor

see respective classes for details !!  E.g. actual flux out is not always
all out fluxes but sometimes only one of two.



For particle filtering (not yet tested)
-------------------------------------------

To create observed discharge, run createObservedDischarge.sh, it creates maps in the folder 'observations'. This is only required
for particle filtering

todo and changes are in changes.txt
at the bottom of changes is also the names of files used

file name conventions
~~~~~~~~~~~~~~~~~~~~~~

# reports as numpy arrays

Got.npy   self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour, from subsurface
Gxt.npy   self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,from subsurface
Vot.npy   self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour, from canopy
Rq.npy    discharge m3/h
Rqs.npy   discharge m3/h, averaged over 2 hours

RPic.npy  self.maximumInterceptionCapacityPerLAI,
RPks.npy  self.ksat,
RPmm.npy  self.multiplierMaxStomatalConductance
RPrt.npy  self.regolithThicknessHomogeneous,
RPsc.npy  self.saturatedConductivityMetrePerDay,
