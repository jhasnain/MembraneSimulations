import freeEnergyAnalysis as FE
import plotmembranedatafiles as pltmem

toptag=[ 
         ["Sungmin", "Run_4"]
       ]

InDir="/Data/membranesims/DataFiles/"
OutDirs = ["./output/"]*len(toptag)

lstart=500;
binnums=100
maxpoints=1000000

for index, tag in enumerate(toptag):
  OutDir = OutDirs[index]
  SimRun = FE.ObtainData(InDir, OutDir, tag, binnums, lstart, maxpoints)
  FE.PerformFreeEnergyCalculations(SimRun, binnums)
  pltmem.PlotDatafiles(OutDir)
  

