import sys
sys.path.append("../datanalysis/")
import os 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pickle
from membranesimhelperfuncs import *
from effectivepotentials import *
from scipy import integrate


def Outname(Params, prefix, suffix):return Params.OutDir+'/'+prefix+'_'+Params.filename+'_'+suffix

def GetValFromStr(string, keyword): return string.split(keyword+'_')[1].split('_')[0]

def ListUniquePercentages(InDir):
  wLlist      = filter(lambda x: 'Data' in x, os.listdir(InDir+'/'))
  preamb      = wLlist[0].split('percent')[0]  
  percentlist = list(set( [GetValFromStr( wLfile, 'percent' ) for wLfile in wLlist] ))
  percentlist = sorted( percentlist, key=lambda x: float(x) )
  return [ [preamb+'percent_'+f] for f in percentlist]

def MeanThetaFEComputation(SimRun, OutDir):
  #[data, params]
  Simnums = len(SimRun)
  FEDat   = np.zeros((len(SimRun), 6))
  
  for i in range(Simnums):
    p = SimRun[i][1];dat = SimRun[i][0]
    FEDat[i][0] = np.nanmean(dat['Thetaval'])
    FEDat[i][1] = -0.5*p.biasstr*(FEDat[i][0] - p.theta0)/p.A
    FEDat[i][3] = np.nanmean(dat['NBot1'] + dat['NTop1']) /p.nodes
    FEDat[i][4] = np.nanmean(dat['NBot3'] + dat['NTop3']) /p.nodes
    FEDat[i][5] = np.nanmean(dat['NBot2'] + dat['NTop2']) /p.nodes
    
  FEDat[::-1, 2] = integrate.cumtrapz(FEDat[::-1, 1], FEDat[::-1, 0], initial=0)
  Outname = OutDir+'/ThetaAveFE_percent_'+GetValFromStr(p.filename, 'percent')+'_.dat'
  np.savetxt(Outname, FEDat, delimiter=' ', header='z(1) MeanFor(2) FE(3) f_S(4) f_L(5) f_B(6)\n')
  
  

InDir="Harmrawdata/DataFiles/"
OutDir = "intermediatedatafiles"

lstart=1000;
binnums=100
maxpoints=100000

for runtags in ListUniquePercentages(InDir)[:]:
  SimRun = ObtainData(InDir, OutDir, runtags, binnums, lstart, maxpoints)
  MeanThetaFEComputation(SimRun, OutDir)
  #MakeFreeEnergy(SimRun, binnums)
  #MakePartialFreeEnergies('Thetaval', SimRun, np.linspace(10.5, 120., 100) )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

#rhotot     = 8000
#perclist   = [ 1, 10, 25, 50, 75, 100]
#USB        = 10.0
#ULB        = 1.0
#LSBinder   = 4.5
#LNonBinder = 17.5
#LLBinder   = 19.0

#def InitParams():
  #Params             = paramssmall()
  #Params.numptypes   = 3
  #Params.lengths     = [LSBinder, LNonBinder, LLBinder, LSBinder, LNonBinder, LLBinder]
  #Params.bindingstrs = [USB, 0.0, ULB]
  #Params.diams       = [Params.dx for i in range(2*Params.numptypes)]
  #Params.rads        = [d/2       for d in Params.diams]
  #return Params
  
#def MakeDensities(percent):
  #rhoBind    = rhotot*float(percent)/100
  #rhoNonBind = max(0.0, rhotot - rhoBind)
  #return [200.0, rhoNonBind, rhoBind, 200.0, rhoNonBind, rhoBind]
