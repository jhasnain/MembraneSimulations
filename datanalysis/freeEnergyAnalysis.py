import os 
import numpy as np
import pickle
import plotmembranedatafiles as memplt

#Routines to store contents of data file and for data processing
def ExtractValWithKeyword(mylist, keyword):
  try:
    a=[elem for elem in mylist if keyword in elem][0]
    a=a.split(':')[1].split()
    a=[ int(elem) if elem.lstrip('+-').isdigit() else float(elem) for elem in a]
    if len(a)==1 and ( ('Binding str' != keyword) and ('Binding ranges' != keyword) ):
      a=a[0]
    return a
  except:return -1

def ExtractIndexofKeyword(mylist, keyword):return [i for i, elem in enumerate(mylist) if keyword in elem][0]

def ComputeChemPot(densities, alpha):
  zbar = np.sum(densities)
  zbar = zbar/(alpha-zbar)
  muvals = [np.log((1.0 + zbar)*d/alpha) for d in densities]
  return muvals

def GetkVecs(N, dx):
  kvecs = np.meshgrid( 
           np.fft.fftshift(np.fft.fftfreq(N, d=dx)),
           np.fft.fftshift(np.fft.fftfreq(N, d=dx))
          )
  kvecs[0]*=2.0*np.pi
  kvecs[1]*=2.0*np.pi
  return kvecs

def GetHeightVar(kc, A, kvecs):
  fk = kc*(kvecs[0]**2 + kvecs[1]**2)**2
  return np.sum(1.0/fk[fk>0.0])/A

class params:
  def __init__(self, data, tag, binnums, InFile, OutDir):
    self.binnums = binnums
    
    self.Lx  = ExtractValWithKeyword(data, "Lx")
    self.Ly  = ExtractValWithKeyword(data, "Ly")
    self.A   = self.Lx*self.Ly
    self.kc  = ExtractValWithKeyword(data, "Modulus")
    
    self.theta0  = ExtractValWithKeyword(data, "Bias location")
    self.biasstr = ExtractValWithKeyword(data, "Bias strength")
    
    self.nx      = ExtractValWithKeyword(data, "in x")
    self.ny      = ExtractValWithKeyword(data, "in y")
    self.nodes   = self.nx*self.ny
    
    self.numptypes = ExtractValWithKeyword(data, "Ptype num")
    self.alpha     = 1e6*self.nodes/self.A
    
    self.densities = ExtractValWithKeyword(data, "Densities")
    if self.densities==-1:self.muvals = ExtractValWithKeyword(data, "ChemPot vals")
    else:                 self.muvals = np.array(ComputeChemPot(self.densities[:self.numptypes], self.alpha) + ComputeChemPot(self.densities[self.numptypes:], self.alpha))
    
    self.dx = self.Lx/self.nx
    self.dy = self.Ly/self.ny
    
    if self.numptypes>0:
      self.lengths        = ExtractValWithKeyword(data, "Lengths")
      self.bindingstrs    = ExtractValWithKeyword(data, "Binding str")
      self.bindingran     = ExtractValWithKeyword(data, "Binding ranges")
      
  
    self.kvecs     = GetkVecs(self.nx, self.Lx/self.nx)
    self.heightvar = GetHeightVar(self.kc, self.A, self.kvecs)
    
    self.infile   = InFile
    self.OutDir   = OutDir
    self.filename = "{}_Nodes_{}_NumPTypes_{}_kc_{:.2f}_A_{}".format(
                      tag,
                      self.nodes,
                      self.numptypes,
                      self.kc,
                      self.A
                    )
    #self.OutDir=self.OutDir+self.filename
    #if os.path.isdir(self.OutDir)==False:os.mkdir(self.OutDir)

#Read data file and create both parameter object and datafiles
def PrepFileData(filedata, Params, lstart, maxpoints):
  keys = [k.split('(')[0] for k in filedata.pop(0).replace('#','').split()]
  limit = len(filedata[0].split())
  filedata = np.array([ [int(y) if y.lstrip('+-').isdigit() else float(y) for y in x.split() ] for x in filedata[lstart:min(maxpoints+lstart, len(filedata))] if len(x.split())==limit ])
  
  data={};
  for k in range(len(keys)):data[keys[k]]=filedata[:, k]
  data['Thetaval'] = data['ThetaValTop'] - data['ThetaValBot'] + Params.theta0
  data['TotPop']   = data['TotPopTop']   + data['TotPopBot']
  
  return data
  
def FormatData(filedata, d, binnums, OutDir, lstart, maxpoints):
  index = ExtractIndexofKeyword(filedata, "#Sweep")
  Params=params(filedata[:index], d.split('Data_')[1].split('_Nodes')[0], binnums, d, OutDir)
  data = PrepFileData(filedata[index:], Params, lstart, maxpoints)
  return [data, Params] 
 
def ObtainData(InDir, OutDir, runtag, binnums, lstart, maxpoints):
  runtag.append("Data")
  datfiles=filter(lambda x: all( r+'_' in x for r in runtag), os.listdir(InDir)  )
  datfiles=[InDir+i for i in datfiles]
  datfiles = sorted(datfiles, key= lambda x: float(x.split('theta0_')[1].split('_')[0]), reverse=False)
  
  if datfiles==[]:
    print "Error: this tag {} does not find any hits".format(runtag)
    exit(1)
  
  alldata=[]
  for d in datfiles[:]:
    data = open(d, "r").readlines()
    print d, len(data)
    data = FormatData(data, d, binnums, OutDir, lstart, maxpoints)
    alldata.append(data)
  
  print alldata[0][0].keys()
  
  for k in range(len(alldata)):
    theta0 = alldata[k][1].theta0/2
    data   = alldata[k][0]
    data['ThetaValTop'] += theta0
    data['ThetaValBot'] +=-theta0
    #Thetaval is a derived quantity and has already been shifted appropriately
    
  return alldata

#Perform data analysis and prepare to store results 
def Outname(Params, prefix, suffix): return Params.OutDir+'/'+prefix+'_'+Params.filename+'_'+suffix

def Get_Theta_Bins_Edges(SimRun, binnums):
  obsmax=-1000000
  obsmin= 1000000
  for Sim in SimRun:
    obsmax=max( np.nanmax(Sim[0]["Thetaval"]), obsmax )
    obsmin=min( np.nanmin(Sim[0]["Thetaval"]), obsmin )
  thetabins  = np.linspace(obsmin, obsmax, binnums, endpoint=True)
  dtheta     = thetabins[-1] - thetabins[-2]
  thetaedges = thetabins - dtheta/2
  thetaedges = np.append(thetaedges, thetaedges[-1] + dtheta)
  
  return thetabins, thetaedges, dtheta

def MakeFreeEnergies(SimRun, thetabins, thetaedges, dtheta):
  memplt.InitPlot(r'step', r'$\theta$', SimRun[0][1])
  for sim in SimRun[::-1]:memplt.PlotPlot(np.arange(sim[0]["Thetaval"].shape[0]), sim[0]["Thetaval"], '', '-', None, 2.2, None)
  memplt.PlotSave(SimRun[0][1], "thetatraj" )
  
  prev   = []
  Full   = []
  memplt.InitPlot(r'membrane gap [nm]', r'$\beta \Delta F$', SimRun[0][1])
  for sim in SimRun[::-1]:
    Params = sim[1]
    Dat    = sim[0]
    
    ave = np.mean(Dat["BiasE"])
    hist, edges=np.histogram(Dat["Thetaval"], bins=thetaedges, weights=np.exp(Dat["BiasE"]-ave), density=True )
    partial = np.ones(hist.shape)*np.nan;mask=hist>0.0
    partial[mask] = -np.log(hist[mask]) + ave
    
    if prev==[]:av=np.nanmean(partial[np.isnan(partial)==False])
    else:
      av = partial-prev;nanmean=np.nanmean(av)
      if np.isnan(nanmean)==True:
        av=partial[np.where(np.isnan(partial)==False)[0][-1]] - prev[np.where(np.isnan(prev)==False)[0][0]]
      else:av = nanmean;
    
    partial-=av
    memplt.PlotPlot(thetabins, partial, '', '-', None, 2.2, r'$\theta_0={}$'.format(Params.theta0))
    prev = partial
    Full.append(partial)
    
  memplt.PlotLengths(SimRun[0][1])
  memplt.PlotSave(SimRun[0][1], 'PartialFE_thetaval', lg=False)
  
  Full=np.array(Full)
  Final=np.ones(thetabins.shape)
  for i in range(Full.shape[1]):Final[i]=np.nanmean(Full[:,i])
  return Final
 
def MakeExpectations(SimRun, thetabins, thetaedges, dtheta):
  skipkeys = ["TAcc", "NAcc", "SHAcc", "ThetaValTop", "ThetaValBot", "Sweep", "Thetaval"]
  keys = filter(lambda x: x not in skipkeys  , SimRun[0][0].keys())
  flatdat = [{} for i in range(thetabins.shape[0])]
  for i in range(len(flatdat)):
    for k in keys:flatdat[i][k] = np.array([])
  
  for sim in SimRun[::-1]:
    Params = sim[1]
    Dat    = sim[0]
    thetavals = Dat["Thetaval"]
    
    for t, edge in np.ndenumerate(thetaedges[:-1]):
      mask = (thetavals >= edge ) & (thetavals < edge + dtheta)
      t_ind = t[0]
      if np.any(mask)==True: 
        for k in keys: flatdat[t_ind][k] = np.append( flatdat[t_ind][k], Dat[k][mask] )
  
  Averages = np.nan*np.ones( (thetabins.shape[0], len(keys)) )
  for i in range(thetabins.shape[0]):
    bvals = flatdat[i]["BiasE"]
    if bvals.shape[0]>0:
      for ind, k in enumerate(keys):
        Averages[i, ind] = np.average(flatdat[i][k], weights = bvals)
  
  ExpectationVals = {}
  for ind, k in enumerate(keys):ExpectationVals[k] = Averages[:, ind]
  
  return ExpectationVals

def WriteParams(Params):
  exceptions = ["kvecs", "OutDir", "filename", "infile", "theta0"]
  footer=''
  for prop, value in vars(Params).iteritems():
    if prop not in exceptions:
      footer+=str(prop)+': {}'.format(value)+'\n'
  return footer
  
def StoreFEData(FEData, Params):
  fh = open(Outname(Params, 'FEData', '.obj'), 'wb')
  pickle.dump([FEData, Params], fh)
  fh.close()
  
  keys = FEData.keys()
  keys.insert(0, keys.pop(keys.index("thetavals")))
  keys.insert(1, keys.pop(keys.index("FE")))
  
  nparray_FEData = np.nan*np.ones( (FEData["thetavals"].shape[0], len(keys)) ) 
  for ind, k in enumerate(keys):nparray_FEData[:, ind] = FEData[k][:]
  
  np.savetxt( Outname(Params, 'FullDat', '.dat'), nparray_FEData, 
              header = ' '.join(keys), 
              footer = WriteParams(Params)
            )
   
def PerformFreeEnergyCalculations(SimRun, binnums):
  thetabins, thetaedges, dtheta = Get_Theta_Bins_Edges(SimRun, binnums)
  FECurve                       = MakeFreeEnergies(SimRun, thetabins, thetaedges, dtheta) 
  ExpectationVals               = MakeExpectations(SimRun, thetabins, thetaedges, dtheta) 
  FEData                        = ExpectationVals;
  FEData["thetavals"]           = thetabins;
  FEData["FE"]                  = FECurve
  StoreFEData(FEData, SimRun[0][1])

    
