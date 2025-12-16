import os 
import numpy as np
import matplotlib.pyplot as plt
#import pymbar
import pickle


Codes={
       'TotPop'  : r'$\langle N_{{\mathrm{{Tot}}}}\rangle$',
       'MeshE'   : r'$\langle E_{{\mathrm{{Mem}}}}\rangle$', 
       'TotArea' : r'$\langle A_{{\mathrm{{Mem}}}}\rangle/A_0$', 
       'PartE'   : r'$\langle \beta U_{{B}}\rangle$',
       'etaAve'  : r'$\langle \eta \rangle$',
       'rhoAve'  : r'$\langle \rho \rangle [1/ \mu m^{{2}}]$',
      }

maxpoints=50000
lstart=1000;

def InitPlot(xlabel, ylabel, Params):
  plt.clf()
  plt.rc('text', usetex=True);plt.rc('font', family='serif', size=18);
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(r'PTypes$={}$ $\beta\kappa_B={} A_0 ={}\,\mu m^2$'.format(Params.numptypes, Params.kc, Params.A/1e6) )
  #return plt

def PlotPlot(x, y, mark, line, msize, lw, label):plt.plot(x, y, marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label)
  
def PlotLengths(Params):
  for i in range(Params.numptypes):
    plt.axvline(x=2.0*(Params.lengths[i] + Params.rads[i]), 
                color='k', linestyle='--', 
                label=r'$2 L_{}={}$'.format(i+1, 2.0*(Params.lengths[i] + Params.rads[i]))
               )
  plt.axvline(x=2.0*np.sqrt(Params.heightvar), color='k', linestyle='-', label=r'$2\sigma_h = {:.1f}$'.format(2.0*np.sqrt(Params.heightvar)))

def PlotSave(Params, name):
    plt.legend()
    plt.tight_layout()
    plt.savefig(Outname(Params, name, '.png'))
   
def GetParamFilename(Dir):
  return Dir+'/'+filter(lambda x: 'Params_'in x and '.obj' in x, os.listdir(Dir))[0]

def Outname(Params, prefix, suffix):return Params.OutDir+'/'+prefix+'_'+Params.filename+'_'+suffix

def StoreParams(Params):
  filehandler = open(Outname(Params, 'Params', '.obj'), 'wb')
  pickle.dump(Params, filehandler)
  filehandler.close()
  
def LoadParams(Dir):
  filehandler = open(GetParamFilename(Dir), 'rb')
  Params=pickle.load(filehandler)
  return Params
  
def GetkVecs(N, dx):
  kvecs = np.meshgrid( 
           np.fft.fftfreq(N, d=dx),
           np.fft.fftfreq(N, d=dx)
          )
  kvecs[0]*=2.0*np.pi
  kvecs[1]*=2.0*np.pi
  return kvecs

def GetHeightVar(kc, A, kvecs):
  fk = kc*(kvecs[0]**2 + kvecs[1]**2)**2
  print np.sum(1.0/fk[fk>0.0])/A
  exit(1)
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
    self.muvals    = ExtractValWithKeyword(data, "ChemPot")
    self.Diameters = ExtractValWithKeyword(data, "Diameters")
    
    if self.numptypes>0:
      self.mubot          = self.muvals[:self.numptypes]
      self.mutop          = self.muvals[self.numptypes:2*self.numptypes]
      
      if self.numptypes==1:
        self.rads           = [self.Diameters/2]
        self.lengths        = [ExtractValWithKeyword(data, "Lengths")]
        self.packingfactors = [np.pi*self.rads[0]**2]
      else:
        self.rads           = [d/2 for d in self.Diameters]
        self.lengths        = ExtractValWithKeyword(data, "Lengths")
        self.packingfactors = [np.pi*r**2 for r in self.rads]
    
    self.kvecs     = GetkVecs(self.nx, self.Lx/self.nx)
    self.heightvar = GetHeightVar(self.kc, self.A, self.kvecs)
    
    self.infile   = InFile
    self.OutDir   = OutDir
    self.filename = "{}_Nodes_{}_NumPTypes_{}_kc_{}_A_{}".format(
                      tag,
                      self.nodes,
                      self.numptypes,
                      self.kc,
                      self.A
                    )
    self.OutDir=self.OutDir+self.filename
    if os.path.isdir(self.OutDir)==False:os.mkdir(self.OutDir)
  
def Get_Theta_Bins_Edges(SimRun):
  thetavals=np.array([]);binnums=SimRun[0][1].binnums
  for Sim in SimRun:thetavals=np.hstack((thetavals, Sim[0]['Thetaval']))
  thetabins  = np.linspace(np.min(thetavals), np.max(thetavals), num=binnums, endpoint=True)
  dtheta = (thetabins[-1]-thetabins[-2])
  thetaedges = thetabins-dtheta/2;thetaedges=np.append(thetaedges, thetaedges[-1]+dtheta)
  return thetavals, thetabins, thetaedges, dtheta

def ExtractValWithKeyword(mylist, keyword):
  try:
    a=[elem for elem in mylist if keyword in elem][0]
    a=a.split(':')[1].split()
    a=[ int(elem) if elem.lstrip('+-').isdigit() else float(elem) for elem in a]
    if len(a)==1:a=a[0]
    else: a=[float(b) for b in a]  #This is a hack and you should remove it for next gen simulation data, oh and dont pass command line floats as ints in your submit scripts plx
    return a
  except:return -1

def ExtractIndexofKeyword(mylist, keyword):return [i for i, elem in enumerate(mylist) if keyword in elem][0]

def PrepFileData(filedata, Params):
  keys = [k.split('(')[0] for k in filedata.pop(0).replace('#','').split()]
  limit = len(filedata[0].split())
  
  filedata = np.array([ [int(y) if y.lstrip('+-').isdigit() else float(y) for y in x.split() ] for x in filedata[lstart:min(maxpoints+lstart, len(filedata))] if len(x.split())==limit ])
  data={};
  
  for k in range(len(keys)):data[keys[k]]=filedata[:, k]
    
  data['Thetaval'] = data['ThetaValTop']-data['ThetaValBot'] + Params.theta0
  
  if Params.numptypes>0:
    for p in range(1, Params.numptypes+1):
      ind=str(p)
      data['etaBot'+ind] = Params.packingfactors[p-1]*data['NBot'+ind]/data['TotAreaBot']
      data['etaTop'+ind] = Params.packingfactors[p-1]*data['NTop'+ind]/data['TotAreaBot']
      
    data['TotetaBot'] = data['etaBot1']
    data['TotetaTop'] = data['etaTop1']  
    for p in range(2, Params.numptypes+1):
      ind=str(p)
      data['TotetaBot']+=data['etaBot'+ind]
      data['TotetaTop']+=data['etaTop'+ind]  
  
  return data
  
def FormatData(filedata, d, binnums, OutDir):
  index = ExtractIndexofKeyword(filedata, "#Sweep")
  Params=params(filedata[:index], d.split('Data_')[1].split('_Nodes')[0], binnums, d, OutDir)
  data = PrepFileData(filedata[index:], Params)
  return [data, Params]

def ObtainData(InDir, OutDir, runtag, binnums):
  datfiles=filter(lambda x: all( r in x for r in runtag), os.listdir(InDir)  )
  datfiles=[InDir+i for i in datfiles]
  datfiles = sorted(datfiles, key= lambda x: float(x.split('theta0_')[1].split('_')[0]), reverse=False)
  alldata=[]
  for d in datfiles[:]:
    print d
    data = open(d, "r").readlines()
    data = FormatData(data, d, binnums, OutDir)
    alldata.append(data)
  print alldata[0][0].keys()
  return alldata

def ShiftTheta(SimRun):
  for k in range(len(SimRun)):
    theta0 = SimRun[k][1].theta0/2
    data   = SimRun[k][0]
    data['ThetaValTop'] += theta0
    data['ThetaValBot'] +=-theta0
    #Since Thetaval is a derived quantity, it has already been shifted appropriately earlier
  return SimRun

def QuantitiesOfInterest(SimRun):
  TheData={}
  for key in Codes.keys():
    if os.path.isfile(Outname(SimRun[0][1], 'Expectation_'+key, '.dat'))==False:
      compileddata=np.array([]);
      
      if key=='MeshE':
        for Sim in SimRun:compileddata=np.hstack((compileddata, Sim[0]["MeshE"]))
      if key=='TotArea':
        for Sim in SimRun:compileddata=np.hstack( (compileddata, (Sim[0]["TotAreaBot"] + Sim[0]["TotAreaTop"])/(2.0*Sim[1].A)) )
      
      if SimRun[0][1].numptypes>0:
        if key=='TotPop':
          for Sim in SimRun:compileddata=np.hstack(( compileddata, Sim[0]["TotPop"] ))
        if key=='etaAve':
          for Sim in SimRun:compileddata=np.hstack(( compileddata, 0.5*(Sim[0]["TotetaTop"] + Sim[0]["TotetaBot"]) ))
        if key=='rhoAve':
          for Sim in SimRun:compileddata=np.hstack(( compileddata, 1e6*0.5*( Sim[0]["TotBotPop"]/Sim[0]["TotAreaBot"] + Sim[0]["TotTopPop"]/Sim[0]["TotAreaTop"]) ))
        if key=='PartE':
          for Sim in SimRun:compileddata=np.hstack(( compileddata, (Sim[0]["PartE"]) ))
      if compileddata.shape[0]>0:TheData[key]=[Codes[key], compileddata]
      
  return TheData

def CompileExpectations(SimRun):
  xlabel=r'membrane gap [nm]'
  thetavals, thetabins, thetaedges, dtheta = Get_Theta_Bins_Edges(SimRun)
  for name, [ylabel, fulldata] in QuantitiesOfInterest(SimRun).iteritems():
    hist=np.ones(thetabins.shape)*np.nan
    for i, binval in enumerate(thetaedges):
      mask = (thetavals >= binval ) & (thetavals < binval + dtheta)
      if np.any(mask)!=False:hist[i] = np.mean(fulldata[mask])
    np.savetxt( Outname(SimRun[0][1], 'Expectation_'+name, '.dat'), np.column_stack((thetabins, hist)) )
    
    InitPlot(xlabel, ylabel, SimRun[0][1])
    delta=max(1, int(fulldata.shape[0]/1000))
    PlotPlot(thetavals[::delta], fulldata[::delta], '', '-', None, 1, None)
    PlotPlot(thetabins, hist, 'x', '-', 2, 2, None)
    PlotLengths(SimRun[0][1])
    PlotSave(SimRun[0][1], 'FullDat_Expectation_'+name)

def DistDat(SimRun):
  BondLengths=[];meangap=[];vargap=[];
  for Sim in SimRun:
    BondLengths.append( 2.0*(Sim[1].lengths[0] + Sim[1].rads[0]) - (Sim[1].lengths[1] + 2*Sim[1].rads[1]) )
    tmp = BondLengths[-1]*np.ones(Sim[0]["Thetaval"].shape)
    PlotPlot(tmp, Sim[0]["Thetaval"], 'x', '', None, None, None)
    meangap.append(Sim[0]["Thetaval"].mean())
    vargap.append(Sim[0]["Thetaval"].var())
  np.savetxt( Outname(SimRun[0][1], 'DistDat', '.dat'), np.column_stack((BondLengths, meangap, vargap)) )

def PopulationDat(SimRun):
  BondLengths=[];meanp1=[];varp1=[];meanp2=[];varp2=[];
  for Sim in SimRun:
    p1dat = Sim[0]['NBot1'] + Sim[0]['NTop1']
    p2dat = Sim[0]['NBot2'] + Sim[0]['NTop2']
    BondLengths.append( 2.0*(Sim[1].lengths[0] + Sim[1].rads[0]) - (Sim[1].lengths[1] + 2*Sim[1].rads[1]) )
    meanp1.append(np.mean(p1dat))
    varp1.append(np.var(p1dat))
    meanp2.append(np.mean(p2dat))
    varp2.append(np.var(p2dat))
    #print BondLengths[-1], np.mean(p1dat), np.mean(p2dat), 2*Sim[1].rads[1], Sim[1].rads[0]
  np.savetxt( Outname(SimRun[0][1], 'PopDat', '.dat'), np.column_stack((BondLengths, meanp1, meanp2, varp1, varp2)) )

def PlotExpectations(OutDir):
  for dirpath, dirnames, filenames in os.walk(OutDir):
    for f in [h for h in filenames if ("_.dat" in h and "Expectation_" in h)]:
      name=dirpath+'/'+f;
      dat = np.loadtxt(name)
      Params = LoadParams(dirpath)
      
      keyword=f.split("Expectation_")[1].split('_')[0]
      InitPlot(r'membrane gap [nm]', Codes[keyword], Params)
      PlotPlot(dat[:,0], dat[:,1], 'x', '-', 7, 2.2, None)
      PlotLengths(Params)
      PlotSave(Params, 'Expectation'+keyword)

def PlotDistDat(OutDir):
  for dirpath, dirnames, filenames in os.walk(OutDir):
    for f in [h for h in filenames if ("_.dat" in h and "DistDat_" in h)]:
      name=dirpath+'/'+f;
      dat = np.loadtxt(name)
      Params = LoadParams(dirpath)
      InitPlot(r'Protein Size [nm]', 'Gap [nm]', Params)
      PlotPlot(dat[:,0], dat[:,1], 'x', '-', 7, 2.2, None)
      PlotSave(Params, 'DistDat')

def PlotPopDat(OutDir):
  for dirpath, dirnames, filenames in os.walk(OutDir):
    for f in [h for h in filenames if ("_.dat" in h and "PopDat_" in h)]:
      name=dirpath+'/'+f;
      dat = np.loadtxt(name)
      Params = LoadParams(dirpath)
      InitPlot(r'Protein Size [nm]', 'nums', Params)
      PlotPlot(dat[:,0], dat[:,1], 'x', '-', 7, 2.2, None)
      PlotPlot(dat[:,0], dat[:,2], 'o', '-', 7, 2.2, None)
      PlotSave(Params, 'PopDat')

def MyPlots(OutDir):
  #PlotExpectations(OutDir)
  PlotDistDat(OutDir)
  PlotPopDat(OutDir)
  
binnums=100
tags=[ ["rod"] ]

InDir="/Data/FletcherMembrane/DataFiles/"
OutDir = "../plots/exclusion/"

for tag in tags:
  SimRun     = ObtainData(InDir, OutDir, tag, binnums)
  SimRun     = ShiftTheta(SimRun)
  DistDat(SimRun)
  PopulationDat(SimRun)
  #CompileExpectations(SimRun)
  StoreParams(SimRun[0][1])
  
MyPlots(OutDir)

