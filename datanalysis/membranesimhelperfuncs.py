import os 
import numpy as np
import matplotlib.pyplot as plt
import pickle
from   effectivepotentials import *

Codes={
       'MeshE'   : r'$\langle E_{{\mathrm{{Mem}}}}\rangle$', 
       'TotArea' : r'$\langle A_{{\mathrm{{Mem}}}}\rangle/A_0$', 
       'PartE'   : r'$\langle \beta U_{{B}}\rangle$',
       'TotPop'  : r'$\langle N_{{\mathrm{{Tot}}}}\rangle$',
       'etaAve'  : r'$\langle \eta \rangle$',
       'rhoAve'  : r'$\langle \rho \rangle [1/ \mu m^{{2}}]$',
       'Pops' : r'$\langle N_{i}\rangle$',
       'eta'  : r'$\langle \eta_i \rangle$',
       'rho'  : r'$\langle \rho_i \rangle [1/ \mu m^{{2}}]$',
      }

def ExtractValWithKeyword(mylist, keyword):
  try:
    a=[elem for elem in mylist if keyword in elem][0]
    a=a.split(':')[1].split()
    a=[ int(elem) if elem.lstrip('+-').isdigit() else float(elem) for elem in a]
    if len(a)==1 and ( ('Binding str' != keyword) and ('Binding ranges' != keyword) ):
      a=a[0]
    #else: a=[float(b) for b in a]  #This is a hack and you should remove it for next gen simulation data, oh and dont pass command line floats as ints in your submit scripts plx
    return a
  except:return -1

def ExtractIndexofKeyword(mylist, keyword): return [i for i, elem in enumerate(mylist) if keyword in elem][0]

def InitPlot(xlabel, ylabel, Params):
  plt.clf()
#  plt.rc('text', usetex=True);plt.rc('font', family='serif', size=18);
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if Params!=None:plt.title(r'PTypes$={}$ $\beta\kappa_B={} A_0 ={}\,\mu m^2$'.format(Params.numptypes, Params.kc, Params.A/1e6) )
  
def PlotPlot(x, y, mark, line, msize, lw, label, c, error = None):
  mask=np.where(np.isnan(y)==False)
  option =          (2 if error is not None else 0)
  option = option + (1 if     c is not None else 0)
  
  if option == 0: g = plt.plot(x[mask], y[mask], marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label)
  if option == 1: g = plt.plot(x[mask], y[mask], marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label, color=c)
  if option == 2: g = plt.errorbar(x[mask], y[mask], yerr = error[mask], marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label)
  if option == 3: g = plt.errorbar(x[mask], y[mask], yerr = error[mask], marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label, color=c)
  
  #if c==None:g= plt.plot(x[mask], y[mask], marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label)
  #else:      g= plt.plot(x[mask], y[mask], marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label, color=c)
  
  return g

def PlotLengths(Params):
  for i in range(Params.numptypes):
    plt.axvline(x=(Params.lengths[i] + Params.rads[i]), 
                color='k', linestyle=':', 
                label=r'$L_{}={}$'.format(i+1, (Params.lengths[i] + Params.rads[i]))
               )
    #plt.axvline(x=2.0*(Params.lengths[i] + Params.rads[i]), 
                #color='k', linestyle='--', 
                #label=r'$2 L_{}$'.format(i+1)
               #)
    top = (Params.lengths[i] + Params.rads[i]) + np.sqrt(Params.heightvar)
    bot = (Params.lengths[i] + Params.rads[i]) - np.sqrt(Params.heightvar)
    plt.axvline(x= top, color='green', lw=0, label=r'$2\sigma_h = {:.1f}$'.format(2.0*np.sqrt(Params.heightvar)))
    #plt.axvline(x=bot, color='green', alpha=0.5)
    ymin, ymax = plt.ylim()
    plt.fill_between(np.linspace(bot, top, 100), ymin, ymax, alpha=0.3, facecolor='green')
    plt.ylim((ymin, ymax))
    
def PlotSave(Params, name, cols, loc):
  #plt.legend()
  if loc!=None: plt.legend(loc='lower left', prop={'size': 12}, ncol=cols, bbox_to_anchor=loc )
  else:         plt.legend(loc='best', prop={'size': 12}, ncol=cols )
  plt.tight_layout()
  if "Summary" not in name:plt.savefig(Outname(Params, name, '.png'))
  else:plt.savefig(Params.OutDir+"/"+name+"_.png")
    
def GetParamFilename(Dir):return Dir+'/'+filter(lambda x: 'Params_'in x and '.obj' in x, os.listdir(Dir))[0]

def Outname(Params, prefix, suffix):return Params.OutDir+'/'+prefix+'_'+Params.filename+'_'+suffix

def StoreParams(Params):
  filehandler = open(Outname(Params, 'Params', '.obj'), 'wb')
  pickle.dump(Params, filehandler)
  filehandler.close()
  
def LoadParams(name):
  filehandler = open(name, 'rb')
  Params=pickle.load(filehandler)
  return Params

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
    
    self.numptypes  = ExtractValWithKeyword(data, "Ptype num")
    self.densities  = ExtractValWithKeyword(data, "Densities")
    
    self.dx = self.Lx/self.nx
    self.dy = self.Ly/self.ny
    
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
    
    if self.numptypes>0:
      self.Diameters      = [4.0 for i in range(2*self.numptypes)]
      self.rads           = [d/2 for d in self.Diameters]
      self.lengths        = ExtractValWithKeyword(data, "Lengths")
      self.packingfactors = [np.pi*r**2 for r in self.rads]
      self.bindingstrs    = ExtractValWithKeyword(data, "Binding str")
      self.bindingran     = ExtractValWithKeyword(data, "Binding ranges")
      self.fugacities = self.convertdensitiestoFugacity(self.densities)
    
    
  def convertdensitiestoFugacity(self, densities):
    alpha  = 1e6/(self.Diameters[0]**2) 
    zbar   = np.sum(densities[:self.numptypes]);
    zbar   = zbar/(alpha - zbar)
    muvals = [ d*(1.0 + zbar)/alpha for d in self.densities[:self.numptypes] ]
    
    zbar   = np.sum(densities[self.numptypes:]);
    zbar   = zbar/(alpha - zbar)
    muvals +=[ d*(1.0 + zbar)/alpha for d in self.densities[self.numptypes:] ]
          
    return muvals

def PrepFileData(filedata, Params, lstart, maxpoints, d):
  keys = [k.split('(')[0] for k in filedata.pop(0).replace('#','').split()]
  limit = len(filedata[0].split())

  try:
    filedata = np.array([ [int(y) if y.lstrip('+-').isdigit() else float(y) for y in x.split() ] for x in filedata[lstart:min(maxpoints+lstart, len(filedata))] if len(x.split())==limit ])
  
  except:
    for i, x in enumerate(filedata[lstart:len(filedata)]):
      y = x.split()
      print i, y, d
      a = [int(y) if y.lstrip('+-').isdigit() else float(y) for y in x.split() ]
      

  data={};
  for k in range(len(keys)):data[keys[k]]=filedata[:, k]
  data['Thetaval'] = data['ThetaValTop']-data['ThetaValBot'] + Params.theta0
  data['TotPop']   = data['TotPopTop'] + data['TotPopBot']
  
  
  #if Params.numptypes>0:
    #for p in range(1, Params.numptypes+1):
      #ind=str(p)
      #data['etaBot'+ind] = Params.packingfactors[p-1]*data['NBot'+ind]/data['TotAreaBot']
      #data['etaTop'+ind] = Params.packingfactors[p-1]*data['NTop'+ind]/data['TotAreaTop']
      
    #data['TotetaBot'] = np.zeros(data['etaBot1'].shape)
    #data['TotetaTop'] = np.zeros(data['etaTop1'].shape)
    
    #for p in range(1, Params.numptypes+1):
      #ind=str(p)
      #data['TotetaBot']+=data['etaBot'+ind]
      #data['TotetaTop']+=data['etaTop'+ind]  
  return data
  
def FormatData(filedata, d, binnums, OutDir, lstart, maxpoints):
  try: 
    index  = ExtractIndexofKeyword(filedata, "#Sweep")
    Params = params(filedata[:index], d.split('Data_')[1].split('_Nodes')[0], binnums, d, OutDir)
    data   = PrepFileData(filedata[index:], Params, lstart, maxpoints, d)
    return [data, Params] 
  except: 
    fp = open("error_tmp", "a")
    fp.write(d+'\n')
    fp.close()
    return None
  
 
def ObtainData(InDir, OutDir, runtag, binnums, lstart, maxpoints):
  runtag.append("Data")
  x        = os.listdir(InDir)
  
  datfiles = filter(lambda x: all( r+'_' in x for r in runtag), os.listdir(InDir)  )
  datfiles = [InDir+i for i in datfiles]
  datfiles = sorted(datfiles, key= lambda x: float(x.split('theta0_')[1].split('_')[0]), reverse=False)
  alldata  = []
  
  if datfiles==[]:
    print "Error: this tag {} does not find any hits".format(runtag)
    exit(1)
  
  for d in datfiles[:]:
    data = open(d, "r").readlines()
    print d, len(data)
    data = FormatData(data, d, binnums, OutDir, lstart, maxpoints)
    if data is not None: alldata.append(data)
  
  print alldata[0][0].keys()
  
  for k in range(len(alldata)):
    theta0 = alldata[k][1].theta0/2
    data   = alldata[k][0]
    #Since Thetaval is a derived quantity, it has already been shifted appropriately earlier
    
  StoreParams(alldata[0][1])
  return alldata

def MakePartialFreeEnergies(keyword, SimRun, bin_centers):
  xlabel=r'membrane gap [nm]';ylabel=r'$\beta \Delta F$'
  dtheta = (bin_centers[-1] - bin_centers[-2])
  bins   = bin_centers - dtheta/2
  bins   = np.concatenate( (bin_centers, np.array([bin_centers[-1]+dtheta])) )
  prev   = []
  Full   = []
  
  InitPlot(r'step', r'$\theta$', SimRun[0][1])
  for sim in SimRun[::-1]:
    PlotPlot(np.arange(sim[0][keyword].shape[0]), sim[0][keyword], '', '-', None, 2.2, None)
  PlotSave(SimRun[0][1], 'thetatraj_'+keyword)
  
  InitPlot(xlabel, ylabel, SimRun[0][1])
  for sim in SimRun[::-1]:
    Params = sim[1]
    Dat    = sim[0]
    
    ave = np.mean(Dat["BiasE"])
    hist, edges=np.histogram(Dat[keyword], bins=bins, weights=np.exp(Dat["BiasE"]-ave), density=True )
    partial = np.ones(hist.shape)*np.nan;mask=hist>0.0
    partial[mask] = -np.log(hist[mask]) + ave
    
    if prev==[]:av=np.nanmean(partial[np.isnan(partial)==False])
    else:
      av = partial-prev;nanmean=np.nanmean(av)
      if np.isnan(nanmean)==True:
        av=partial[np.where(np.isnan(partial)==False)[0][-1]] - prev[np.where(np.isnan(prev)==False)[0][0]]
      else:av = nanmean;
    
    partial-=av
    PlotPlot(bin_centers, partial, '', '-', None, 2.2, r'$\theta_0={}$'.format(Params.theta0))
    prev = partial
    Full.append(partial)
  
  PlotLengths(SimRun[0][1])
  PlotSave(SimRun[0][1], 'PartialFE_'+keyword)
  
  Full=np.array(Full)
  Final=np.ones(bin_centers.shape)
  for i in range(Full.shape[1]):Final[i]=np.nanmean(Full[:,i])
  
  return bin_centers, Final
  
def MakeFreeEnergy(SimRun, binnums):
  keywords=["Thetaval"]
  obsmax=-1000000
  obsmin= 1000000
  for keyword in keywords:
    for Sim in SimRun:
      obsmax=max( np.nanmax(Sim[0][keyword]), obsmax )
      obsmin=min( np.nanmin(Sim[0][keyword]), obsmin )
    
    thetabinvalues = np.array(np.linspace(obsmin, obsmax, binnums))
    
    FullFE=MakePartialFreeEnergies(keyword, SimRun, thetabinvalues) 
    np.savetxt( Outname(SimRun[0][1], 'FE_'+keyword, '.dat'), np.column_stack((FullFE[0], FullFE[1])) )
    
    InitPlot(r'membrane gap [nm]', r'$\beta \Delta F$', SimRun[0][1])
    PlotPlot(FullFE[0], FullFE[1], '', '-', None, 2.2, None)
    
    #veff = nomemflux.GetFe(SimRun[0][1], np.linspace(0.0,15.0, 100))
    
    #PlotPlot(FullFE[0], veff/(np.amax(veff)), '', '--', None, 2.2, None)
    PlotLengths(SimRun[0][1])
    
    PlotSave(SimRun[0][1], 'FE_'+keyword)
    
def ReturnFreeEnergy(SimRun, binnums):
  keywords=["Thetaval"]
  obsmax=-1000000
  obsmin= 1000000
  for keyword in keywords:
    for Sim in SimRun:
      obsmax=max( np.nanmax(Sim[0][keyword]), obsmax )
      obsmin=min( np.nanmin(Sim[0][keyword]), obsmin )
    
    FullFE=MakePartialFreeEnergies(keyword, SimRun, np.array(np.linspace(obsmin, obsmax, binnums))) 
    np.savetxt( Outname(SimRun[0][1], 'FE_'+keyword, '.dat'), np.column_stack((FullFE[0], FullFE[1])) )
    
    #InitPlot(r'membrane gap [nm]', r'$\beta \Delta F$', SimRun[0][1])
    #PlotPlot(FullFE[0], FullFE[1], '', '-', None, 2.2, None)
    #veff = nomemflux.GetFe(SimRun[0][1], np.linspace(0.0,15.0, 100))
    #PlotPlot(FullFE[0], veff/(np.amax(veff)), '', '--', None, 2.2, None)
    #PlotLengths(SimRun[0][1])
    
    #PlotSave(SimRun[0][1], 'FE_'+keyword)
    
    return FullFE
    
def Get_Theta_Bins_Edges(SimRun):
  thetavals=np.array([]);binnums=SimRun[0][1].binnums
  for Sim in SimRun:thetavals=np.hstack((thetavals, Sim[0]['Thetaval']))
  thetabins  = np.linspace(np.min(thetavals), np.max(thetavals), num=binnums, endpoint=True)
  dtheta = (thetabins[-1]-thetabins[-2])
  thetaedges = thetabins-dtheta/2;thetaedges=np.append(thetaedges, thetaedges[-1]+dtheta)
  return thetavals, thetabins, thetaedges, dtheta

def QuantitiesOfInterest(SimRun):
  TheData={}
  for key in Codes.keys():
    #if os.path.isfile(Outname(SimRun[0][1], 'Expectation_'+key, '.dat'))==False:
    compileddata=np.array([]);
    if key=='MeshE':
      for Sim in SimRun:compileddata=np.hstack((compileddata, Sim[0]["MeshE"]))
    if key=='TotArea':
      for Sim in SimRun:compileddata=np.hstack(( compileddata, (Sim[0]["TotAreaBot"] + Sim[0]["TotAreaTop"])/(2.0*Sim[1].A) ))
    
    if SimRun[0][1].numptypes>0:
      if key=='TotPop':
        for Sim in SimRun:compileddata=np.hstack(( compileddata, Sim[0]["TotPop"] ))
      if key=='etaAve':
        for Sim in SimRun:compileddata=np.hstack(( compileddata, 0.5*(Sim[0]["TotetaTop"] + Sim[0]["TotetaBot"]) ))
      if key=='rhoAve':
        for Sim in SimRun:compileddata=np.hstack(( compileddata, 1e6*0.5*( Sim[0]["TotPopBot"]/Sim[0]["TotAreaBot"] + Sim[0]["TotPopTop"]/Sim[0]["TotAreaTop"]) ))
      
      if SimRun[0][1].numptypes>1:
        popnums=SimRun[0][1].numptypes;
        
        if key=='Pops':
          simlen = SimRun[0][0]["TotPop"].shape[0]
          compileddata = np.ones(( simlen, popnums+1 ))
          compileddata[:, 0] = 0.5*SimRun[0][0]["TotPop"][:];
          for i in range(1, popnums+1):compileddata[:, i]=0.5*(SimRun[0][0]["NBot"+str(i)] + SimRun[0][0]["NTop"+str(i)])[:]
                      
          for Sim in SimRun[1:]:
            simlen = Sim[0]["TotPop"].shape[0]
            tmp = np.ones(( simlen, popnums+1 ))
            tmp[:, 0] = 0.5*Sim[0]["TotPop"][:];
            for i in range(1, popnums+1):tmp[:, i]=0.5*(Sim[0]["NBot"+str(i)] + Sim[0]["NTop"+str(i)])[:]
            compileddata=np.vstack(( compileddata, tmp ))
        
        if key=='eta':
          simlen = SimRun[0][0]["TotetaTop"].shape[0]
          compileddata = np.ones(( simlen, popnums+1 ))
          compileddata[:, 0] =  0.5*(SimRun[0][0]["TotetaTop"] + SimRun[0][0]["TotetaBot"])[:];
          for i in range(1, popnums+1):
            compileddata[:, i]=0.5*(SimRun[0][0]["etaTop"+str(i)] + SimRun[0][0]["etaBot"+str(i)])[:]
                      
          for Sim in SimRun[1:]:
            simlen = Sim[0]["TotetaTop"].shape[0]
            tmp = np.ones(( simlen, popnums+1 ))
            tmp[:, 0] = 0.5*(Sim[0]["TotetaTop"] + Sim[0]["TotetaBot"])[:];
            for i in range(1, popnums+1):tmp[:, i]=0.5*(Sim[0]["etaTop"+str(i)] + Sim[0]["etaBot"+str(i)])[:]
            compileddata=np.vstack(( compileddata, tmp ))
        #if key=='eta':
          #for Sim in SimRun:compileddata=np.hstack(( compileddata, 0.5*(Sim[0]["TotetaTop"] + Sim[0]["TotetaBot"]) ))
        #if key=='rho':
          #for Sim in SimRun:compileddata=np.hstack(( compileddata, 1e6*0.5*( Sim[0]["TotBotPop"]/Sim[0]["TotAreaBot"] + Sim[0]["TotTopPop"]/Sim[0]["TotAreaTop"]) ))   
      if key=='Veff':
        for Sim in SimRun:compileddata=np.hstack(( compileddata, (Sim[0]["Veff"]) ))
    if compileddata.shape[0]>0:TheData[key]=[Codes[key], compileddata]
  return TheData
    
def ParticleProps(SimRun):
  xlabel=r'membrane gap [nm]'
  thetavals, thetabins, thetaedges, dtheta = Get_Theta_Bins_Edges(SimRun)
  for name, [ylabel, fulldata] in QuantitiesOfInterest(SimRun).iteritems():
    if len(fulldata.shape)>1:hist=np.ones(( thetabins.shape[0], fulldata.shape[1] ))*np.nan
    else:hist=np.ones(thetabins.shape)*np.nan
    
    for i, binval in enumerate(thetaedges):
      mask = (thetavals >= binval ) & (thetavals < binval + dtheta)
      if np.any(mask)!=False:hist[i] = np.nanmean(fulldata[mask], axis=0)
    np.savetxt( Outname(SimRun[0][1], 'Expectation_'+name, '.dat'), np.column_stack((thetabins, hist)) )

def PlotExpectations(OutDir):
  alldat={}
  for dirpath, dirnames, filenames in os.walk(OutDir):
    for f in [h for h in filenames if ("_.dat" in h and "Expectation_" in h)]:
      name=dirpath+'/'+f;paramname=dirpath+"/Params_"+"_".join(f.split("_")[2:]).replace(".dat",".obj")
      dat = np.loadtxt(name);Params = LoadParams(paramname)
      keyword=f.split("Expectation_")[1].split('_')[0]
      print name
      try:alldat[keyword].append([dat, Params])
      except:alldat[keyword]=[[dat, Params]]
      
      InitPlot(r'membrane gap [nm]', Codes[keyword], Params)
      if keyword=="TotPop":
        NPopsofz = nomemflux.ConstructPnums(Params)
        Ntot=np.array([NPopsofz(z) for z in dat[:, 0]])
        PlotPlot(dat[:, 0], Params.nodes*Ntot[:, -1], '', '--', 7, 2.2, r'$V_{\mathrm{eff}}$')
      
      for i in range(1, dat.shape[1]): PlotPlot(dat[:, 0], dat[:, i], '', '-', 7, 2.2, None)
      PlotLengths(Params)
      PlotSave(Params, 'Expectation_'+keyword)
  
  for k in alldat.keys():
    newcode=Codes[k][:-1]+'/'+Codes[k][1:-1]+r'_{{z=\infty}}$'
    InitPlot(r'membrane gap [nm]', newcode, None)
    for a in alldat[k]:
      dat=a[0];Params=a[1]
      for i in range(1, dat.shape[1]): 
        label=r'$\mu={:.1f},\, \kappa='.format(Params.densities[0])
        if Params.kc>100.0:
          label=label+'\infty$';lt='--'
        else:
          label=label+'{}$'.format(Params.kc);lt='-'
        norm=max(1.0, np.average(dat[-5:,i]))
        PlotPlot(dat[:, 0], dat[:, i]/norm, '', lt, 7, 2.2, label)
    PlotLengths(Params)
    PlotSave(Params, "SummaryRel_"+k)
    
  for k in alldat.keys():
    InitPlot(r'membrane gap [nm]', Codes[k], None)
    for a in alldat[k]:
      dat=a[0];Params=a[1]
      for i in range(1, dat.shape[1]): 
        label=r'$\mu={:.1f},\, \kappa='.format(Params.densities[0])
        if Params.kc>100.0:
          label=label+'\infty$';lt='--'
        else:
          label=label+'{}$'.format(Params.kc);lt='-'
        norm=max(1.0, np.nanmax(dat[:,i]))
        PlotPlot(dat[:, 0], dat[:, i], '', lt, 7, 2.2, label)
    PlotLengths(Params)
    PlotSave(Params, "Summary_"+k)
    
    
    
    
    
    
    
    
    
    
    
    
    
    

