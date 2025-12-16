import matplotlib.pyplot as plt
import numpy as np
import pickle
import effectivepot as effpot
import os

Codes={
       'MeshE'   : r'$\langle E_{{\mathrm{{Mem}}}}\rangle$', 
       'TotArea' : r'$\langle A_{{\mathrm{{Mem}}}}\rangle/A_0$', 
       'PartE'   : r'$\langle \beta U_{{B}}\rangle$',
       'TotPop'  : r'$\langle N_{{\mathrm{{Tot}}}}\rangle$',
       'rhoAve'  : r'$\langle \rho \rangle [1/ \mu m^{{2}}]$',
       'Pops'    : r'$\langle N_{i}\rangle$',
       'rho'     : r'$\langle \rho_i \rangle [1/ \mu m^{{2}}]$',
      }

def Outname(Params, prefix, suffix):return Params.OutDir+'/'+prefix+'_'+Params.filename+'_'+suffix

def InitPlot(xlabel, ylabel, Params):
  plt.clf()
  plt.rc('text', usetex=True);plt.rc('font', family='serif', size=18);
  plt.xlabel(xlabel);plt.ylabel(ylabel)
  if Params!=None:plt.title(r'PTypes$={}$ $\beta\kappa_B={} A_0 ={}\,\mu m^2$'.format(Params.numptypes, Params.kc, Params.A/1e6) )
  
def PlotPlot(x, y, mark, line, msize, lw, label): plt.plot(x, y, marker=mark, markersize=msize, linestyle=line, linewidth=lw, label=label)

def PlotLengths(Params):
  for i in range(Params.numptypes):
    if Params.bindingstrs[i]<1e-6:
      plt.axvline(x=(Params.lengths[i] ), 
                  color='k', linestyle=':', 
                  label=r'$L_{}={}$'.format(i+1, (Params.lengths[i] ))
                )
      #plt.axvline(x=2.0*(Params.lengths[i] + Params.rads[i]), #color='k', linestyle='--', #label=r'$2 L_{}$'.format(i+1) )
      top = Params.lengths[i] + np.sqrt(Params.heightvar)
      bot = Params.lengths[i] - np.sqrt(Params.heightvar)
      plt.axvline(x= top, color='green', lw=0, label=r'$2\sigma_h = {:.1f}$'.format(2.0*np.sqrt(Params.heightvar)))
      #plt.axvline(x=bot, color='green', alpha=0.5)
      ymin, ymax = plt.ylim()
      plt.fill_between(np.linspace(bot, top, 100), ymin, ymax, alpha=0.3, facecolor='green')
      plt.ylim((ymin, ymax))
    
    else:
      length = Params.lengths[i] + Params.lengths[i+Params.numptypes] + Params.bindingran[i]/2
      plt.axvline(x=length, 
                  color='k', linestyle=':', 
                  label=r'$L_{}^{{\mathrm{{rest}}}}={}$'.format(i+1, length)
                )
    
def PlotSave(Params, name, lg=True):
  if lg==True: plt.legend(prop={'size': 12})
  plt.tight_layout()
  plt.savefig(Outname(Params, name, ".png"))

#Variety of plots
  
def PlotFreeEnergy(Data):
  FEDat, Params = Data
  InitPlot(r'membrane gap [nm]', r'$\beta \Delta F$', Params)
  PlotPlot(FEDat["thetavals"], FEDat["FE"], '', '-', None, 2.2, None)
  veff = effpot.GetFe(Params, FEDat["thetavals"])
  PlotPlot(FEDat["thetavals"], np.nanmax(np.abs(FEDat["FE"]))*veff/np.amax(np.abs(veff)), '', '--', None, 2.2, r'$V_{\mathrm{eff}}$')
  PlotLengths(Params)
  PlotSave(Params, 'FE_thetaval')

def PlotPopulations(Data, Area):
  FEDat, Params = Data
  thetavals = FEDat["thetavals"]
  
  InitPlot(r'membrane gap [nm]', r'$N$', Params)
  PlotPlot(thetavals, FEDat["TotPop"], '', '-', None, 2.2, r'$N_T$')
  if Params.numptypes>1:
    for i in range(1, Params.numptypes+1):
      t = "NTop"+str(i);b = "NBot"+str(i)
      PlotPlot(thetavals, FEDat[t] + FEDat[b], '', '-', None, 2.2, r'$N_{}$'.format(i))
  PlotLengths(Params)
  
  #computepnums = effpot.ConstructPnums(Params)
  #pnums = []
  #for t in thetavals:pnums.append(computepnums(t))
  #pnums = np.array(pnums)
  #plt.gca().set_prop_cycle(None)
  #PlotPlot(thetavals, pnums[:, -1], '', '--', None, 2.2, r'$N_T$')
  #if Params.numptypes>1:
    #for i in range(Params.numptypes):
      #t = Params.numptypes + i + 1
      #PlotPlot(thetavals, pnums[:, i] + pnums[:, t], '', '--', None, 2.2, r'$N_{}$'.format(i+1))
  PlotSave(Params, 'Pops')

def PlotDatafiles(OutDir):
  for Dname in filter(lambda x: '_.obj' in x, os.listdir(OutDir)):
    Data = pickle.load( open(OutDir+'/'+Dname, 'rb') )
    PlotFreeEnergy(Data)
    PlotPopulations(Data)
  



















