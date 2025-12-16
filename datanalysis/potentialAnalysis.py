import matplotlib.pyplot as plt
import numpy as np
import pickle
import effectivepot as effpot
import os
plt.rc('text', usetex=True);plt.rc('font', family='serif', size=16);

def PlotLenghts(Params):
  for i in range(Params.numptypes):
      if Params.bindingstrs[i]<1e-6:
        plt.axvline(x=(Params.lengths[i] ), 
                    color='k', linestyle='--', 
                    label=r'$L_{}={}$'.format("B", (Params.lengths[i] ))
                  )
      else:
        length = Params.lengths[i] + Params.lengths[i+Params.numptypes] + Params.bindingran[i]/2
        plt.axvline(x=length, 
                    color='k', linestyle=':', 
                    label=r'$L_{}^{{\mathrm{{rest}}}}={}$'.format(i+1, length)
                  )
        top = length + Params.bindingran[i]/2
        bot = length - Params.bindingran[i]/2
        plt.axvline(x= top, color='green', lw=0)
        #plt.axvline(x=bot, color='green', alpha=0.5)
        ymin, ymax = plt.ylim()
        plt.fill_between(np.linspace(bot, top, 100), ymin, ymax, alpha=0.3, facecolor='green')
        #plt.ylim((ymin, ymax))

def ComputeChemPot(densities, alpha):
  zbar = np.sum(densities)
  zbar = zbar/(alpha-zbar)
  muvals = [np.log((1.0 + zbar)*d/alpha) for d in densities]
  return muvals

class params:
    def __init__(self):
      self.Lx        = 320.0
      self.Ly        = 320.0
      self.nx        = 80
      self.ny        = 80
      
      self.numptypes = 2
      if self.numptypes>0:
        self.densities      = [1000, 100, 1000, 100]
        self.lengths        = [17.0, 6.0, 17.0, 6.0]
        self.bindingstrs    = [0.0, 6.0]
        self.bindingran     = [4.0, 4.0]
      
      #self.numptypes = 1
      #if self.numptypes>0:
        #self.densities      = [1000, 1000]
        #self.lengths        = [6.0, 6.0]
        #self.bindingstrs    = [0.0]
        #self.bindingran     = [4.0]
      
      self.nodes     = self.nx*self.ny
      self.dx        = self.Lx/self.nx
      self.dy        = self.Ly/self.ny
      self.A         = self.Lx*self.Ly
      self.alpha     = 1e6*self.nodes/self.A
      self.muvals    = np.array(ComputeChemPot(self.densities[:self.numptypes], self.alpha) + ComputeChemPot(self.densities[self.numptypes:], self.alpha))
      
    def newdense(self, densities):
      self.densities = densities
      self.muvals    = np.array(ComputeChemPot(self.densities[:self.numptypes], self.alpha) + ComputeChemPot(self.densities[self.numptypes:], self.alpha))
      
def PlotFreeEnergies(Params, zvals, blockerdensities):
  Params.lengths = [17.0, 6.0, 17.0, 6.0]
  plt.clf()
  plt.xlabel(r"\textrm{membrane gap [nm]}");plt.ylabel(r"$\beta V_{\textrm{eff}}$")
  for d in blockerdensities:
    Params.newdense( [d, 100, d, 100] );
    veff = effpot.GetFe(Params, zvals)  
    plt.plot(zvals, veff, label = r"$\rho_B = {:.0f} $".format(d))
  PlotLenghts(Params)
  plt.legend()
  plt.tight_layout()
  plt.savefig("Longblockers.pdf")


  Params.lengths = [14.0, 6.0, 14.0, 6.0]
  plt.clf()
  plt.xlabel(r"\textrm{membrane gap [nm]}");plt.ylabel(r"$\beta V_{\textrm{eff}}$")
  for d in blockerdensities:
    Params.newdense( [d, 100, d, 100] );
    veff = effpot.GetFe(Params, zvals)  
    plt.plot(zvals, veff, label = r"$\rho_B = {:.0f} $".format(d))
  PlotLenghts(Params)
  plt.legend()
  plt.tight_layout()
  plt.savefig("Sameblockers.pdf")

  Params.lengths = [10.0, 6.0, 10.0, 6.0]
  plt.clf()
  plt.xlabel(r"\textrm{membrane gap [nm]}");plt.ylabel(r"$\beta V_{\textrm{eff}}$")
  for d in blockerdensities:
    Params.newdense( [d, 100, d, 100] );
    veff = effpot.GetFe(Params, zvals)  
    plt.plot(zvals, veff, label = r"$\rho_B = {:.0f} $".format(d))
  PlotLenghts(Params)
  plt.legend()
  plt.tight_layout()
  plt.savefig("Shortblockers.pdf")

def PlotParticlePopulations(Params, zvals, blockerdensities):
  Params.lengths = [17.0, 6.0, 17.0, 6.0]
  plt.clf()
  plt.xlabel(r"\textrm{membrane gap [nm]}");plt.ylabel(r"$\rho [\mathrm{\mu m ^{-1}}]$")
  for d in blockerdensities:
    Params.newdense( [d, 100, d, 100] );
    p = effpot.ConstructPnums(Params)
    particlepops  = np.array([ p(z) for z in zvals])
    plt.plot(zvals, Params.alpha*(particlepops[:, 1] + particlepops[:, 4])/2 , label = r"$\rho_B = {:.0f} $".format(d))
  PlotLenghts(Params)
  plt.legend()
  plt.tight_layout()
  plt.savefig("N_Longblockers.pdf")

  Params.lengths = [14.0, 6.0, 14.0, 6.0]
  plt.clf()
  plt.xlabel(r"\textrm{membrane gap [nm]}");plt.ylabel(r"$\rho [\mathrm{\mu m ^{-1}}]$")
  for d in blockerdensities:
    Params.newdense( [d, 100, d, 100] );
    p = effpot.ConstructPnums(Params)
    particlepops  = np.array([ p(z) for z in zvals])
    plt.plot(zvals, Params.alpha*(particlepops[:, 1] + particlepops[:, 4])/2 , label = r"$\rho_B = {:.0f} $".format(d))
  PlotLenghts(Params)
  plt.legend()
  plt.tight_layout()
  plt.savefig("N_Sameblockers.pdf")

  Params.lengths = [10.0, 6.0, 10.0, 6.0]
  plt.clf()
  plt.xlabel(r"\textrm{membrane gap [nm]}");plt.ylabel(r"$\rho [\mathrm{\mu m ^{-1}}]$")
  for d in blockerdensities:
    Params.newdense( [d, 100, d, 100] );
    p = effpot.ConstructPnums(Params)
    particlepops  = np.array([ p(z) for z in zvals])
    plt.plot(zvals, Params.alpha*(particlepops[:, 1] + particlepops[:, 4])/2 , label = r"$\rho_B = {:.0f} $".format(d))
  PlotLenghts(Params)
  plt.legend()
  plt.tight_layout()
  plt.savefig("N_Shortblockers.pdf")
        
Params=params()  

zvals    = np.linspace(25.0, 11.0, 200)
blockerdensities = np.linspace(100.0, 2000.0, 5)

PlotFreeEnergies(Params, zvals, blockerdensities)
PlotParticlePopulations(Params, zvals, blockerdensities)










