import numpy as np
import matplotlib.pyplot as plt

def StateProbs(Params):
  Veff = ConstructVeff(Params)
  pnums = Params.numptypes
  
  lengthsbot=np.array(Params.lengths[:pnums])
  zbot = np.exp(np.array(Params.muvals[:pnums]))
  
  lengthstop=np.array(Params.lengths[pnums:])
  ztop = np.exp(np.array(Params.muvals[pnums:]))
  
  bstr = np.array(Params.bindingstrs)
  brange = np.array(Params.bindingran)
  
  lensums = lengthsbot[:, np.newaxis] + lengthstop[np.newaxis, :]
  
  VeffZero = np.log((1.0 + np.sum(zbot))*(1.0 + np.sum(ztop)))
  expVeff = lambda z: np.exp(Veff(z) - VeffZero)
  
  DoubleSpins = lambda sbot, stop, z: np.where( (sbot==stop) & (bstr[sbot]>1e-6),
                                                np.where( z < lensums [sbot, stop], 
                                                          0.0, zbot[sbot]*ztop[stop]*np.where((z - lensums [sbot, stop]) < brange[sbot],  np.exp(bstr[sbot]), 1.0)
                                                        ),
                                                np.where( z > lensums[sbot, stop],
                                                         zbot[sbot]*ztop[stop], 0.0
                                                        )
                                              )
  
  SingleSpinsBot = lambda sbot, z: np.where( z > lengthsbot[sbot], zbot[sbot], 0.0)
  SingleSpinsTop = lambda stop, z: np.where( z > lengthstop[stop], ztop[stop], 0.0)
  
  return lambda z: [ SingleSpinsBot(sbot, z)*expVeff(z) for sbot in range(pnums) ] + \
                   [ SingleSpinsTop(stop, z)*expVeff(z) for stop in range(pnums) ] + \
                   [ DoubleSpins(sbot, stop, z)*expVeff(z) for sbot in range(pnums) for stop in range(pnums) ] + \
                   [ expVeff(z) ]

def ConstructPnums(Params):
  Veff = ConstructVeff(Params)
  pnums = Params.numptypes
  
  lengthsbot=np.array(Params.lengths[:pnums])
  zbot = np.exp(np.array(Params.muvals[:pnums]))
  
  lengthstop=np.array(Params.lengths[pnums:])
  ztop = np.exp(np.array(Params.muvals[pnums:]))
  
  bstr = np.array(Params.bindingstrs)
  brange = np.array(Params.bindingran)
  
  lensums = lengthsbot[:, np.newaxis] + lengthstop[np.newaxis, :]
  
  VeffZero = np.log((1.0 + np.sum(zbot))*(1.0 + np.sum(ztop)))
  expVeff = lambda z: np.exp(Veff(z) - VeffZero)
  
  doublens = lambda z, i, opt: np.sum(np.where(
                                  opt==0,
                                  np.where(z>lensums[i, :], ztop, 0.0),
                                  np.where(z>lensums[:, i], zbot, 0.0)
                              ))
  
  bindterms = lambda z, i, opt: np.sum(np.where(bstr[i]>1e-6,
                                   np.where(
                                      ((z - lensums[i,i]) > 0) & ( (z - lensums[i,i])< brange[i]),
                                      np.exp(bstr[i])*np.where( opt==0, ztop, zbot),
                                      0.0,
                                  ),
                                  0.0
                                ))
  
  singlens = lambda z, i, opt: np.where( 
                                    opt==0, 
                                    np.where(z>lengthsbot[i], 1.0, 0.0),
                                    np.where(z>lengthstop[i], 1.0, 0.0)       
                                  )
  
  ConstructN_i=lambda z, i, opt: np.where( 
                                    opt==0, 
                                    zbot[i]*(singlens(z, i, opt) + doublens(z, i, opt) + bindterms(z, i, opt))*expVeff(z),
                                    ztop[i]*(singlens(z, i, opt) + doublens(z, i, opt) + bindterms(z, i, opt))*expVeff(z)
                                  )
  
  PopnumsBot = lambda z: [ConstructN_i(z, i, 0) for i in range(Params.numptypes)] 
  PopnumsTop = lambda z: [ConstructN_i(z, i, 1) for i in range(Params.numptypes)] 
                   
  return lambda z: PopnumsBot(z) + [1.0 - sum(PopnumsBot(z)) ] +\
                   PopnumsTop(z) + [1.0 - sum(PopnumsTop(z)) ] +\
                   [ np.sum( PopnumsBot(z) + PopnumsTop(z)) ]
                  

def GetFe(Params, zvals):
  Veff=ConstructVeff(Params)
  return Params.nodes*np.array([Veff(z) for z in zvals])

def ConstructVeff(Params):
  pnums = Params.numptypes
  lengthsbot=np.array(Params.lengths[:pnums])
  zbot = np.exp(np.array(Params.muvals[:pnums]))
  
  lengthstop=np.array(Params.lengths[pnums:])
  ztop = np.exp(np.array(Params.muvals[pnums:]))
  
  bstr = np.array(Params.bindingstrs)
  brange = np.array(Params.bindingran)
  
  lensums = lengthsbot[:, np.newaxis] + lengthstop[np.newaxis, :]
  zprods  = zbot   [:, np.newaxis] * ztop   [np.newaxis, :]
  
  singlespins = lambda z: np.sum(zbot[z>lengthsbot]) + np.sum(ztop[z>lengthstop])
  
  doublens =  lambda z: np.sum(np.where( z < lensums , 0.0, 
                                  zprods*np.where( z > lensums, 1.0, 0.0 )
                                )
                              )
  
  bindterms = lambda z: np.sum( 
                                np.where(bstr > 1e-6,
                                   np.where(
                                      ((z - lensums.diagonal())>0) & ((z - lensums.diagonal())< brange),
                                      np.exp(bstr)*zprods.diagonal(),
                                      0.0,
                                  ),
                                  0.0
                                )
                              )
                                   
                                   
  Veff = lambda z: -np.log(1.0 + singlespins(z) + doublens(z) + bindterms(z)) + np.log((1.0 + np.sum(zbot))*(1.0 + np.sum(ztop)))
  return Veff
  
  
if __name__ == "__main__":
  class params:
    def __init__(self):
      self.numptypes = 3
      self.muvals    = [-3.34109346, -3.34109346, -3.34109346, -3.34109346, -3.34109346, -3.34109346]
      self.Lx        = 320.0
      self.Ly        = 320.0
      self.nx        = 80
      self.ny        = 80
      
      self.nodes     = self.nx*self.ny
      self.dx        = self.Lx/self.nx
      self.dy        = self.Ly/self.ny
      
      if self.numptypes>0:
        self.lengths        = [4.0, 6.0, 8.0, 4.0, 6.0, 8.0]
        self.bindingstrs    = [6.0, 4.0, 2.0]
        self.bindingran     = [4.0, 4.0, 4.0]
    
  Params=params()  
    
  zvals = np.linspace(20.0, 0.0, 1000)
  
  veffFunc=ConstructVeff(Params)
  
  pot=np.array([veffFunc(z) for z in zvals])
  plt.clf()
  plt.plot(zvals, pot, marker='x')
  plt.show()


