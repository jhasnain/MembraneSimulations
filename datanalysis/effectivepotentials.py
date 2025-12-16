import numpy as np
import matplotlib.pyplot as plt

class paramssmall:
  def __init__(self):
    self.numptypes = 1
    self.Lx        = 320.0
    self.Ly        = 320.0
    self.nx        = 80
    self.ny        = 80

    self.nodes     = self.nx*self.ny
    self.dx        = self.Lx/self.nx
    self.dy        = self.Ly/self.ny

    self.densities = [0.0, 0.0]
    self.fugacities= self.convertdensitytofugacity(self.densities)

    if self.numptypes>0:
      self.lengths        = [0.0, 0.0]
      self.bindingstrs    = [8.0]
      self.bindingran     = [4.0]

      self.diams = [self.dx for i in range(2*self.numptypes)]
      self.rads  = [d/2 for d in self.diams]

  def convertdensitytofugacity(self, densities):
    alpha  = 1e6*self.nodes/(self.Lx*self.Ly)
    zbar   = np.sum(densities[:self.numptypes]);
    zbar   = zbar/(alpha - zbar)
    fugacities = [ d*(1.0 + zbar)/alpha for d in self.densities[:self.numptypes] ]
    
    zbar   = np.sum(densities[self.numptypes:]);
    zbar   = zbar/(alpha - zbar)
    fugacities +=[ d*(1.0 + zbar)/alpha for d in self.densities[self.numptypes:] ]
    return fugacities
    

def StateProbs(Params, Diams=True):
  val = (1 if Diams==True else 0)
  Veff = ConstructVeff(Params, Diams)
  pnums = Params.numptypes
  
  lengthsbot = np.array(Params.lengths[:pnums])
  radsbot    = np.array(Params.rads[:pnums])
  zbot       = np.array(Params.fugacities[:pnums])
  
  lengthstop = np.array(Params.lengths[pnums:])
  radstop    = np.array(Params.rads[pnums:])
  ztop       = np.array(Params.fugacities[pnums:])
  
  bstr = np.array(Params.bindingstrs)
  brange = np.array(Params.bindingran)
  
  lensums = lengthsbot[:, np.newaxis] + lengthstop[np.newaxis, :]
  radsums = radsbot[:, np.newaxis] + radstop[np.newaxis, :]
  
  VeffZero = np.log((1.0 + np.sum(zbot))*(1.0 + np.sum(ztop)))
  expVeff = lambda z: np.exp(Veff(z) - VeffZero)
  
  DoubleSpins = lambda sbot, stop, z: np.where( (sbot==stop) & (bstr[sbot]>1e-6),
                                                np.where( z < lensums [sbot, stop], 
                                                          0.0, zbot[sbot]*ztop[stop]*np.where((z - lensums [sbot, stop]) < brange[sbot],  np.exp(bstr[sbot]), 1.0)
                                                        ),
                                                np.where( z > (lensums + val*radsums)[sbot, stop],
                                                          zbot[sbot]*ztop[stop], 0.0
                                                        )
                                              )
  
  SingleSpinsBot = lambda sbot, z: np.where( z > (lengthsbot + val*radsbot)[sbot], zbot[sbot], 0.0)
  SingleSpinsTop = lambda stop, z: np.where( z > (lengthstop + val*radstop)[stop], ztop[stop], 0.0)
  
  return lambda z: [ SingleSpinsBot(sbot, z)*expVeff(z)    for sbot in range(pnums) ] + \
                   [ SingleSpinsTop(stop, z)*expVeff(z)    for stop in range(pnums) ] + \
                   [ DoubleSpins(sbot, stop, z)*expVeff(z) for sbot in range(pnums) for stop in range(pnums) ] + \
                   [ expVeff(z) ]

def ConstructPnums(Params, Diams=True):
  val = (1 if Diams==True else 0)
  Veff = ConstructVeff(Params, Diams)
  pnums = Params.numptypes
  
  lengthsbot=np.array(Params.lengths[:pnums])
  radsbot = np.array(Params.rads[:pnums])
  zbot = np.exp(np.array(Params.muvals[:pnums]))
  
  lengthstop=np.array(Params.lengths[pnums:])
  radstop = np.array(Params.rads[pnums:])
  ztop = np.exp(np.array(Params.muvals[pnums:]))
  
  bstr = np.array(Params.bindingstrs)
  brange = np.array(Params.bindingran)
  
  lensums = lengthsbot[:, np.newaxis] + lengthstop[np.newaxis, :]
  radsums = radsbot[:, np.newaxis] + radstop[np.newaxis, :]
  
  VeffZero = np.log((1.0 + np.sum(zbot))*(1.0 + np.sum(ztop)))
  expVeff = lambda z: np.exp(Veff(z) - VeffZero)
  
  doublens = lambda z, i, opt: np.sum(np.where(
                                  opt==0,
                                  np.where(z>(lensums + val*radsums)[i, :], ztop, 0.0),
                                  np.where(z>(lensums + val*radsums)[:, i], zbot, 0.0)
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
                                    np.where(z>(lengthsbot[i] + val*radsbot[i]), 1.0, 0.0),
                                    np.where(z>(lengthstop[i] + val*radstop[i]), 1.0, 0.0)       
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

def GetFe(Params, zvals, Diams=True):
  Veff=ConstructVeff(Params, False)
  return Params.nodes*np.array([Veff(z) for z in zvals])

def ConstructVeff(Params, Diams=True):
  val = (1 if Diams==True else 0)
  pnums = Params.numptypes
  
  lengthsbot = np.array(Params.lengths[:pnums])
  radsbot    = np.array(Params.rads[:pnums])
  zbot       = np.array(Params.fugacities[:pnums])
  
  lengthstop = np.array(Params.lengths[pnums:])
  radstop    = np.array(Params.rads[pnums:])
  ztop       = np.array(Params.fugacities[pnums:])
  
  bstr = np.array(Params.bindingstrs)
  brange = np.array(Params.bindingran)
  
  lensums = lengthsbot[:, np.newaxis] + lengthstop[np.newaxis, :]
  radsums = radsbot[:, np.newaxis] + radstop[np.newaxis, :]
  zprods  = zbot   [:, np.newaxis] * ztop   [np.newaxis, :]
  
  singlespins = lambda z: np.sum(zbot[z>(lengthsbot+val*radsbot)]) + np.sum(ztop[z>(lengthstop+val*radstop)])
  
  doublens =  lambda z: np.sum(np.where( z < lensums , 0.0, 
                                  zprods*np.where( z > (lensums + val*radsums), 1.0, 0.0 )
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
                                   
  Veff = lambda z: np.log((1.0 + np.sum(zbot))*(1.0 + np.sum(ztop))/(1.0 + singlespins(z) + doublens(z) + bindterms(z)))
  return Veff
  
if __name__ == "__main__":
  Params=paramssmall()  
    
  zvals = np.linspace(20.0, 0.0, 1000)
  
  veffFunc=ConstructVeff(Params)
  
  pot=np.array([veffFunc(z) for z in zvals])
  plt.clf()
  plt.plot(zvals, pot, marker='x')
  plt.show()


