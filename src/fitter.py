import numpy as np
from Tools import Tools
from iminuit import Minuit
from Meas import settings
from Meas import Theory

class Fitter:
     
     def __init__(self):
         
         m: Minuit

         mb: float
         chisq: float

         theo: Theory

         NumbOfPar: int

         Lambda: float

    
    #NumbPar: How many parameters will be fitted, so how many of start_pars will be used
    #with_minos: Calculate the fit with or without minos
     def DoSingleFit(self, NumbPar: int = 3, with_minos: bool = True):

        self.NumbOfPar = NumbPar

        print('Fit for: '+ str(settings.KeyOrder) + ' using %d parameters' % NumbPar)

        print('b')

        mes = Theory()

        n = Minuit(mes.Chisq, settings.config["StartValues"][0:NumbPar], name = settings.config["FitVars"][0:NumbPar])
        #value = mes.Chisq(np.array(n.values))
        #n.print_level = 1

        #n.tol = 0.01 #1000 #sets EDM_max = 0.1
        
        n.errordef = 0.0001
        #n.simplex()
        n.tol = 0.1 #500 # sets EDM_max = 0.0001
        n.migrad(ncall= 100000, use_simplex= True)
        n.hesse()
        
        self.chisq = mes.chisq
        self.mb = mes.mb

        if with_minos == True:
            n.minos()
        
        #value = mes.Chisq(np.array(n.values))

        self.m  = n
        self.theo = mes

        self.Lambda = Tools.StrToLambda(settings.BasisExpansion)
            
        return