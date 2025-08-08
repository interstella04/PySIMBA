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

        ########################
        print('Fit for: '+ str(settings.KeyOrder) + ' using %d parameters' % NumbPar)
        ########################

        mes = Theory()

        n = Minuit(mes.Chisq, settings.config["StartValues"][0:NumbPar], name = settings.config["FitVars"][0:NumbPar])
        value = mes.Chisq(np.array(n.values))
        n.print_level = 1

        n.tol = 1000 #sets EDM_max = 0.1
        
        n.errordef = 0.0001
        #n.simplex()
        n.tol = 500 # sets EDM_max = 0.0001
        n.migrad(ncall= 100000, use_simplex= True)
        n.hesse()
        
        self.chisq = mes.chisq
        self.mb = mes.mb

        if with_minos == True:
            n.minos()
        
        value = mes.Chisq(np.array(n.values))

        self.m  = n
        self.theo = mes

        self.Lambda = Tools.StrToLambda(settings.BasisExpansion)
            
        return
     
     '''
     #Calculates fits with several numbers of fitting Parameters for a specific number of measurements
     def CollectFitData(start_pars, start_numb_pars, end_numb_pars, no_meas = 2, with_minos = True):
        results = {}
        for i in range(start_numb_pars, end_numb_pars+1, 1):
            temp_dict = {}
            m, mes, chisq, mb = Fitter.DoFit(start_pars, i, no_meas=no_meas, with_minos=with_minos)
            temp_dict['an'] = np.array(m.values)
            temp_dict['errors'] = np.array(m.errors)
            temp_dict['Chisq'] = chisq
            temp_dict['mb'] = mb


            results['%d' % (i)] = temp_dict

        overall = {
            'results': results,
        }

        return overall
     

start_pars = np.array([1, 0.00506919, 0.0, 0.0798100, 0.0870341, 0.0250290, 0.0])

collected_fits ={
    #'4': Fitter.CollectFitData(start_pars, 3,7,no_meas = 4, with_minos = True),
    '3': Fitter.CollectFitData(start_pars, 3,7,no_meas = 3, with_minos=True),
    '2': Fitter.CollectFitData(start_pars, 3,7,no_meas = 2, with_minos=True)
}


Tools.store_in_pickle(collected_fits, 'fit/collected_fits')

'''