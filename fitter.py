import numpy as np
from Tools import Tools
import traceback
from dataclasses import dataclass
from iminuit import Minuit, cost
import iminuit
from Meas import settings
from Meas import Meas

class Fitter:
    
    #start_pars: Which start parameters will be used
    #number_pars: How many parameters will be fitted, so how many of start_pars will be used
    #no_meas: how many measurements will be included in the fit
    #with_minos: Calculate the fit with or without minos
     def DoFit(start_pars, number_pars = 4, no_meas = 2, with_minos = True):
        
        sets = settings()
        sets.KeyOrder = sets.KeyOrder[0:no_meas]
        print('Fit for: '+ str(sets.KeyOrder))

        mes = Meas()

        m = Minuit(mes.Chisq, start_pars[0:number_pars])
        value = mes.Chisq(np.array(m.values))
        m.print_level = 1

        m.tol = 1000 #sets EDM_max = 0.1
        
        m.errordef = 0.0001
        #m.simplex()
        m.tol = 500 # sets EDM_max = 0.0001
        m.migrad(ncall= 100000, use_simplex= True)
        m.hesse()
        

        chisq = mes.chisq
        mb = mes.mb

        if with_minos == True:
            m.minos()
        
        value = mes.Chisq(np.array(m.values))
            
        return m, mes, chisq, mb
     
     #Calculates fits with several numbers of fitting Parameters for a specific number of measurements
     def collect_fit_data(start_pars, start_numb_pars, end_numb_pars, no_meas = 2, with_minos = True):
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
'''
start_pars = np.array([1, 0.00506919, 0.0, 0.0798100, 0.0870341, 0.0250290, 0.0])

collected_fits ={
    #'4': Fitter.collect_fit_data(start_pars, 3,7,no_meas = 4, with_minos = True),
    '3': Fitter.collect_fit_data(start_pars, 3,7,no_meas = 3, with_minos=True),
    '2': Fitter.collect_fit_data(start_pars, 3,7,no_meas = 2, with_minos=True)
}


Tools.store_in_pickle(collected_fits, 'fit/collected_fits')

'''