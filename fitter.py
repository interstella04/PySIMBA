import numpy as np
from Tools import Tools
import traceback
from dataclasses import dataclass
from iminuit import Minuit, cost
import iminuit
from Meas import settings
from Meas import Meas

class Fitter:
    #def __init__(self):
        
     def DoFit(start_pars, number_pars = 4, with_sub = True, no_meas = 2, with_minos = False):
        settings.KeyOrder = settings.KeyOrder[0:no_meas]
        print(settings.KeyOrder)
        mes = Meas()
        if with_sub == True:
            m = Minuit(mes.Chisq, start_pars[0:number_pars])
        else:
            m = Minuit(mes.Chisq_no_sub, start_pars[0:number_pars])

        m.migrad()
        m.hesse()

        if with_minos == True:
            m.minos()

        return m, mes
     
     def collect_fit_data(start_pars, start_numb_pars, end_numb_pars, no_meas = 2, with_minos = False):
        subleading = {}
        for i in range(start_numb_pars, end_numb_pars+1, 1):
            temp_dict = {}
            m, mes = Fitter.DoFit(start_pars, i, with_sub=True, no_meas=no_meas, with_minos=with_minos)
            temp_dict['an'] = np.array(m.values)
            temp_dict['errors'] = np.array(m.errors)
            temp_dict['Chisq'] = mes.chisq
            temp_dict['mb'] = mes.mb

            # New line to check global covariance matrix, can be deleted
            temp_dict['glob_cov'] = mes.glob_cov

            subleading['%d' % (i)] = temp_dict

        leading = {}
        for i in range(start_numb_pars, end_numb_pars+1, 1):
            other_dict = {}
            m, mes = Fitter.DoFit(start_pars, i, with_sub=False, no_meas=no_meas, with_minos=with_minos)
            other_dict['an'] = np.array(m.values)
            other_dict['errors'] = np.array(m.errors)
            other_dict['Chisq'] = mes.chisq
            other_dict['mb'] = mes.mb

            # New line to check global covariance matrix, can be deleted
            other_dict['glob_cov'] = mes.glob_cov

            leading['%d' % (i)] = other_dict

        overall = {
            'subleading': subleading,
            'leading': leading
        }

        #print(overall['subleading']['5']['an'])

        return overall

start_pars = np.array([1, 0.00506919, 0.0, 0.0798100, 0.0870341, 0.0250290, 0.0])

collected_fits ={
    '2': Fitter.collect_fit_data(start_pars, 4,7,no_meas = 4),
    '3': Fitter.collect_fit_data(start_pars, 4,7,no_meas = 3),
    '4': Fitter.collect_fit_data(start_pars, 4,7,no_meas = 2)
}

Tools.store_in_pickle(collected_fits, 'fit/collected_fits')

'''
just_sem = {
    '1': Fitter.collect_fit_data(start_pars, 4,7, no_meas=1)
}


Tools.store_in_pickle(just_sem, 'fit/just_sem')
'''