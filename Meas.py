import numpy as np
import pickle
import pandas as pd
from Tools import Tools
import traceback
from dataclasses import dataclass

@dataclass
class settings:
    la: list
    cn: list
    SubLeadCoefficients: list
    TheoryOrder: list
    
    # Constants from fit.config
    rho2: float = -0.05
    mB  : float = 5.279
    La2 : float = 0.12
    N0  : float = 794.0705566209606
    VtbVts: float = 0.04129300965
    C2C7: float = -2.109974309809888
    C2C2: float = 1.1129978970144285
    C8C7: float = 0.2760118962325344
    C8C8: float = 0.019045641715469835
    C2C8: float = -0.29118900512628015

class Meas:
    def __init__(self):

        #Theory Dictionaries
        self.theory_expx3 = Tools.pickle_to_dict('theory/theory_dictionary_expx3')
        self.theory_SSF = Tools.pickle_to_dict('theory/theory_dictionary_expx3')

        #Dictionary with experimental data
        self.exp_data = Tools.pickle_to_dict('data/exp_data')
        
        #Histogram Dictionary from Root Data
        self.hist_nom = Tools.pickle_to_dict("theory/hist_nom")

        #Moments Dictionary
        self.Fmn_moments = Tools.pickle_to_dict("theory/Fmn_moments")

        self.mb

        return
    
    # Prefactor for leading and subleading theory
    def TheoryPrefactor(self, TheoryOrder:str, norm:float):

        value = settings.N0

        if ('22' in TheoryOrder):
            value *= settings.C2C2
            value /= norm
            value *= np.pow(settings.VtbVts * self.mb, 2.0)
        elif ('SSF27' in TheoryOrder):
            value *= settings.C2C7 / np.sqrt(norm) * settings.VtbVts * self.mb
            value *= settings.La2 / self.mb
        elif ('27' in TheoryOrder):
            value *= settings.C2C7 / np.sqrt(norm) * settings.VtbVts * self.mb
        elif ('28' in TheoryOrder):
            value *= settings.C2C8 / norm * np.pow(settings.VtbVts * self.mb, 2.0)
        elif ('88' in TheoryOrder):
            value *= settings.C8C8
            value /= norm
            value *= np.pow(settings.VtbVts * self.mb, 2.0)
        elif ('78'  in TheoryOrder):
            value *= settings.C8C7
            value /= np.sqrt(norm)
            value *= settings.VtbVts * self.mb
        
        return value



    # Returns a prediction array for one specific experiment with one specific lambda
    #TODO: mid is specific for Prediction?
    def BsgPrediction(self, key:str, mid:str, end:str, c_n_params, norm):      
        pred = np.array([])
        for i in range(np.size(self.exp_data[key]['dBFs'])):
            value = 0
            for j in range(np.size(c_n_params)):
                for k in range(np.size(c_n_params)):
                    value += self.theory_expx3[key][mid]['la'+end]['Values'][i][j][k] * c_n_params[j] * c_n_params[k] 
                    value *= 1/(self.theory_expx3[key][mid]['la'+end]['lambda']*norm)
            pred = np.append(pred, value)
        return pred 


    #TODO: Mid always the same with subleading and so leading?
    def BsgSubLeadingPrediction(self, key:str, mid:str, end:str, c_n_params, norm):
        pred = np.array([])
        for i in range(np.size(self.exp_data[key]['dBFs'])):
            value = 0
            for j in range(np.size(c_n_params)):
                value +=  self.theory_expx3[key][mid]['la'+end]['Values'][i][0][j] * c_n_params[j]
                value *= 1/(self.theory_expx3[key][mid]['la'+end]['lambda']*norm)
            pred = np.append(pred, value)

        return pred
    
    def SubLeadPars_read(self, key:str, end:str): #TODO: eigentlich Quatsch diese Funktion
        return self.theory_SSF[key]['SSF27'][end]['Values']
    
    # dds Zoltan's SF into the prediction
    # converts c_n's into subleading coefficients

    #TODO: c_n_params not used, also in c++ code?
    def SubLeadPars(self, c_n_params, d2: float, opt:int):
        
        #TODO: Find the real values and how/when/where to read them in EDIT: Found them in fit.config, are they correct? Is there another file
        Rho2 = -0.05
        mB = 5.279
        mb = self.mb1SPrediction(c_n_params)
        la = self.lambda11SPrediction(c_n_params)
        Lambda2 = 0.12 


        if(opt == 1):
            x = (0.6514810199386504 * (mB - mb) - 0.8686413599182005 * la + 0.3257405099693252 * Rho2 / Lambda2) / (1.8531015678254943 * la - d2 * la)
        elif(opt == 2):
            x = (0.4722982832332954 * (mB - mb) - 0.5667579398799545 * la + 0.2361491416166477 * Rho2 / Lambda2) / (1.3695121740357048 * la - d2 * la)
        else:
            print('Unknown subleading shape function -- please specify Meas.SubLeadPars()')

        #calculating d0, d1 and d2 and returning them in a numpy array
        return np.array([1.-x, x* (1-d2), x*d2])
    
    
    # Calculate shape fuction moment of order 'order'
    def Moment(self, c_n_params, order: int): #TODO: Test with Values for c_n
        if (order > np.size(self.Fmn_moments['expx3']['Moment'])): print('Meas.Moment(): requested an order of moment, that cannot be calculated')
        value = 0
        for i in range(np.size(c_n_params)):
            for j in range(np.size(c_n_params)):
                value += self.Fmn_moments['expx3']['Moment'][order][i][j] * c_n_params[i] * c_n_params[j] #TODO: Dont know if indices are correct
        return value
    
    # Calculate mb^1S for a given set of cn's
    def mb1SPrediction(self, c_n_params): #TODO: Test it
        mb = 0
        M1 = self.Moment(c_n_params, 1)
        M2 = self.Moment(c_n_params, 2)
        M3 = self.Moment(c_n_params, 3)

        Lambda2 = settings.La2
        mB = settings.mB
        Rho2 = settings.rho2

        p = -18 * Lambda2 - 25 * M1 * M1 + 36 * M2 - 22 * M1 * mB + 11 * mB * mB
        q = 5 * M1 * (27 * Lambda2 + 25 * M1 * M1 - 54 * M2) - 3 * (45 * Lambda2 - 55 * M1 * M1 + 72 * M2) * mB + 51 * M1 * mB * mB - 17 * mB * mB * mB + 162 * (M3 + Rho2)
        u = np.pow((-q + np.sqrt(q * q + p * p * p)), 1. / 3.)
        mb = 1. / 6. * (5 * (mB - M1) - p / u + u)

        return mb
    
    # Calculate lambda for a given set of cn's
    def lambda11SPrediction(self, c_n_params):
        M1 = self.Moment(c_n_params, 1)
        lambda2 = settings.La2
        mB = settings.mB

        mb = self.mb1SPrediction(c_n_params)
        lambda1 = 3. * lambda2 + 2. * mb * (mB - M1 - mb)

        return lambda1


    

#print(Meas.SubLeadPars('babar_hadtag', '105')['Values'][0])

#print(h.BsgPrediction('babar_hadtag','NNLLNNLO', '055', test_fit_results, test_norm))

#print(h.theory_expx3['babar_hadtag']['NNLLNNLO']['la'+'055']['Values'][0][0][0]) #NOTE: ...[#bin][line][column]


#Schlachtplan
# Read in moment files EDIT: DONE, BUT TO LESS NUMBERS 
# Make Moment function EDIT: DONE, BUT NOT TESTED
# Make mb1SPrediction EDIT: DONE, BUT NOT TESTED
# Make lambda11Prediction EDIT: DONE, BUT NOT TESTED
# Finish SubLeadPars EDIT: DONE, BUT NOT TESTED
# Implement TheoryPrefactor
# _mb must be a changable variable in BsgPrediction
# Find out what FindIndex does, and which data it uses
# Implement large BsgPrediction function 