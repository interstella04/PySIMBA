import numpy as np
import pickle
import pandas as pd
from Tools import Tools
import traceback
from dataclasses import dataclass
from iminuit import Minuit, cost

@dataclass
class settings:
    la: list
    cn: list
    SubLeadCoefficients = [0.] # Don't really know, because in fit.config it is only SSF27_1055 and SSF27_2055 both are zero for d2
    TheoryOrder = ['NNLLNNLO', 'NS22NNLO', 'NS27NNLO', 'NS28NNLO', 'NS78NNLO', 'NS88NNLO']
    SubLeadTheoryOrder = 'SSF27_10575' # What is it? 
    FitVars = ['norm', 'a1', 'a2'] #List of Strings according to fit.config
    KeyOrder = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]
    BasisExpansion = '0575'
    SubLeadBasisExpansion = '10575' # ATTENTION, is in SubLeadingPred used

    DOMomentConstraints: bool #Should the Constraints be in the calculation or not
    
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
        self.theory_SSF = Tools.pickle_to_dict('theory/theory_dictionary_SSF27')

        #Dictionary with experimental data
        self.exp_data = Tools.pickle_to_dict('data/exp_data')

        #Dictionary with fit Configurations specific on dataset
        self.fit_config = Tools.pickle_to_dict('fit/fit_config')
        
        #Histogram Dictionary from Root Data
        self.hist_nom = Tools.pickle_to_dict("theory/hist_nom")

        #Moments Dictionary
        self.Fmn_moments = Tools.pickle_to_dict("theory/Fmn_moments")

        # Variables that will get a value in the process
        self.chisq: float
        self.mb: float 
        self.M1: float
        self.M2: float
        self.M3: float
        self.la: float

        return

    def BsgPrediction_full(self, key: str, end: str, c_n_params, norm):

        pred = np.zeros(np.size(self.exp_data[key]['dBFs']))
        self.mb = self.mb1SPrediction(c_n_params) #HERE mb should change values in the object

        for order in range(np.size(settings.TheoryOrder)):
            pred += self.BsgPrediction(key, settings.TheoryOrder[order], end, c_n_params, norm ) * self.TheoryPrefactor(settings.TheoryOrder[order], norm)
        
        for order in range(np.size(settings.SubLeadTheoryOrder)):
            pred += self.BsgSubLeadingPrediction(key, self.SubLeadPars(c_n_params, settings.SubLeadCoefficients[order]), norm) * self.TheoryPrefactor(settings.SubLeadTheoryOrder[order], norm)

        pred = np.matmul(self.exp_data[key]['Smear'],pred) 

        return pred
        
    
    # Prefactor for leading and subleading theory
    def TheoryPrefactor(self, TheoryOrder:str, norm:float):

        value = settings.N0

        if ('22' in TheoryOrder):
            value *= settings.C2C2
            value /= norm
            value *= np.power(settings.VtbVts * self.mb, 2.0)
        elif ('SSF27' in TheoryOrder):
            value *= settings.C2C7 / np.sqrt(norm) * settings.VtbVts * self.mb
            value *= settings.La2 / self.mb
        elif ('27' in TheoryOrder):
            value *= settings.C2C7 / np.sqrt(norm) * settings.VtbVts * self.mb
        elif ('28' in TheoryOrder):
            value *= settings.C2C8 / norm * np.power(settings.VtbVts * self.mb, 2.0)
        elif ('88' in TheoryOrder):
            value *= settings.C8C8
            value /= norm
            value *= np.power(settings.VtbVts * self.mb, 2.0)
        elif ('78'  in TheoryOrder):
            value *= settings.C8C7
            value /= np.sqrt(norm)
            value *= settings.VtbVts * self.mb
        
        return value


    # Returns a prediction array for one specific experiment with one specific lambda
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


    # mid is here the SSF27 Files --> SubLeadingTheoryOrder, called in BsgPrediction_full()
    def BsgSubLeadingPrediction(self, key:str, c_n_params, norm):
        pred = np.array([])
        for i in range(np.size(self.exp_data[key]['dBFs'])):
            value = 0
            for j in range(np.size(c_n_params)):
                try:
                    value +=  self.theory_SSF[key]['la'+ settings.SubLeadBasisExpansion]['Values'][i][0][j] * c_n_params[j]
                    value *= 1/(self.theory_SSF[key]['la'+ settings.SubLeadBasisExpansion]['lambda']*norm)
                except IndexError:
                    #print('SubLeadingTheory for %s and %s has only %d values per bin, so we just used the first %d Parameters' % (key,settings.SubLeadBasisExpansion,j, j))
                    break
            pred = np.append(pred, value)

        return pred
    
    def SubLeadPars_read(self, key:str, end:str): #TODO: eigentlich Quatsch diese Funktion
        return self.theory_SSF[key]['SSF27'][end]['Values']
    
    # puts Zoltan's SF into the prediction
    # converts c_n's into subleading coefficients

    #TODO: c_n_params not used, also in c++ code?
    def SubLeadPars(self, c_n_params, d2: float):
        
        #TODO: Find the real values and how/when/where to read them in EDIT: Found them in fit.config, are they correct? Is there another file
        Rho2 = -0.05
        mB = 5.279
        mb = self.mb1SPrediction(c_n_params)
        la = self.lambda11SPrediction(c_n_params)
        Lambda2 = 0.12 


        if('SSF27_1' in settings.SubLeadTheoryOrder):
            x = (0.6514810199386504 * (mB - mb) - 0.8686413599182005 * la + 0.3257405099693252 * Rho2 / Lambda2) / (1.8531015678254943 * la - d2 * la)
        elif('SSF27_2' in settings.SubLeadTheoryOrder):
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
                value += self.Fmn_moments['expx3']['Values'][order][i][j] * c_n_params[i] * c_n_params[j] #TODO: Dont know if indices are correct
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
        u = np.cbrt((-q + np.sqrt(q * q + p * p * p)))
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
    
    #converts an's into cn's
    def ConvertPars(pars):
        length = np.size(pars)-1
        if(length < 0):
            length = np.size(settings.FitVars)-1
        cn = np.array([])
        annorm = 1.0
        for i in range(length):
            annorm += pars[i+1]**2
        np.append(cn, 1/np.sqrt(annorm))
        for i in range(length):
            np.append(cn, pars[i+1]/np.sqrt(annorm))

        return cn

    #The given parameters pars are probaply the an's which should be converted into cn's
    def Chisq(self, pars):
        cn = Meas.ConvertPars(pars)
        norm = pars[0]

        pred_glob = np.array([])
        meas_glob = np.array([])
        Cov_glob = np.zeros((52,52)) #FIXME Oddly specific

        ntot = 0

        for key in settings.KeyOrder:

            nbins = np.size(self.exp_data[key]['dBFs'])
            min = self.fit_config[key]['min']
            max = nbins # TODO:In the C++ Code, there is an if else statement, don't really know why we need it

            full_pred = Meas.BsgPrediction_full(self,key, settings.BasisExpansion, pars, norm) # My BsgPrediction Function is for one data set
            pred =  full_pred[min:max+1]
            meas = self.exp_data[key]['dBFs'][min:max+1]
            Cov = self.exp_data[key]['Cov'][min:max+1,min:max+1]
            pred_glob = np.append(pred_glob, pred)
            meas_glob = np.append(meas_glob, meas)
            Tools.set_sub(Cov_glob, ntot, ntot, Cov)
            ntot += max-min
        
        Cov_glob = np.linalg.inv(Cov_glob) 


        #Calculate Chi^2
        Chisq = np.dot(meas_glob-pred_glob, np.matmul(Cov_glob, meas_glob-pred_glob)) # NOTE: Correct like this ?

        self.chisq = Chisq
        self.mb = self.mb1SPrediction(cn)
        self.M1 = self.Moment(cn,1)
        self.M2 = self.Moment(cn,2)
        self.M3 = self.Moment(cn,3)
        self.la = settings.BasisExpansion # TODO: in C++ via a function, which returns a string with the used _expansion, I guess it is my 'end'
        return Chisq


'''
start_pars = np.array([1, 0.00506919, 0.0, 0.0798100, 0.0870341, 0.0250290, 0.0])

mes = Meas()

m = Minuit(mes.Chisq, start_pars)
m.migrad()
m.hesse()


with open('results.csv', 'a') as f:
    #f.write('la;x0;x1;x2;x3;x4;x5;x6\n')
    f.write('%s;%f;%f;%f;%f;%f;%f;%f\n' % (settings.BasisExpansion,m.values[0],m.values[1],m.values[2],m.values[3],m.values[4],m.values[5],m.values[6]))

'''


#print(Meas.SubLeadPars('babar_hadtag', '105')['Values'][0])

#print(h.BsgPrediction('babar_hadtag','NNLLNNLO', '055', test_fit_results, test_norm))

#print(h.theory_expx3['babar_hadtag']['NNLLNNLO']['la'+'055']['Values'][0][0][0]) #NOTE: ...[#bin][line][column]
