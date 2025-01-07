import numpy as np
import pickle
import pandas as pd
from Tools import Tools
import traceback


class Meas:
    def __init__(self):

        #Theory Dictionaries
        self.theory_expx3 = Tools.pickle_to_dict('theory/theory_dictionary_expx3')
        self.theory_SSF = Tools.pickle_to_dict('theory/theory_dictionary_expx3')

        #Dictionary with experimental data
        self.exp_data = Tools.pickle_to_dict('data/exp_data')
        
        #Histogram Dictionary from Root Data
        self.hist_nom = Tools.pickle_to_dict("theory/hist_nom")

        return

    # Returns a prediction array for one specific experiment with one specific lambda
    def BsgPrediction(self, key, mid, end, c_n_params, norm):      
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
    def BsgSubLeadingPrediction(self, key, mid, end, c_n_params, norm):
        pred = np.array([])
        for i in range(np.size(self.exp_data[key]['dBFs'])):
            value = 0
            for j in range(np.size(c_n_params)):
                value +=  self.theory_expx3[key][mid]['la'+end]['Values'][i][0][j] * c_n_params[j]
                value *= 1/(self.theory_expx3[key][mid]['la'+end]['lambda']*norm)
            pred = np.append(pred, value)

        return pred
    
    def SubLeadPars_read(self, key, end):
        return self.theory_SSF[key]['SSF27'][end]['Values']
    
    # dds Zoltan's SF into the prediction
    # converts c_n's into subleading coefficients 
    def SubLeadPars(c_n_params, d2, la, opt):
        
        #TODO: Find the real values and how/when to read them in
        Rho2 = 0
        mB = 5.279
        mb = 1
        Lambda2 = 0.105 #Found in lam2rho2_lain.pdf in GeV^2

        if(opt == 1):
            x = (0.6514810199386504 * (mB - mb) - 0.8686413599182005 * la + 0.3257405099693252 * Rho2 / Lambda2) / (1.8531015678254943 * la - d2 * la)
        elif(opt == 2):
            x = (0.4722982832332954 * (mB - mb) - 0.5667579398799545 * la + 0.2361491416166477 * Rho2 / Lambda2) / (1.3695121740357048 * la - d2 * la)
        else:
            print('Unknown subleading shape function -- please specify Meas.SubLeadPars()')

        #calculating d0, d1 and d2 and returning them in a numpy array
        return np.array([1.-x, x* (1-d2), x*d2])



#print(Meas.SubLeadPars('babar_hadtag', '105')['Values'][0])

#print(h.BsgPrediction('babar_hadtag','NNLLNNLO', '055', test_fit_results, test_norm))

#print(h.theory_expx3['babar_hadtag']['NNLLNNLO']['la'+'055']['Values'][0][0][0]) #NOTE: ...[#bin][line][column]