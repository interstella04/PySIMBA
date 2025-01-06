import numpy as np
import pickle
import pandas as pd
from Tools import Tools


class Meas:
    def __init__(self):

        #Theory Dictionaries
        self.theory_expx3 = Tools.pickle_to_dict('theory/theory_dictionary_expx3')
        self.theory_SSF = Tools.pickle_to_dict('theory/theory_dictionary_expx3')

        #Experimental Dictionary
        #TODO put it in one file to open
        self.exp_data = {
            "babar_incl": Tools.pickle_to_dict("data/babar_incl"),
            "babar_hadtag": Tools.pickle_to_dict("data/babar_hadtag"),
            "babar_sem": Tools.pickle_to_dict("data/babar_sem"),
            "belle": Tools.pickle_to_dict("data/belle")
            }
        
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
    
    def SubLeadPars(self, key, end):
        return self.theory_SSF[key]['SSF27'][end]['Values']



    

#h = Meas()   
    

#print(h.exp_data['belle']['Smear'])

#test_fit_results = np.array([0.9956, 0.0641, 0.0624, 0.0267])
#test_norm = 4.925



#print(Meas.SubLeadPars('babar_hadtag', '105')['Values'][0])

#print(h.BsgPrediction('babar_hadtag','NNLLNNLO', '055', test_fit_results, test_norm))

#print(h.theory_expx3['babar_hadtag']['NNLLNNLO']['la'+'055']['Values'][0][0][0]) #NOTE: ...[#bin][line][column]