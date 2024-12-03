import numpy as np
import pickle
import pandas as pd


class Meas:
    def __init__(self):
        #Lambdas and endings of the filenames
        self.end_lambda = ['03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                    '035', '045', '055', '065', '075', '085', '095', '0475',
                    '0525', '0575' , '0625']
        # Possible Middle Terms of the filenames
        self.mids = ['NNLLNNLO', 'NS22NNLO', 'NS27NNLO', 'NS28NNLO', 'NS78NNLO', 'NS88NNLO']
        
        self.theory_expx3 = Meas.pickle_to_dict('theory/theory_dictionary_expx3')

        self.exp_data = {
            "babar_incl": Meas.pickle_to_dict("data/babar_incl"),
            "babar_hadtag": Meas.pickle_to_dict("data/babar_hadtag"),
            "babar_sem": Meas.pickle_to_dict("data/babar_sem"),
            "belle": Meas.pickle_to_dict("data/belle")
            }
        return
    
    # These two following functions are tools
    def store_in_pickle(dictionary, tag):
        with open(tag+".pkl", "wb") as file:
            pickle.dump(dictionary,file)
        return
            
    def pickle_to_dict(pkl_file):
        with open(pkl_file + ".pkl", "rb") as file:
            data = pickle.load(file)
        return data
    
    def grab_theory(self, key, mid):
        
        # Different start lines of the data
        if mid == self.mids[0]:
            start_line = 17
        elif mid == self.mids[5] or mid == self.mids[4]:
            start_line = 12
        else:
            start_line = 13
        
        # Different Strings to be deletet in the Data
        if key == 'babar_sem':
            strip_string = '{fmX2, }'
        else:
            strip_string = '{fEgY, }'

        # Dictionary with all data for a specific experiment and a specific mid 
        collected_data = {}

        #Reading out the data, currently works only with a specific structur of folders
        for j, end in enumerate(self.end_lambda):
            count = 0
            curr_block = []
            bins = []
            width = []
            with open('../simba/Cpp/share/simba/theory/mb47_mc13_nf3_as207_expx3/'+ key + '_'+ mid +'_la'+ end +'.txt', 'r') as file:
                for i, line in enumerate(file):
                    if i < start_line:  # jumps to line 18, only correct for NNLLNNLO Data
                        continue
                    line = line.strip()  # removes empty lines
                    if line.startswith("# bin ") :  # skips line starting with # bin
                        line = line.strip('# bin = %d' % (count,))
                        line = line.strip('%s' % strip_string)
                        width.append(list(map(float, line.split(','))))
                        count+=1
                        continue
                    if line == "":
                        if curr_block:  # if current block not empty
                            bins.append(np.array(curr_block, dtype=float))
                            curr_block = []
                    else:
                        # extract numbers out of line and put them into an array
                        curr_block.append(list(map(float, line.split())))

            if curr_block:
                bins.append(np.array(curr_block, dtype=float))
        
            theo_dictio = {
                "lambda": int(end) / (10 ** (len(end) - 1)),
                "bins": bins,
                "width": width
            }

            collected_data['la'+end] = theo_dictio

        return collected_data
    
    # Reads out every file of a given experiment name and returns them in a dictionary
    def grab_mids(self, key):
        final_data = {}

        for i,mid in enumerate(self.mids):
            final_data['%s' % (mid,)] = self.grab_theory(key, mid)
        return final_data
    
    
    # Returns a prediction array for one specific experiment with one specific lambda
    def BsgPrediction(self, key, mid, end, c_n_params, norm):      
        val_arr = np.array([])
        for i in range(np.size(self.exp_data[key]['dBFs'])):
            value = 0
            for j in range(np.size(c_n_params)):
                for k in range(np.size(c_n_params)):
                    value += self.theory_expx3[key][mid]['la'+end]['bins'][i][j][k] * c_n_params[j] * c_n_params[k] 
                    value *= 1/(self.theory_expx3[key][mid]['la'+end]['lambda']*norm)
            val_arr = np.append(val_arr, value)
        return val_arr 




#test_fit_results = np.array([0.9956, 0.0641, 0.0624, 0.0267])
#test_norm = 4.925


#h = Meas()
#print(h.BsgPrediction('babar_hadtag','NNLLNNLO', '055', test_fit_results, test_norm))


#print(h.theory_expx3['babar_hadtag']['NNLLNNLO']['la'+'055']['bins'][0][0][0]) #NOTE: [bin][line][column]

# Theory Dictionary in pickle with all data of expx3 
'''
theory_dictionary = {
    "babar_hadtag": h.grab_mids('babar_hadtag'),
    "babar_incl": h.grab_mids('babar_incl'),
    "babar_sem": h.grab_mids('babar_sem'),
    "belle": h.grab_mids('belle')
    }

Meas.store_in_pickle(theory_dictionary, 'theory/theory_dictionary_expx3')

print(theory_dictionary['babar_hadtag']['NNLLNNLO']['la03']['lambda'])
'''