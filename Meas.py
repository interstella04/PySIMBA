import numpy as np
import pickle
import pandas as pd

class Meas:

    def store_in_pickle(dictio):
        with open(dictio["Label"]+".pkl", "wb") as file:
            pickle.dump(dictio,file)
        return
            
    def pickle_to_dict(pkl_file):
        with open("data/"+pkl_file + ".pkl", "rb") as file:
            data = pickle.load(file)
        return data
    
    '''Function to grab theory, currently just for one NNLLNNLO data'''
    def grab_theory(key):
        '''end_nll = ['03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                    '035', '045', '055', '065', '075', '085', '095', '0475',
                    '0525', '0625']'''
        bins = []
        curr_block = []
        #for j in range(np.size(end_nll)):
        with open('theory/'+ key + '_NNLLNNLO_la07.txt', 'r') as file:
            for i, zeile in enumerate(file):
                if i < 18:  # jumps to line 19, only correct for NNLLNNLO Data
                    continue
                zeile = zeile.strip()  # removes empty lines
                if zeile.startswith("# bin"):  # skips line starting with # bin
                    continue
                if zeile == "":
                    if curr_block:  # if current block not empty
                        bins.append(np.array(curr_block, dtype=float))
                        curr_block = []
                else:
                    # extract numbers out of line and put them into an array
                    curr_block.append(list(map(float, zeile.split())))

        if curr_block:
            bins.append(np.array(curr_block, dtype=float))
        
        return bins
        
        
'''
bins = Meas.grab_theory('babar_hadtag')

for i, arr in enumerate(bins):
    print(f"Bin = {i}:\n{arr}\n")
'''