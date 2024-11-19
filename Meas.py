import numpy as np
import pickle
import pandas as pd

class Meas:

    def store_in_pickle(dictio, tag):
        with open(tag+".pkl", "wb") as file:
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
        count = 0
        #for j in range(np.size(end_nll)):
        with open('theory/'+ key + '_NNLLNNLO_la07.txt', 'r') as file:
            for i, line in enumerate(file):
                if i < 18:  # jumps to line 19, only correct for NNLLNNLO Data
                    continue
                line = line.strip()  # removes empty lines
                if line.startswith("# bin"):  # skips line starting with # bin
                    #line = line.strip("# bin = %d {fEg}" % count)
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
        
        return bins
        
string = '# bin = 0 {fEg, 1.9, 2.}'
string = string.strip('# bin = 0 {fEg, }')
print(string)


'''
bins = Meas.grab_theory('babar_hadtag')

for i, arr in enumerate(bins):
    print(f"Bin = {i}:\n{arr}\n")
'''