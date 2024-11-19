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
    def grab_theory(key, mid):

        if mid == 'NNLLNNLO':
            end_NLL = ['03', '07']
            '''['03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                    '035', '045', '055', '065', '075', '085', '095', '0475',
                    '0525', '0625']'''
        else:
            print('wurde noch nicht implementiert ^^ ')
            return
        
        collected_data = {}

        for j, end in enumerate(end_NLL):
            count = 0
            curr_block = []
            bins = []
            width = []
            with open('theory/'+ key + '_'+mid+'_la'+end+'.txt', 'r') as file:
                for i, line in enumerate(file):
                    if i < 17:  # jumps to line 18, only correct for NNLLNNLO Data
                        continue
                    line = line.strip()  # removes empty lines
                    if line.startswith("# bin"):  # skips line starting with # bin
                        line = line.strip('# bin = %d {fEgY, }' % (count,))
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
                "name": 'la'+end,
                "bins": bins,
                "width": width
            }

            collected_data['la'+end] = theo_dictio

        return collected_data

d = Meas.grab_theory('babar_hadtag', 'NNLLNNLO')

print(d['la07']['width'])
'''
for i, arr in enumerate(bins):
    print(f"Bin = {i} in range {width[i]}:\n{arr}\n")
'''