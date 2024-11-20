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

        end_NLL = ['03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                    '035', '045', '055', '065', '075', '085', '095', '0475',
                    '0525', '0575' , '0625']
        
        mids = ['NNLLNNLO', 'NS22NNLO', 'NS27NNLO', 'NS28NNLO', 'NS78NNLO', 'NS88NNLO']

        if mid == mids[0]:
            start_line = 17
        elif mid == mids[5] or mid == mids[4]:
            start_line = 12
        else:
            start_line = 13

        if key == 'babar_sem':
            strip_string = '{fmX2, }'
        else:
            strip_string = '{fEgY, }'

        collected_data = {}

        for j, end in enumerate(end_NLL):
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
                "name": 'la'+end,
                "bins": bins,
                "width": width
            }

            collected_data['la'+end] = theo_dictio

        return collected_data
    
    def grab_mids(key):
        mids = ['NNLLNNLO', 'NS22NNLO', 'NS27NNLO', 'NS28NNLO', 'NS78NNLO', 'NS88NNLO']
        final_data = {}

        for i,mid in enumerate(mids):
            print('%s' % mid)
            final_data['%s' % (mid,)] = Meas.grab_theory(key, mid)
        return final_data


theory_dictionary = {
    "babar_hadtag_theo": Meas.grab_mids('babar_hadtag'),
    "babar_incl_theo": Meas.grab_mids('babar_incl'),
    "babar_sem_theo": Meas.grab_mids('babar_sem'),
    "belle_theo": Meas.grab_mids('belle')
    }

Meas.store_in_pickle(theory_dictionary, 'theory/theory_dictionary_expx3')


#d = Meas.grab_theory('belle', 'NNLLNNLO')
#print(d['la03']['width'])

'''
s = '# bin = 12 {fEgY, 2., 2.05}'
s = s.strip('# bin = 12')
print(s)
s = s.strip('{fEgY, }')
print(s)
'''
'''
for i, arr in enumerate(bins):
    print(f"Bin = {i} in range {width[i]}:\n{arr}\n")
'''