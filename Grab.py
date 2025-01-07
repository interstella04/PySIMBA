#from ROOT import *
import uproot
import numpy as np
import pickle
from Tools import Tools

class Grab:

    def __init__(self):
        #Lambdas and endings of the filenames for expx folder and SSF folder
        self.end_lambda = ['03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                    '035', '045', '055', '065', '075', '085', '095', '0475',
                    '0525', '0575' , '0625']
        
        self.end_SSF = ['105', '105', '107', '205', '206', '207', '1045', '1055', '1065', '2045',
                         '2055', '2065', '10525', '10575', '10625', '20525', '20575']
        # Possible Middle Terms of the filenames
        self.mids = ['NNLLNNLO', 'NS22NNLO', 'NS27NNLO', 'NS28NNLO', 'NS78NNLO', 'NS88NNLO']


    '''
    #TODO: Compile again because of the label names
    def GrabMeasurement(Tag, Label):

        # Our Measurement Dictionary
        measDict = {
            "Label": Label,
            "Bins" : np.array([]),
            "dBFs" : np.array([]),
            "dBFErrors" : np.array([]),
            "Cov" : np.array([[]]),
            "Cor" : np.array([[]]),
            "Smear" : np.array([[]])
        }

        # Open File and grab measurement histogram
        f = TFile("/home/s08shoff/physik/shk_bernlochner/simba/Cpp/share/simba/measurements/" + Tag + ".root")
        h = f.Get(Tag)
        # Determine the number of bins
        nbins = h.GetNbinsX()
        # Grab the measurement covariance and correlations
        cov = f.Get("Cov")
        cor = f.Get("Corr")
        smear = f.Get("Smear")

        # Loop over all bins
        for i in range(1,nbins+1):

            # Set measured differential / partial branching fractions of a given bin and their errors
            measDict["dBFs"] = np.append( measDict["dBFs"], h.GetBinContent(i))
            measDict["dBFErrors"] = np.append( measDict["dBFErrors"], h.GetBinError(i))
            measDict["Bins"] = np.append( measDict["Bins"], h.GetBinLowEdge(i))

            # Set Covariance matrix of the measured bins
            c = np.array([])
            d = np.array([]) 
            s = np.array([])
            for j in range(1,nbins+1):
                c = np.append( c, cov[i-1][j-1] )
                d = np.append( d, cor[i-1][j-1] )
                s = np.append( s, smear[i-1][j-1] )

            # If its the first row, we just set it equial, otherwise we stack the items
            measDict["Cov"] = c if i == 1 else np.vstack((measDict["Cov"], c))
            measDict["Cor"] = d if i == 1 else np.vstack((measDict["Cor"], d))
            measDict["Smear"] = s if i == 1 else np.vstack((measDict["Smear"], s))

        # Add the last bin (as we have one more bin boundary than # of measurements)
        measDict["Bins"] = np.append( measDict["Bins"], h.GetBinLowEdge(nbins+1))    

        return measDict

    '''

    #key: Which experiment. possibilities: ["babar_incl", "babar_hadtag", "babar_sem", "belle"]
    #mid: Which Order. possibilities: see self.mids
    #
    def grab_theory(self, key, mid):
        
        # Different start lines of the data
        if mid == self.mids[0]:
            start_line = 17
            end_array = self.end_lambda
        elif mid == self.mids[5] or mid == self.mids[4]:
            start_line = 12
            end_array = self.end_lambda
        elif mid == 'SSF27':
            start_line = 12
            end_array = self.end_SSF
        else:
            start_line = 13
            end_array = self.end_lambda
        
        # Different Strings to be deletet in the Data
        if key == 'babar_sem':
            strip_string = '{fmX2, }'
        else:
            strip_string = '{fEgY, }'

        # Dictionary with all data for a specific experiment and a specific mid 
        collected_data = {}

        #Reading out the data, currently works only with a specific structur of folders
        for j, end in enumerate(end_array):
            count = 0
            curr_block = []
            values = []
            bins = []
            if mid == 'SSF27':
                file = open('../simba/Cpp/share/simba/theory/mb47_mc13_nf3_as207_SSF27/'+key+'_SSF27_'+end+'.txt', 'r')
            else:
                file = open('../simba/Cpp/share/simba/theory/mb47_mc13_nf3_as207_expx3/'+ key + '_'+ mid +'_la'+ end +'.txt', 'r')
            for i, line in enumerate(file):
                if i < start_line:  # jumps to start line, where the data begins
                    continue
                line = line.strip()  # removes useless lines
                if line.startswith("# bin ") :  # skips line starting with # bin
                    line = line.strip('# bin = %d' % (count,))
                    line = line.strip('%s' % strip_string)
                    bins.append(list(map(float, line.split(','))))
                    count+=1
                    continue
                if line == "":
                    if curr_block:  # if current block not empty
                        values.append(np.array(curr_block, dtype=float))
                        curr_block = []
                else:
                    # extract numbers out of line and put them into an array
                    curr_block.append(list(map(float, line.split())))

            if curr_block:
                values.append(np.array(curr_block, dtype=float))

            file.close()

            if mid == 'SSF27':
                new = end.strip('1')
                new = new.strip('2')
            else:
                new = end

            theo_dictio = {
                    "lambda": int(new) / (10 ** (len(new) - 1)),
                    "Values": values,
                    "Bins": bins
                }

            collected_data['la'+end] = theo_dictio

        return collected_data
    
    
    # Reads out every file of a given experiment name and returns them in a dictionary
    def grab_mids(self, key):
        final_data = {}

        for i,mid in enumerate(self.mids):
            final_data['%s' % (mid,)] = self.grab_theory(key, mid)
        return final_data

    def GrabNomfit(Tag, Label):
        hist_dict = {
            'Label': Label, 
            'Bins': np.array([]),
            'Values': np.array([])
        }
        file = uproot.open('../simba/Cpp/share/simba/measurements/script/nomfit_mom_3001_NNLLNNLO_la06_'+Tag+'_fit.root')
        hist_dict['Bins'] = file['fit'].axis().edges()
        hist_dict['Values'] = file['fit'].values()
        return hist_dict



g = Grab()

#Experimental dictionary from 'measurements/%key%.root'

# We can change the titles of the Plots here
data_label_list = ["Babar\\Inclusive\\Spectra", "Babar\\Hadronic\\Tag", "Babar\\Semileptonic", "Belle"]

data_tag_list = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]

exp_data = {
            "babar_incl": g.GrabMeasurement(data_tag_list[0],data_label_list[0]),
            "babar_hadtag": g.GrabMeasurement(data_tag_list[1],data_label_list[1]),
            "babar_sem": g.GrabMeasurement(data_tag_list[2],data_label_list[2]),
            "belle": g.GrabMeasurement(data_tag_list[3],data_label_list[3])
            }

Tools.store_in_pickle(exp_data, 'data/exp_data')


# Theory Dictionary in pickle with all data of '%key%_NNLLNNLO_la%end%.txt' from mb47_mc13_nf3_as207_expx3

theory_dictionary_expx3 = {
    "babar_hadtag": g.grab_mids('babar_hadtag'),
    "babar_incl": g.grab_mids('babar_incl'),
    "babar_sem": g.grab_mids('babar_sem'),
    "belle": g.grab_mids('belle')
    }

Tools.store_in_pickle(theory_dictionary_expx3, 'theory/theory_dictionary_expx3')

print(theory_dictionary_expx3['babar_hadtag']['NNLLNNLO']['la03']['Bins'])

#Theory Dictionary for SSF27 files (the ones calculated with simba c++ code)

theory_dictionary_SSF27 = {
    "babar_hadtag": g.grab_theory('babar_hadtag', 'SSF27'),
    "babar_incl": g.grab_theory('babar_incl', 'SSF27'),
    "babar_sem": g.grab_theory('babar_sem', 'SSF27'),
    "belle": g.grab_theory('belle', 'SSF27')
    }

Tools.store_in_pickle(theory_dictionary_SSF27, 'theory/theory_dictionary_SSF27')

# Histogram dictionary from root data 'nomfit_mom_3001_NNLLNNLO_la06_%key%.root'
hist_nom = {
    "babar_hadtag": g.GrabNomfit('babar_hadtag'),
    "babar_incl": g.GrabNomfit('babar_incl'),
    "babar_sem": g.GrabNomfit('babar_sem'),
    "belle": g.GrabNomfit('belle')
    }

Tools.store_in_pickle(hist_nom, 'histogram/hist_nom')