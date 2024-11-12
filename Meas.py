from ROOT import *
import numpy as np
import pickle
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import pandas as pd

class Meas:
        
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


    def store_in_pickle(dictio):
        with open(dictio["Label"]+".pkl", "wb") as file:
            pickle.dump(dictio,file)
        return
            
    def pickle_to_dict(pkl_file):
        with open(pkl_file + ".pkl", "rb") as file:
            data = pickle.load(file)
        return data



data_label_list = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]
data_tag_list = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]


for i in range(4):
    Meas.store_in_pickle(Meas.GrabMeasurement(data_tag_list[i],data_label_list[i]))
    #print(data_tag_list[i])
    #print(Meas.GrabMeasurement(data_tag_list[i],data_label_list[i])["dBFs"])
    
