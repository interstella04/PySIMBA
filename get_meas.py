from ROOT import *
import numpy as np
import pickle
from Meas import Meas
 
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

# We can change the titles of the Plots here
data_label_list = ["Babar Inclusive Spectra", "Babar Hadronic Tag", "Babar Semileptonic", "Belle"]

data_tag_list = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]


for i in range(4):
    Meas.store_in_pickle(GrabMeasurement(data_tag_list[i],data_label_list[i]), data_tag_list[i])
