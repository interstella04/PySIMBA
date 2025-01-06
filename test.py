import uproot
import numpy as np
from Meas import Meas
import matplotlib.pyplot as plt

def GrabNomfit(Tag, Label):
    hist_dict = {
        'Label': Label, 
        'Bins': np.array([]),
        'Values': np.array([]),
        'Errors': np.array([])
    }
    file = uproot.open('../simba/Cpp/share/simba/measurements/script/nomfit_mom_3001_NNLLNNLO_la06_'+Tag+'_fit.root')
    hist_dict['Bins'] = file['fit'].axis().edges()
    hist_dict['Values'] = file['fit'].values()
    hist_dict['Errors'] = file['fit'].errors()
    return hist_dict


data_tag_list = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]
data_label_list = ["Babar Inclusive Spectra", "Babar Hadronic Tag", "Babar Semileptonic", "Belle"]

hist_nom = {}

for i in range(4):
    hist_nom[data_tag_list[i]] = GrabNomfit(data_tag_list[i], data_label_list[i])

Meas.store_in_pickle(hist_nom, 'theory/hist_nom')

print(hist_nom['babar_hadtag']['Bins'])