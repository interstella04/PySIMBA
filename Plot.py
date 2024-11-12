import numpy as np
import pickle
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import pandas as pd
from pathlib import Path
from Meas import Meas

class Plot:
    
    def __init__(self):
        self.dictio_of_dictios = {
            "Babar_incl": Meas.pickle_to_dict("babar_incl"),
            "Babar_hadtag": Meas.pickle_to_dict("babar_hadtag"),
            "Babar_sem": Meas.pickle_to_dict("babar_sem"),
            "Belle": Meas.pickle_to_dict("belle")
            }
        return
        
    def simple_plot(self, key):
        
        x = np.zeros(np.size(self.dictio_of_dictios[key]['dBFs']))
        dx = np.zeros(np.size(self.dictio_of_dictios[key]['dBFs']))
        
        for i in range(np.size(self.dictio_of_dictios[key]['dBFs'])):
            dx[i] = (self.dictio_of_dictios[key]['Bins'][i+1] - self.dictio_of_dictios[key]['Bins'][i])/2
            x[i] = dx[i]+self.dictio_of_dictios[key]['Bins'][i]
        
        mean = x.mean()
        dev = x.std()
        entries = np.size(x)
        #Textbox settings not the right vakues
        textstr = '\n'.join((
            r'Entries: %d' % (entries,),
            r'Mean: %0.3f' % (mean,),
            r'Std Dev: %0.4f' % (dev,)))
        
        #Font settings (LATEX NOT FOUND)
        #plt.rcParams.update({
                #"text.usetex": True,
                #"font.family": "sans-serif",
                #"font.sans-serif": "Helvetica"})
            
        if key == 'Belle':
            plt.errorbar(x, self.dictio_of_dictios[key]['dBFs']*10**(-3), self.dictio_of_dictios[key]['dBFErrors']*10**(-3),dx,fmt = 'k.', markersize = 5)
            plt.ylabel(r'Events $[10^{3}  / 50 \mathrm{MeV}]$', fontsize = 16)
        else:
            plt.errorbar(x, self.dictio_of_dictios[key]['dBFs'],self.dictio_of_dictios[key]['dBFErrors'],dx, fmt = 'k.', markersize = 5)
            plt.ylabel(r'$\Delta B \left( B \rightarrow X_S \gamma \right) [10^{-4}$  / 0.1 $\mathrm{GeV}]$', fontsize = 16)
            
        plt.xlabel(r'E[GeV]', fontsize = 16, loc = "right")
        plt.title(r''+ self.dictio_of_dictios[key]["Label"],fontsize = 20)
        #plt.text(0.05, 0.95, textstr, fontsize = 14, verticalalignment = 'top') not working
        plt.savefig(''+self.dictio_of_dictios[key]["Label"], dpi = 500)

h = Plot()
h.simple_plot('Belle')
