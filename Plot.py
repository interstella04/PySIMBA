import numpy as np
import pickle
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import pandas as pd
from pathlib import Path
from Meas import Meas
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 13})

class Plot:
    
    def __init__(self):
        self.dictio_of_dictios = {
            "Babar_incl": Meas.pickle_to_dict("babar_incl"),
            "Babar_hadtag": Meas.pickle_to_dict("babar_hadtag"),
            "Babar_sem": Meas.pickle_to_dict("babar_sem"),
            "Belle": Meas.pickle_to_dict("belle")
            }
        return
        
    def simple_plot(self, key, div_bin = False, box_opt = False):

        fig, ax = plt.subplots()

        x = np.zeros(np.size(self.dictio_of_dictios[key]['dBFs']))
        dx = np.zeros(np.size(self.dictio_of_dictios[key]['dBFs']))
        
        for i in range(np.size(self.dictio_of_dictios[key]['dBFs'])):
            dx[i] = (self.dictio_of_dictios[key]['Bins'][i+1] - self.dictio_of_dictios[key]['Bins'][i])/2
            x[i] = dx[i]+self.dictio_of_dictios[key]['Bins'][i]
        
        # Divide by the width of the bin or not
        if(div_bin == False):
            y = self.dictio_of_dictios[key]['dBFs']
            dy = self.dictio_of_dictios[key]['dBFErrors']
        else:
            y = self.dictio_of_dictios[key]['dBFs']/dx
            dy = self.dictio_of_dictios[key]['dBFErrors']/dx

        #Set Box options
        if (box_opt):
            textstr = '\n'.join((
                '$\\mathrm{Entries:} \\: %d$' % (np.size(x),),
                '$\\mathrm{Mean:} \\: %0.3f$' % (x.mean(),),
                '$\\mathrm{Std Dev:} \\: %0.4f$' % (x.std(),)))
            ax.text(0.8,0.95,textstr,fontsize=12,color="black",bbox=dict(facecolor="white", alpha=1, edgecolor="black"), transform = ax.transAxes)
        
        if key == 'Belle':
            ax.errorbar(x, y*10**(-3), dy*10**(-3),dx,fmt = 'k.', markersize = 5)
            ax.set_ylabel('$\\mathrm{Events} \\: [10^{3}  / 50 \\mathrm{MeV}]$', fontsize = 16)
        else:
            ax.errorbar(x, y,dy,dx, fmt = 'k.', markersize = 5)
            ax.set_ylabel('$\\Delta B \\left( B \\rightarrow X_S \\gamma \\right) \\: [10^{-4}  / 0.1 \\mathrm{GeV}]$', fontsize = 16)
              
        ax.set_xlabel('$E[\\mathrm{GeV}]$', fontsize = 16, loc = "right")
        ax.set_title('$\\mathrm{%s}$' % self.dictio_of_dictios[key]["Label"],fontsize = 20)
        plt.savefig('data/'+key, dpi = 500)

h1 = Plot()
h2 = Plot()
h3 = Plot()
h4 = Plot()
h1.simple_plot('Babar_hadtag', box_opt= True)
h2.simple_plot('Belle', box_opt=True)
h3.simple_plot('Babar_sem', box_opt= True)
h4.simple_plot('Babar_incl', box_opt= True)