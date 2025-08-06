import numpy as np
import matplotlib.pyplot as plt 
from Meas import Theory
from Meas import settings
from matplotlib import rc
from Tools import Tools


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 13})

class result:



    def __init__(self):
        #Dictionary with experimental data
        self.ExpData = Tools.PickleToDict(settings.config["MeasurementPath"])
        
        #Histogram Dictionary from Root Data
        self.HistNomFit = Tools.PickleToDict("theory/hist_nom")

        pass

    def PrepareExpData(self, key: str, div_bin: bool = False):

        # Calculate the middle of the bin
        # x error (dx) is bin width
        x = np.zeros(np.size(self.ExpData[key]['dBFs']))
        dx = np.zeros(np.size(self.ExpData[key]['dBFs']))
        for i in range(np.size(self.ExpData[key]['dBFs'])):
            dx[i] = (self.ExpData[key]['Bins'][i+1] - self.ExpData[key]['Bins'][i])/2
            x[i] = dx[i]+self.ExpData[key]['Bins'][i]
        
        # Divide by the width of the bin or not
        if(div_bin == False):
            y = self.ExpData[key]['dBFs']
            dy = self.ExpData[key]['dBFErrors']
        else:
            y = self.ExpData[key]['dBFs']/(2*dx)
            dy = self.ExpData[key]['dBFErrors']/(2*dx)

        ''' Catch Scale for Measurements '''
        #if key == 'belle': #TODO: BELLE SCALE WHAAAT?
        #    ax.errorbar(x, y, dy, dx, fmt = 'k.', markersize = 5, label = '$\\mathrm{Exp. Data}$')
        #    ax.set_ylabel('$\\mathrm{Events} \\: [10^{3}  / 50 \\mathrm{MeV}]$', fontsize = 16)
        
        # return the Label of the Measurement to put in the title
        label = self.ExpData[key]['Label']

        return x, y, dx, dy, label 


    def SimplePlot(ax, x, y, dx = None, dy = None, label: str = None , box_opt: bool = False):

        # Set Box options
        if (box_opt):
            textstr = '\n'.join((
                '$\\mathrm{Entries:} \\: %d$' % (np.size(x),),
                '$\\mathrm{Mean:} \\: %0.3f$' % (x.mean(),),
                '$\\mathrm{Std Dev:} \\: %0.4f$' % (x.std(),)))
            ax.text(0.8,0.95,textstr,fontsize=12,color="black",bbox=dict(facecolor="white", alpha=1, edgecolor="black"), transform = ax.transAxes)
        
        
        ax.errorbar(x, y, dy, dx, fmt = 'k.', markersize = 5, label = '$\\mathrm{Exp. Data}$')
        ax.set_ylabel('$\\Delta B \\left( B \\rightarrow X_S \\gamma \\right) \\: [10^{-4}  / 0.1 \\mathrm{GeV}]$', fontsize = 16)
              
        ax.set_xlabel('$E[\\mathrm{GeV}]$', fontsize = 16, loc = "right")
        ax.set_title('$\\mathrm{%s}$' % label, fontsize = 20)

        return
    
    #def CalculatePrediction():


    #def calc_pred(self, key, an):
    #    cn = Meas.ConvertPars(an)
    #    norm = an[0]
    #    pred  = self.measurement.FullBsgPrediction(key,settings.BasisExpansion, cn, norm)
    #    #if key == 'belle': pred*= 1/1000

    #    return pred
    
    def Plot(self, key, div_bin = False, box_opt = False):
        fig, ax = plt.subplots(figsize=(8, 6))

        x, y, dx, dy, label = result.PrepareExpData(self, key, div_bin)
        result.SimplePlot(ax, x, y, dx, dy, label, box_opt)
        fig.savefig('data/'+key, dpi = 500)

        return
    


res = result()  
res.Plot('babar_hadtag', False, False)