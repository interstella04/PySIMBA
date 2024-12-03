import numpy as np
import matplotlib.pyplot as plt 
from Meas import Meas
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 13})

class Plot:
    
    def __init__(self):
        self.measurement = Meas()
        return
        
    def simple_plot(self, key, div_bin = False, box_opt = False):

        self.measurement.exp_data
        fig, ax = plt.subplots()

        x = np.zeros(np.size(self.measurement.exp_data[key]['dBFs']))
        dx = np.zeros(np.size(self.measurement.exp_data[key]['dBFs']))
        
        for i in range(np.size(self.measurement.exp_data[key]['dBFs'])):
            dx[i] = (self.measurement.exp_data[key]['Bins'][i+1] - self.measurement.exp_data[key]['Bins'][i])/2
            x[i] = dx[i]+self.measurement.exp_data[key]['Bins'][i]
        
        # Divide by the width of the bin or not
        if(div_bin == False):
            y = self.measurement.exp_data[key]['dBFs']
            dy = self.measurement.exp_data[key]['dBFErrors']
        else:
            y = self.measurement.exp_data[key]['dBFs']/dx
            dy = self.measurement.exp_data[key]['dBFErrors']/dx

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
        ax.set_title('$\\mathrm{%s}$' % self.measurement.exp_data[key]["Label"],fontsize = 20)
        return x

    #NOTE: Could also be an if-statement in simple Plot
    def plot_exp(self, key, div_bin = False, box_opt = False):
        Plot.simple_plot(self, key, div_bin, box_opt)
        plt.savefig('data/'+key, dpi = 500)
        return

    def check_pred(self, key, mid, end, div_bin = False, box_opt = False):
        x = Plot.simple_plot(self, key, div_bin, box_opt)
        test_fit_results = np.array([0.9956, 0.0641, 0.0624, 0.0267])
        test_norm = 1.6#4.925 #FIXME: Correct normalization? Not working for lambda with c_n arrays smaller than the data in the files

        h = Meas()
        y = h.BsgPrediction(key,mid,end, test_fit_results, test_norm)
        plt.errorbar(x,y, fmt= 'r.', markersize = 5)
        plt.savefig('test_'+key)
         



h = Plot()
h.check_pred('babar_hadtag','NNLLNNLO', '055', box_opt= True)
h.check_pred('babar_sem','NNLLNNLO', '055', box_opt= True)
h.check_pred('babar_incl','NNLLNNLO', '055', box_opt= True)

'''
h1 = Plot()
h2 = Plot()
h3 = Plot()
h4 = Plot()
h1.plot_exp('babar_hadtag', box_opt= True)
h2.plot_exp('belle', box_opt=True)
h3.plot_exp('babar_sem', box_opt= True)
h4.plot_exp('babar_incl', box_opt= True)
'''