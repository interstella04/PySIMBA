import numpy as np
import matplotlib.pyplot as plt 
from Meas import Meas
from Meas import settings
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 13})

class Plot:
    
    def __init__(self):
        self.measurement = Meas()
        return
        
    def simple_plot(self, key, fig, ax, div_bin = False, box_opt = False):

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
            y = self.measurement.exp_data[key]['dBFs']/(2*dx)
            dy = self.measurement.exp_data[key]['dBFErrors']/(2*dx)

        #Set Box options
        if (box_opt):
            textstr = '\n'.join((
                '$\\mathrm{Entries:} \\: %d$' % (np.size(x),),
                '$\\mathrm{Mean:} \\: %0.3f$' % (x.mean(),),
                '$\\mathrm{Std Dev:} \\: %0.4f$' % (x.std(),)))
            ax.text(0.8,0.95,textstr,fontsize=12,color="black",bbox=dict(facecolor="white", alpha=1, edgecolor="black"), transform = ax.transAxes)
        
        if key == 'Belle':
            ax.errorbar(x, y, dy,dx,fmt = 'k.', markersize = 5, label = '$\\mathrm{Exp. Data}$')
            ax.set_ylabel('$\\mathrm{Events} \\: [10^{3}  / 50 \\mathrm{MeV}]$', fontsize = 16)
        else:
            ax.errorbar(x, y,dy,dx, fmt = 'k.', markersize = 5, label = '$\\mathrm{Exp. Data}$')
            ax.set_ylabel('$\\Delta B \\left( B \\rightarrow X_S \\gamma \\right) \\: [10^{-4}  / 0.1 \\mathrm{GeV}]$', fontsize = 16)
              
        ax.set_xlabel('$E[\\mathrm{GeV}]$', fontsize = 16, loc = "right")
        ax.set_title('$\\mathrm{%s}$' % self.measurement.exp_data[key]["Label"],fontsize = 20)
        return x,dx


    def plot_histogram(self, bin_edges, values, fig, ax, title = '', xLabel = '', yLabel = '',div_bin = False, box_opt = False):

        #Option to divide the values by the width of the bin
        if(div_bin == True):
            for i in range(np.size(values)):
                values[i] = values[i]/(bin_edges[i+1]-bin_edges[i])

        ax.step(
        bin_edges,          
        [*values, 0],       # Values + additional 0 to close histogram
        where="post",       # Linienverlauf nach rechts
        color="darkgreen",       
        lw=1,
        label = '$\\mathrm{Fit Data}$'                
        )

        
        if (box_opt):
            textstr = '\n'.join((
                '$\\mathrm{Entries:} \\: %d$' % (np.size(values),),
                '$\\mathrm{Mean:} \\: %0.3f$' % (bin_edges.mean(),), #FIXME: Mean is incorrect compared to root data
                '$\\mathrm{Std Dev:} \\: %0.4f$' % (bin_edges.std(),)))
            ax.text(0.8,0.95,textstr,fontsize=12,color="black",bbox=dict(facecolor="white", alpha=1, edgecolor="black"), transform = ax.transAxes)
        
        
        ax.set_xlabel('$\\mathrm{%s}$' % xLabel, fontsize = 16)
        ax.set_ylabel('$\\mathrm{%s}$' % yLabel, fontsize = 16)
        ax.set_title('$\\mathrm{%s}$' % title, fontsize = 20)
        plt.axis([bin_edges[0], bin_edges[-1], None, None])
        plt.grid(axis="y", linestyle="--", alpha=0.6)

        return 
    
    def check_pred(self, key, end, fig, ax, div_bin = False, box_opt = False, given_pred = False, full = True,  mid = 'NNLLNNLO', made_up_norm = False):
        
        x,dx = Plot.simple_plot(self, key, fig, ax, div_bin, box_opt)

        if given_pred == True:
            test_fit_results = np.array([0.9956, 0.0641, 0.0624, 0.0267])
            if made_up_norm == True:
                test_norm = 0.6
            else:
                test_norm = 4.925 * 10 ** (3)
        else:
            test_fit_results = np.array([ 0.9707776, -0.1417625, -0.28887093, 0.88105546, -0.15055346, 0.25913407, -0.13474316])
            test_fit_results = Meas.ConvertPars(test_fit_results)
            if made_up_norm == True:
                test_norm = 0.4 
            else:
                test_norm = 0.97077755
            


        if full == True:
            y = self.measurement.BsgPrediction_full(key,end, test_fit_results, test_norm)
        else:
            y = self.measurement.BsgPrediction(key,mid, end, test_fit_results, test_norm) #+ self.measurement.BsgSubLeadingPrediction(key, test_fit_results, test_norm)
        
        y = np.matmul(self.measurement.exp_data[key]['Smear'],y)

        if div_bin == True:
            ax.errorbar(x,y/(2*dx),fmt='.r', label = '$\\mathrm{Pred. Data}$')
        else:
            ax.errorbar(x,y, fmt = '.r', label = '$\\mathrm{Pred. Data}$')
        return
    

    #
    #There are 3 Options:
    #1: Plots the experimental data
    #2: Plots a histogram from fitted root data
    #3: Plots 1 and 2 in one Plot and the calculated prediction function in one plot
    #
    def plot_exp(self, key, div_bin = False, box_opt = False, plt_opt = 1 ):
        fig, ax = plt.subplots()
        if plt_opt == 1:
            Plot.simple_plot(self, key, fig, ax, div_bin, box_opt)
            plt.savefig('data/'+key, dpi = 500)
        elif plt_opt == 2:
            Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'],self.measurement.hist_nom[key]['Values'], fig, ax,self.measurement.hist_nom[key]['Label'])
            plt.savefig('histogram/'+key+'_nomfit_hist')
        elif plt_opt == 3: #FIXME: Very specific with check Pred
            self.check_pred(key, settings.BasisExpansion, fig,ax, box_opt=False, div_bin= div_bin)
            Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, ax, self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
            plt.legend()
            plt.savefig('theory/'+key+'_check_pred')
        return
         

    def compare(self, key, div_bin = False, box_opt = False):
        fig, axs = plt.subplots(4, 3, figsize=(16, 16), sharex=True)

        Plot.simple_plot(self, key, fig, axs[0,0], div_bin, box_opt)

        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'],self.measurement.hist_nom[key]['Values'], fig, axs[1,0],self.measurement.hist_nom[key]['Label'])

        # Given cn's
        self.check_pred(key, settings.BasisExpansion, fig,axs[2,0], box_opt=False, div_bin= div_bin, given_pred= True, full = True)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[2,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[2,1], box_opt=False, div_bin= div_bin, given_pred= True, full = False)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[2,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[2,2], box_opt=False, div_bin= div_bin, given_pred= True, full = True, made_up_norm= True)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[2,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        #Calculated cn's
        self.check_pred(key, settings.BasisExpansion, fig,axs[3,0], box_opt=False, div_bin= div_bin, given_pred= False, full = True)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[3,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[3,1], box_opt=False, div_bin= div_bin, given_pred= False, full = False)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[3,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[3,2], box_opt=False, div_bin= div_bin, given_pred= False, full = True, made_up_norm= True)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[3,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)


        plt.tight_layout()

        for ax in axs:
            for i in ax:
                i.set_title("")

        axs[2,0].set_title("Full prediction, Given Norm")
        axs[2,1].set_title("NNLLNNLO prediction, Given Norm")
        axs[2,2].set_title("Full prediction, Made up Norm")

        axs[2,0].set_ylabel("Given cn's")
        axs[3,0].set_ylabel("Calculated cn's")

        fig.suptitle(r""+self.measurement.exp_data[key]['Label'], fontsize=16, fontweight="bold")

        plt.savefig('compare/'+key+'_compare')

h = Plot()


data_tag_list = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]

for i in data_tag_list:
    h.compare(i)

#FIXME: Binned root histogram data has less bins than the experimental and prediction data, so I can't multiply it with the Smear-Matrix

#TODO: fit.config helpful info for smear matrix and amount of data which is fitted