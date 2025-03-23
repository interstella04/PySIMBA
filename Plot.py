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
    
    def check_pred(self, key, end, fig, ax,an_pars, div_bin = False, box_opt = False):
        
        x,dx = Plot.simple_plot(self, key, fig, ax, div_bin, box_opt)

        cn = Meas.ConvertPars(an_pars)
        
        norm = an_pars[0]
            


        y = self.measurement.BsgPrediction_full(key,end, cn, norm)

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
        fig, axs = plt.subplots(6, 3, figsize=(16, 16), sharex=True)

        as23 = np.array([0.838658,-0.422353,0.106128])
        as24 = np.array([0.914286,-0.504880,0.249589,0.207707])
        as25 = np.array([0.042643,26.878645,24.338372,22.634102,16.416970])
        al23 = np.array([0.359303,3947.576246,80.378891])
        al24 = np.array([0.437288,1.694578,-3.864314,-0.074478])
        al25 = np.array([0.671518,-1.099613,6.571752,0.371206,-0.031937])

        al33 = np.array([0.148084,6016.698111,1252.626045])
        al34 = np.array([0.502049,31806.331219,17092.225881,-1177.898593])
        al35 = np.array([0.493309,-2.817043,20.213729,-23.930707,-0.565804])
        as33 = np.array([0.844866,-0.407248,0.058612])
        as34 = np.array([0.922838,-0.415097,0.152539,0.104518])
        as35 = np.array([0.951869,-0.486525,0.203357,0.195060,0.113255])

        as43 = np.array([0.843780,-0.402125,0.062842])
        as44 = np.array([0.920390,-0.439685,0.177431,0.140803])
        as45 = np.array([0.951448,-0.512677,0.211717,0.213323,0.118997])
        al43 = np.array([0.155194,25038.419489,5364.382903])
        al44 = np.array([0.382562,-0.494960,0.009761,0.127197])
        al45 = np.array([0.520989,-3.558982,30.134934,-25.682246,-0.925172])


        #sub, 2
        self.check_pred(key, settings.BasisExpansion, fig,axs[0,0],as23, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[0,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[0,1],as24, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[0,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        self.check_pred(key, settings.BasisExpansion, fig,axs[0,2],as25, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[0,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        #lead, 2
        self.check_pred(key, settings.BasisExpansion, fig,axs[1,0],al23, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[1,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[1,1],al24, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[1,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        self.check_pred(key, settings.BasisExpansion, fig,axs[1,2],al25, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[1,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        #sub, 3
        self.check_pred(key, settings.BasisExpansion, fig,axs[2,0],as33, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[2,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[2,1],as34, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[2,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        self.check_pred(key, settings.BasisExpansion, fig,axs[2,2],as35, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[2,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        #lead, 3
        self.check_pred(key, settings.BasisExpansion, fig,axs[3,0],al33, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[3,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[3,1],al34, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[3,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        self.check_pred(key, settings.BasisExpansion, fig,axs[3,2],al35, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[3,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        #sub, 4
        self.check_pred(key, settings.BasisExpansion, fig,axs[4,0],as43, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[4,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[4,1],as44, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[4,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        self.check_pred(key, settings.BasisExpansion, fig,axs[4,2],as45, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[4,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        #lead, 4
        self.check_pred(key, settings.BasisExpansion, fig,axs[5,0],al43, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[5,0], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)

        self.check_pred(key, settings.BasisExpansion, fig,axs[5,1],al44, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[5,1], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
        self.check_pred(key, settings.BasisExpansion, fig,axs[5,2],al45, box_opt=False, div_bin= div_bin)
        Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[5,2], self.measurement.hist_nom[key]['Label'], div_bin=div_bin)
        
       
        plt.tight_layout()

        for ax in axs:
            for i in ax:
                i.set_title("")
                i.set_ylabel("")

        axs[0,0].set_title("size 3")
        axs[0,1].set_title("size 4")
        axs[0,2].set_title("size 5")



        #fig.suptitle(r""+self.measurement.exp_data[key]['Label'], fontsize=16, fontweight="bold")

        plt.savefig('compare/'+key+'_compare')

h = Plot()


data_tag_list = ["babar_incl", "babar_hadtag", "babar_sem", "belle"]

for i in data_tag_list:
    h.compare(i)

#FIXME: Binned root histogram data has less bins than the experimental and prediction data, so I can't multiply it with the Smear-Matrix

#TODO: fit.config helpful info for smear matrix and amount of data which is fitted