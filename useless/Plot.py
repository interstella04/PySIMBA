import numpy as np
import matplotlib.pyplot as plt 
from Meas import Meas
from Meas import settings
from matplotlib import rc
from Tools import Tools


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 13})

class Plot:
    
    def __init__(self):
        self.measurement = Meas()
        self.collected_fits = Tools.PickleToDict('fit/collected_fits')
        self.just_sem = Tools.PickleToDict('fit/just_sem')
        
        
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


    def plot_histogram(self, bin_edges, values, fig, ax, title = '', xLabel = '', yLabel = '',div_bin = False, box_opt = False, color = 'darkgreen', legend_label = ''):

        #Option to divide the values by the width of the bin
        if(div_bin == True):
            for i in range(np.size(values)):
                values[i] = values[i]/(bin_edges[i+1]-bin_edges[i])

        ax.step(
        bin_edges,          
        [*values, 0],       # Values + additional 0 to close histogram
        where="post",       # Linienverlauf nach rechts
        color=color,       
        lw=1,
        label = '$\\mathrm{%s}$' % (legend_label)               
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

        return 
    

    
    def check_pred(self, key, end, fig, ax,an_pars, div_bin = False, box_opt = False):
        
        x,dx = Plot.simple_plot(self, key, fig, ax, div_bin, box_opt)

        cn = Meas.ConvertPars(an_pars)
        
        norm = an_pars[0]
            


        y = self.measurement.FullBsgPrediction(key,end, cn, norm)

        if key == 'belle': y*= 1/1000

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
    
    def calc_pred(self, key, an):
        cn = Meas.ConvertPars(an)
        norm = an[0]
        pred  = self.measurement.FullBsgPrediction(key,settings.BasisExpansion, cn, norm)
        #if key == 'belle': pred*= 1/1000

        return pred


    
    def compare(self, key, div_bin = False):
        fig, axs = plt.subplots(2, 5, figsize=(42*0.75, 20*0.75), sharex=True)

        start = 0
        if key == 'babar_incl': start = 1

        for ax in axs:
            for i in ax:
                i.set_ylabel("")

        for i in range(start, 2, 1):
            for j in range(0, 5, 1):

                axs[i,0].set_ylabel('with %d measurements ' % (i+2), fontsize = 20)

                number_meas = i+2
            
                #Plots experimental Data and Prediction
                self.check_pred(key, settings.BasisExpansion, fig,axs[i,j],self.collected_fits['%d' % (number_meas)]['results']['%d' % (j+3)]['an'], box_opt=False, div_bin= div_bin)
                #Plots Histogram Data from Root file
                Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[i,j], self.measurement.hist_nom[key]['Label'], div_bin=div_bin, legend_label= 'Fit Goal')
                #Plots difference between previous fit and my fit
                #Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], abs(self.measurement.hist_nom[key]['Values']-self.calc_pred(key,self.collected_fits['%d' % (number_meas)]['results']['%d' % (j+3)]['an'])), fig, axs[i,j], self.measurement.hist_nom[key]['Label'], div_bin=div_bin, color='blue', legend_label= '|Fit Goal - Pred. Data|')
                #Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], abs(self.measurement.exp_data[key]['dBFs']-self.calc_pred(key,self.collected_fits['%d' % (number_meas)]['results']['%d' % (j+3)]['an'] )), fig, axs[i,j], self.measurement.hist_nom[key]['Label'], div_bin=div_bin, color='hotpink', legend_label = '|Fit Goal - Exp. Data|')

        handles, labels = axs[1,4].get_legend_handles_labels()  # Nur einmal Labels abrufen  
        fig.legend(handles, labels, loc="upper left", ncol=2, fontsize=20)

        #fig.text(0.02, 0.8, "$\\mathrm{2 Measurements}$", va="center", ha="left", fontsize=30, fontweight="bold", rotation=90)
        #fig.text(0.02, 0.5, "$\\mathrm{3 Measurements}$", va="center", ha="left", fontsize=30, fontweight="bold", rotation=90)
        #fig.text(0.02, 0.2, "$\\mathrm{4 Measurements}$", va="center", ha="left", fontsize=30, fontweight="bold", rotation=90)

        #target_ax = axs[5, 2]

        # Grünen Rahmen setzen (spine = Rand der Achse)
        #for spine in target_ax.spines.values():
        #    spine.set_edgecolor("green")
        #    spine.set_linewidth(3)

        plt.tight_layout(rect=[0.05, 0, 1, 0.95])

        for ax in axs:
            for i in ax:
                i.set_title("")
                i.grid(axis="y", linestyle="--", alpha=0.6)

        axs[0,0].set_title("$\\mathrm{3 Parameter}$", fontsize = 30)
        axs[0,1].set_title("$\\mathrm{4 Parameter}$", fontsize = 30)
        axs[0,2].set_title("$\\mathrm{5 Parameter}$", fontsize = 30)
        axs[0,3].set_title("$\\mathrm{6 Parameter}$", fontsize = 30)
        axs[0,4].set_title("$\\mathrm{7 Parameter}$", fontsize = 30)

        fig.suptitle("$\\mathrm{%s}$" % self.measurement.exp_data[key]['Label'], fontsize=30, fontweight="bold")

        plt.savefig('compare/'+key+'_soft_compare')
    
    def just_one(self, key, div_bin = False):
        fig, axs = plt.subplots(2, 4, figsize=(24, 16), sharex=True)

        start = 0

        for ax in axs:
            for i in ax:
                i.set_ylabel("")

        for i in range(start, 2, 1):
            for j in range(0, 4, 1):

                if i%2 == 0:
                    sublead_or_not = 'subleading'
                    with_sub = True
                    axs[i,0].set_ylabel('with '+sublead_or_not, fontsize = 20)
                else:
                    sublead_or_not = 'leading'
                    with_sub = False
                    axs[i,0].set_ylabel('just '+sublead_or_not, fontsize = 20)

                #Plots experimental Data and Prediction
                self.check_pred(key, settings.BasisExpansion, fig,axs[i,j],self.just_sem['1'][sublead_or_not]['%d' % (j+4)]['an'], box_opt=False, div_bin= div_bin, with_sub=with_sub)
                #Plots Histogram Data from Root file
                Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], self.measurement.hist_nom[key]['Values'], fig, axs[i,j], self.measurement.hist_nom[key]['Label'], div_bin=div_bin, legend_label= 'Fit Goal')
                #Plots difference between previous fit and my fit
                Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], abs(self.measurement.hist_nom[key]['Values']-self.calc_pred(key,self.just_sem['1'][sublead_or_not]['%d' % (j+4)]['an'],with_sub=with_sub )), fig, axs[i,j], self.measurement.hist_nom[key]['Label'], div_bin=div_bin, color='blue', legend_label= '|Fit Goal - Pred. Data|')
                #Plot.plot_histogram(self, self.measurement.hist_nom[key]['Bins'], abs(self.measurement.exp_data[key]['dBFs']-self.calc_pred(key,self.collected_fits['%d' % (number_meas)][sublead_or_not]['%d' % (j+4)]['an'],with_sub=with_sub )), fig, axs[i,j], self.measurement.hist_nom[key]['Label'], div_bin=div_bin, color='hotpink', legend_label = '|Fit Goal - Exp. Data|')

        handles, labels = axs[0,0].get_legend_handles_labels()  # Nur einmal Labels abrufen  
        fig.legend(handles, labels, loc="upper left", ncol=2, fontsize=20)

        fig.text(0.02, 0.5, "$\\mathrm{1 Measurement}$", va="center", ha="left", fontsize=30, fontweight="bold", rotation=90)

        plt.tight_layout(rect=[0.05, 0, 1, 0.95])

        for ax in axs:
            for i in ax:
                i.set_title("")
                i.grid(axis="y", linestyle="--", alpha=0.6)

        axs[0,0].set_title("$\\mathrm{4 Parameter}$", fontsize = 30)
        axs[0,1].set_title("$\\mathrm{5 Parameter}$", fontsize = 30)
        axs[0,2].set_title("$\\mathrm{6 Parameter}$", fontsize = 30)
        axs[0,3].set_title("$\\mathrm{7 Parameter}$", fontsize = 30)

        fig.suptitle("$\\mathrm{%s}$" % self.measurement.exp_data[key]['Label'], fontsize=30, fontweight="bold")

        plt.savefig('compare/just_'+key)

h = Plot()

#data_tag_list = ["babar_hadtag", "babar_incl"]
data_tag_list = ["belle"]
#data_tag_list = ["babar_sem"]

for i in data_tag_list:
    h.compare(i)
    print("II")

#h.just_one('babar_sem')

#FIXME: Binned root histogram data has less bins than the experimental and prediction data, so I can't multiply it with the Smear-Matrix

#TODO: fit.config helpful info for smear matrix and amount of data which is fitted