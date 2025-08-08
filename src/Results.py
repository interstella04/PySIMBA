import numpy as np
import matplotlib.pyplot as plt 
from Meas import Theory
from Meas import settings
from matplotlib import rc
from Tools import Tools
from fitter import Fitter


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 13})

class result:

    def __init__(self):
        #Dictionary with experimental data
        self.ExpData = Tools.PickleToDict(settings.config["MeasurementPath"])
        
        #Histogram Dictionary from Root Data
        self.HistNomFit = Tools.PickleToDict("../theory/hist_nom")

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


    def SimplePlot(ax, x, y, dx = None, dy = None, label: str = None, box_opt: bool = False, color: str = 'black', min_y = None):

        # Set Box options
        if (box_opt):
            textstr = '\n'.join((
                '$\\mathrm{Entries:} \\: %d$' % (np.size(x),),
                '$\\mathrm{Mean:} \\: %0.3f$' % (x.mean(),),
                '$\\mathrm{Std Dev:} \\: %0.4f$' % (x.std(),)))
            ax.text(0.8, 0.95, textstr, fontsize=12, color="black", bbox=dict(facecolor="white", alpha=1, edgecolor="black"), transform = ax.transAxes)
        
        
        ax.errorbar(x, y, dy, dx, fmt='o', markersize = 5, label = '$\\mathrm{Exp. Data}$', color = color )
        ax.set_ylabel('$\\Delta B \\left( B \\rightarrow X_S \\gamma \\right) \\: [10^{-4}  / 0.1 \\mathrm{GeV}]$', fontsize = 16)
              
        ax.set_xlabel('$E[\\mathrm{GeV}]$', fontsize = 16, loc = "right")
        ax.set_title('$\\mathrm{%s}$' % label, fontsize = 20)

        ax.set_ylim(bottom = min_y)

        return
    
    def SimpleHistogram(ax, bins, values, err = None, color: str = 'black'):
        
        ax.step(
        bins,          
        [*values, 0],       # Values + additional 0 to close histogram
        where="post",       # Linienverlauf nach rechts
        color=color,       
        lw=1,
        label = '$\\mathrm{FIT}$'               
        )

        ax.set_xlim(left = bins[0], right = bins[-1])

        return
    
    # Calculates Fit for a Number of Parameters
    def CalculateFitObject(NumbPar: int = 3) -> Fitter:
        fit = Fitter()
        fit.DoSingleFit(NumbPar, with_minos= True)
        
        return fit
    
    def CalculatePrediction(self, key: str,  fit: Fitter) -> list[float]:
        an = np.array(fit.m.values)
        norm = an[0]

        cn = Theory.ConvertPars(an)

        pred = fit.theo.FullBsgPrediction(key, settings.BasisExpansion, cn, norm)

        return pred
    
    def AddFitInformation(ax, fit: Fitter):
        # Information
        textstr = '\n'.join((
            '$\\ \\chi^2 = \\: %0.2f$' % (fit.chisq,),
            '$\\ m_b = \\: %0.3f$' % (fit.mb,),
            '$ \\ n_p = \\: %d$' % (fit.NumbOfPar,),
            '$ \\ \\lambda = \\: %0.3f$' % (fit.Lambda,)))
        # Box Placement and Options
        ax.text(0.9, 0.95, 
                textstr, 
                fontsize=12, 
                color="black", 
                bbox=dict(facecolor="white", alpha=1, edgecolor="black"), 
                transform = ax.transAxes)
        
        # Add Legend
        ax.legend()
        
        return
        
         
    def Plot(self, key: str, fit_obj: Fitter,  div_bin: bool = False, box_opt: bool = False):

        print("Plotting %s ... " % key)

        # Plots just experimental Data
        fig, ax = plt.subplots(figsize=(8, 6))

        x, y, dx, dy, label = result.PrepareExpData(self, key, div_bin)
        result.SimplePlot(ax, x, y, dx, dy, label, box_opt)
        fig.savefig('../data/' + key, dpi = 500)
        plt.close()


        # Plots the ExpData with the calculated Prediction
        f, a = plt.subplots(figsize=(8, 6))
        pred = result.CalculatePrediction(self, key, fit_obj)

        result.SimpleHistogram(a, self.ExpData[key]['Bins'], pred, color= 'red')
        result.SimplePlot(a, x, y, dx, dy, label, box_opt = False, min_y= 0)
        result.AddFitInformation(a, fit_obj)
        f.savefig(settings.config["ResultPath"] + settings.config["Tag"] + '_' + key, dpi = 500)
        plt.close()

        return
    
    def PrintResultToTxt(fit_obj: Fitter):
        print('Writing ...')

        with open(settings.config["ResultPath"] + settings.config['Tag'] + '.txt', 'w') as file:
            file.write('Results for %d Parameters\n' % fit_obj.NumbOfPar)

            # Fit results
            file.write('Chisq = %.4f \n' % fit_obj.chisq)
            file.write('mb = %.4f \n' % fit_obj.mb)

            # Write an's to file 
            for i in range(0, fit_obj.NumbOfPar, 1):
                file.write('%s = %.5f \n' % (settings.config["FitVars"][i], fit_obj.m.values[i]))

            file.close()
        return

    
    def ShowResults(self, fit_obj: Fitter, div_bin: bool = False, box_opt: bool = False):

        result.PrintResultToTxt(fit_obj)

        for key in settings.KeyOrder:
            self.Plot(key, fit_obj, div_bin, box_opt)

        return
    
    # Calculate the Fit, depending on config file and show the results in the result folder
    def Run():
        res = result()

        fit_obj = result.CalculateFitObject(settings.config["NumbPar"])

        res.ShowResults(fit_obj, False, False)

        return

result.Run()