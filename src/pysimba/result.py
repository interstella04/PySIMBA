# Other packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# This package
from .theory import Theory
from .theory import settings
from .tool import Tools
from .fitter import Fitter
from . import BASE_DIR

# Setting the text font on the plots
rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
rc("text", usetex=True)
plt.rcParams.update({"font.size": 13})


class Result:
    def __init__(self):
        # Dictionary with experimental data
        self.ExpData = Tools.PickleToDict(settings.config["MeasurementPath"])

    def PrepareExpData(self, key: str, div_bin: bool = False):
        # Calculate the middle of the bin
        # x error (dx) is bin width
        x = np.zeros(np.size(self.ExpData[key]["dBFs"]))
        dx = np.zeros(np.size(self.ExpData[key]["dBFs"]))
        for i in range(np.size(self.ExpData[key]["dBFs"])):
            dx[i] = (
                self.ExpData[key]["Bins"][i + 1] - self.ExpData[key]["Bins"][i]
            ) / 2
            x[i] = dx[i] + self.ExpData[key]["Bins"][i]

        # Divide by the width of the bin or not

        scale = settings.config["Scale"].get(
            key, 1
        )  # If there is a scale defined in settings.yml take that, else scale = 1

        if not div_bin:
            y = self.ExpData[key]["dBFs"]
            dy = self.ExpData[key]["dBFErrors"]
        else:
            y = self.ExpData[key]["dBFs"] / (2 * dx)
            dy = self.ExpData[key]["dBFErrors"] / (2 * dx)

        y *= scale
        dy *= scale

        # return the Label of the Measurement to put in the title
        title = self.ExpData[key]["Label"]

        if key == "belle":
            ylabel = r"$\Delta B \left( B \rightarrow X_S \gamma \right) \: [10^{3}  / 50 \mathrm{ MeV}]$"
        else:
            ylabel = r"$\Delta B \left( B \rightarrow X_S \gamma \right) \: [10^{-4}  / 0.1 \mathrm{ GeV}]$"

        return x, y, dx, dy, title, ylabel

    def SimplePlot(
        ax,
        x,
        y,
        dx=None,
        dy=None,
        title: str = None,
        ylabel: str = None,
        box_opt: bool = False,
        color: str = "black",
        min_y=None,
    ):
        # Set Box options
        if box_opt:
            textstr = "\n".join(
                (
                    r"$\mathrm{Entries:} \\: %d$" % (np.size(x),),
                    r"$\mathrm{Mean:} \\: %0.3f$" % (x.mean(),),
                    r"$\mathrm{Std Dev:} \\: %0.4f$" % (x.std(),),
                )
            )
            ax.text(
                0.8,
                0.95,
                textstr,
                fontsize=12,
                color="black",
                bbox=dict(facecolor="white", alpha=1, edgecolor="black"),
                transform=ax.transAxes,
            )

        ax.errorbar(
            x,
            y,
            dy,
            dx,
            fmt="o",
            markersize=5,
            label=r"$\mathrm{Exp. Data}$",
            color=color,
        )
        ax.set_ylabel(ylabel, fontsize=16)

        ax.set_xlabel(r"$E_{\gamma}[\mathrm{GeV}]$", fontsize=16, loc="right")
        ax.set_title(r"$\mathrm{%s}$" % title, fontsize=20)

        ax.set_ylim(bottom=min_y)

        return

    def SimpleHistogram(ax, bins, values, err=None, color: str = "black"):
        ax.step(
            bins,
            [*values, 0],  # Values + additional 0 to close histogram
            where="post",  # Linienverlauf nach rechts
            color=color,
            lw=1,
            label=r"$\mathrm{FIT}$",
        )

        ax.set_xlim(left=bins[0], right=bins[-1])

        return

    # Calculates Fit for a Number of Parameters
    def CalculateFitObject(NumbPar: int = 3) -> Fitter:
        fit = Fitter()
        fit.DoSingleFit(NumbPar, with_minos=True)

        return fit

    def CalculatePrediction(self, key: str, fit: Fitter) -> list[float]:
        an = np.array(fit.m.values)
        norm = an[0]

        cn = Theory.ConvertPars(an)

        pred = fit.theo.FullBsgPrediction(key, settings.BasisExpansion, cn, norm)

        # Apply scale
        scale = settings.config["Scale"].get(key, 1)
        pred *= scale

        return pred

    def AddFitInformation(ax, fit: Fitter):
        # Information
        textstr = "\n".join(
            (
                "$ \\ %s $" % (Tools.make_c_string(fit.NumbOfPar),),
                "$\\ \\chi^2 = \\: %0.2f$" % (fit.chisq,),
                "$\\ m_b = \\: %0.3f \\ \\mathrm{GeV}$" % (fit.mb,),
                "$ \\ \\lambda = \\: %0.3f \\ \\mathrm{GeV}$" % (fit.Lambda,),
            )
        )
        # Box Placement and Options
        ax.text(
            0.82,
            0.87,
            textstr,
            fontsize=12,
            color="black",
            bbox=dict(facecolor="white", alpha=1, edgecolor="black"),
            transform=ax.transAxes,
        )

        # Add Legend
        ax.legend()

        return

    def Plot(
        self, key: str, fit_obj: Fitter, div_bin: bool = False, box_opt: bool = False
    ):
        print("Plotting %s ... " % key)

        # Plots just experimental Data
        fig, ax = plt.subplots(figsize=(8, 6))

        x, y, dx, dy, title, ylabel = Result.PrepareExpData(self, key, div_bin)
        Result.SimplePlot(ax, x, y, dx, dy, title, ylabel, box_opt)

        filename_exp = settings.config["ExpPath"] + key

        fig.savefig(BASE_DIR / filename_exp, dpi=500)
        plt.close()

        # Plots the ExpData with the calculated Prediction
        f, a = plt.subplots(figsize=(8, 6))
        pred = Result.CalculatePrediction(self, key, fit_obj)

        Result.SimpleHistogram(a, self.ExpData[key]["Bins"], pred, color="red")
        Result.SimplePlot(a, x, y, dx, dy, title, ylabel, box_opt=False, min_y=0)
        Result.AddFitInformation(a, fit_obj)

        filename_fit = (
            settings.config["ResultPath"] + settings.config["Tag"] + "_" + key
        )

        f.savefig(BASE_DIR / filename_fit, dpi=500)
        plt.close()

        return

    def PrintResultToTxt(fit_obj: Fitter):
        print("Writing ...")

        filename = settings.config["ResultPath"] + settings.config["Tag"] + ".txt"
        path = BASE_DIR / filename

        with open(path, "w") as file:
            file.write("Results for %d Parameters\n" % fit_obj.NumbOfPar)

            # Fit Results
            file.write("Chisq = %.4f \n" % fit_obj.chisq)
            file.write("mb = %.4f \n" % fit_obj.mb)

            # Write an's to file
            for i in range(0, fit_obj.NumbOfPar, 1):
                file.write(
                    "%s = %.5f \n"
                    % (settings.config["FitVars"][i], fit_obj.m.values[i])
                )

            file.close()
        return

    def ShowResults(
        self, fit_obj: Fitter, div_bin: bool = False, box_opt: bool = False
    ):
        Result.PrintResultToTxt(fit_obj)

        for key in settings.KeyOrder:
            self.Plot(key, fit_obj, div_bin, box_opt)

        return

    # Calculate the Fit, depending on config file and show the Results in the Result folder
    def Run():
        res = Result()

        fit_obj = Result.CalculateFitObject(settings.config["NumbPar"])

        res.ShowResults(fit_obj, False, False)

        return
