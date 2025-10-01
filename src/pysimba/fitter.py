from .tool import Tools
from iminuit import Minuit
from .theory import Theory
from .theory import settings


class Fitter:
    def __init__(self):
        m: Minuit

        mb: float
        chisq: float

        theo: Theory

        NumbOfPar: int

        Lambda: float

    # NumbPar: How many parameters will be fitted, so how many of start_pars will be used
    # with_minos: Calculate the fit with or without minos
    def DoSingleFit(self, NumbPar: int = 3, with_minos: bool = True):
        self.NumbOfPar = NumbPar

        print("Fit for: " + str(settings.KeyOrder) + " using %d parameters" % NumbPar)

        mes = Theory()  # Initialize Theory Object which will be fitted

        # Set up iminuit fit-object
        n = Minuit(
            mes.Chisq,
            settings.config["StartValues"][0:NumbPar],
            name=settings.config["FitVars"][0:NumbPar],
        )

        # n.tol = 0.1 (standart value)
        n.errordef = 0.0001

        # Calculate fit
        n.migrad(ncall=100000, use_simplex=True)
        n.hesse()

        # Save chisq and mb for later (Result class)
        self.chisq = mes.chisq
        self.mb = mes.mb

        # Apply minos algorithm
        if with_minos:
            n.minos()

        self.m = n
        self.theo = mes

        # Turns Basis Expansion into a float for later usage
        self.Lambda = Tools.StrToLambda(settings.BasisExpansion)

        return
