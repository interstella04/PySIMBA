import numpy as np
import numpy.typing as npt

from Tools import Tools
from dataclasses import dataclass

# For config File
import yaml
from pathlib import Path

#BASE_DIR = Path(__file__).resolve().parent
#settings_path = BASE_DIR.parent / 'settings.yml'

@dataclass
class settings:

    # Open the yaml (config) file and store the information in config (dictionary)
    with open('../settings.yml', "r") as f: 
        config = yaml.safe_load(f)

    # Read out Information of config
    SubLeadCoefficients = config["SubLeadCoefficients"]
    TheoryOrder = config["TheoryOrder"]
    SubLeadTheoryOrder = config["SubLeadTheoryOrder"]
    FitVars = config["FitVars"]
    KeyOrder = config["KeyOrder"]
    BasisExpansion = config["BasisExpansion"]
    SubLeadBasisExpansion = config["SubLeadBasisExpansion"]
    
    rho2: float = config["Constants"]["rho2"]
    mB  : float = config["Constants"]["mB"]
    La2 : float = config["Constants"]["La2"]
    N0  : float = config["Constants"]["N0"]
    VtbVts: float = config["Constants"]["VtbVts"]
    C2C7: float = config["Constants"]["C2C7"]
    C2C2: float = config["Constants"]["C2C2"]
    C8C7: float = config["Constants"]["C8C7"]
    C8C8: float = config["Constants"]["C8C8"]
    C2C8: float = config["Constants"]["C2C8"]


class Theory:

    # Define types
    type Vector = list[float]
    type npVector = npt.NDArray[np.float64]

    def __init__(self):

        # Reading out the dictionaries from pickles at the paths given in the yaml (config)
        #Theory Dictionaries
        self.Theory = Tools.PickleToDict(settings.config["TheoryPath"])
        self.SubleadTheory = Tools.PickleToDict(settings.config["SubleadingTheoryPath"])

        #Dictionary with experimental data
        self.ExpData = Tools.PickleToDict(settings.config["MeasurementPath"])
        
        #Moments Dictionary
        self.FmnMoments = Tools.PickleToDict(settings.config["TheoryMomentsPath"])

        # Variables that will get a value in the process
        self.chisq: float
        self.mb: float 
        self.M1: float
        self.M2: float
        self.M3: float
        self.la: float
    
    # Calculates the BsgPrediction for every given Theory in TheoryOrder and the SubLeadingBsgPrediction for every given SubLeadTheoryOrder. Results into the final Prediction with leading and subleading Theory combined
    def FullBsgPrediction(self, key: str, end: str, c_n_params: npVector, norm: float) -> npVector:

        pred = np.zeros(np.size(self.ExpData[key]['dBFs']))
        self.mb = self.mb1SPrediction(c_n_params)

        # Full leading prediction
        for order in range(np.size(settings.TheoryOrder)):
            pred += self.BsgPrediction(key, settings.TheoryOrder[order], end, c_n_params, norm ) * self.TheoryPrefactor(settings.TheoryOrder[order], norm)

        # Full subleading prediction
        for order in range(np.size(settings.SubLeadTheoryOrder)):
            pred += self.BsgSubLeadingPrediction(key, self.SubLeadPars(c_n_params, settings.SubLeadCoefficients[order]), norm) * self.TheoryPrefactor(settings.SubLeadTheoryOrder, norm)

        # Multiply with Smearing-Matrice
        pred = np.matmul(self.ExpData[key]['Smear'],pred) 

        return pred
        
    
    # Prefactor for leading and subleading theory
    def TheoryPrefactor(self, TheoryOrder: str, norm: float) -> float:

        value = settings.N0

        if ('22' in TheoryOrder):
            value *= settings.C2C2
            value /= norm
            value *= np.power(settings.VtbVts * self.mb, 2.0)
        elif ('SSF27' in TheoryOrder):
            value *= settings.C2C7 / np.sqrt(norm) * settings.VtbVts * self.mb
            value *= settings.La2 / self.mb
        elif ('27' in TheoryOrder):
            try:
                value *= settings.C2C7 / np.sqrt(norm) * settings.VtbVts * self.mb
            except RuntimeWarning:
                print('Norm for 27%s :' % (TheoryOrder) + norm)
        elif ('28' in TheoryOrder):
            value *= settings.C2C8 / norm * np.power(settings.VtbVts * self.mb, 2.0)
        elif ('88' in TheoryOrder):
            value *= settings.C8C8
            value /= norm
            value *= np.power(settings.VtbVts * self.mb, 2.0)
        elif ('78'  in TheoryOrder):
            value *= settings.C8C7
            value /= np.sqrt(norm)
            value *= settings.VtbVts * self.mb
        
        return value


    # Returns a prediction array for one specific experiment with one specific lambda
    def BsgPrediction(self, key: str, mid: str, end: str, c_n_params: npVector, norm: float) -> npVector:      
        pred = np.array([])
        for i in range(np.size(self.ExpData[key]['dBFs'])):
            value = 0
            for j in range(np.size(c_n_params)):
                for k in range(np.size(c_n_params)):
                    value += self.Theory[key][mid]['la'+end]['Values'][i][j][k] * c_n_params[j] * c_n_params[k] 
            
            value *= norm
            pred = np.append(pred, value)
        return pred 


    # Returns a the subleading prediction for a given measurement
    def BsgSubLeadingPrediction(self, key:str, c_n_params: npVector, norm: float) -> npVector:
        pred = np.array([])
        for i in range(np.size(self.ExpData[key]['dBFs'])):
            value = 0
            for j in range(np.size(c_n_params)):
                value +=  self.SubleadTheory[key]['la'+ settings.SubLeadBasisExpansion]['Values'][i][0][j] * c_n_params[j]
            value *= norm
            pred = np.append(pred, value)

        return pred
    
    # puts Zoltan's SF into the prediction
    # converts c_n's into subleading coefficients
    def SubLeadPars(self, c_n_params: npVector, d2: float) -> npVector:
        
        Rho2 = settings.rho2
        mB = settings.mB
        mb = self.mb1SPrediction(c_n_params)
        la = Tools.StrToLambda(settings.SubLeadBasisExpansion)
        Lambda2 = settings.La2


        if('SSF27_1' in settings.SubLeadTheoryOrder):
            x = (0.6514810199386504 * (mB - mb) - 0.8686413599182005 * la + 0.3257405099693252 * Rho2 / Lambda2) / (1.8531015678254943 * la - d2 * la)
        elif('SSF27_2' in settings.SubLeadTheoryOrder):
            x = (0.4722982832332954 * (mB - mb) - 0.5667579398799545 * la + 0.2361491416166477 * Rho2 / Lambda2) / (1.3695121740357048 * la - d2 * la)
        else:
            print('Unknown subleading shape function -- please specify Theory.SubLeadPars()')

        #calculating d0, d1 and d2 and returning them in a numpy array
        return np.array([1.-x, x* (1-d2), x*d2])
    
    
    # Calculate shape fuction moment of order 'order'
    def Moment(self, c_n_params: npVector, order: int) -> float:
        if (order > np.size(self.FmnMoments[settings.config["TheoryTag"]]['Moment'])): print('Theory.Moment(): requested an order of moment, that cannot be calculated')
        value = 0

        # Prepares Moment, using the lambda for which the fit is currently calculated
        prepared_moment = np.pow(Tools.StrToLambda(settings.SubLeadBasisExpansion), order) * self.FmnMoments[settings.config["TheoryTag"]]['Values'][order]
        for i in range(np.size(c_n_params)):
            for j in range(np.size(c_n_params)):
                value += prepared_moment[i][j] * c_n_params[i] * c_n_params[j]
        return value
    
    # Calculate mb for a given set of cn's
    def mb1SPrediction(self, c_n_params: npVector) -> float:
        mb = 0
        M1 = self.Moment(c_n_params, 1)
        M2 = self.Moment(c_n_params, 2)
        M3 = self.Moment(c_n_params, 3)

        Lambda2 = settings.La2
        mB = settings.mB
        Rho2 = settings.rho2

        p = -18 * Lambda2 - 25 * M1 * M1 + 36 * M2 - 22 * M1 * mB + 11 * mB * mB
        q = 5 * M1 * (27 * Lambda2 + 25 * M1 * M1 - 54 * M2) - 3 * (45 * Lambda2 - 55 * M1 * M1 + 72 * M2) * mB + 51 * M1 * mB * mB - 17 * mB * mB * mB + 162 * (M3 + Rho2)
        u = np.cbrt((-q + np.sqrt(q * q + p * p * p)))
        mb = 1. / 6. * (5 * (mB - M1) - p / u + u)

        return mb
    
    # Calculate lambda for a given set of cn's
    def lambda11SPrediction(self, c_n_params: npVector) -> float:
        M1 = self.Moment(c_n_params, 1)
        lambda2 = settings.La2
        mB = settings.mB

        mb = self.mb1SPrediction(c_n_params)
        lambda1 = 3. * lambda2 + 2. * mb * (mB - M1 - mb)

        return lambda1
    
    #converts an's into cn's
    def ConvertPars(par: Vector) -> Vector:
        length = np.size(par)-1
        if(length < 0):
            length = np.size(settings.FitVars)-1

        annorm = 1.0 + np.sum(par[1:]**2)

        cn = [1.0/np.sqrt(annorm)]

        for i in range(length):
            cn.append(par[i+1]/np.sqrt(annorm))

        return cn

    # The given parameters par are an's which will be converted into cn's
    # Returns the chisq of the current prediction. It's the function we minimize to get the 'par' as fit results
    def Chisq(self, par: Vector) -> float:
        cn = Theory.ConvertPars(par)
        norm = par[0]

        pred_glob = np.array([])
        meas_glob = np.array([])

        # Sets Size of the global matrices, depending on the size of the different measurements
        size = 0
        minimum = {}
        for key in settings.KeyOrder:
            # If there is no minimum set in the config file for a measurement, then the minimum is zero
            minimum[key] = settings.config["Minimum"].get(key, 0)

            size += np.size(self.ExpData[key]['dBFs'])-minimum[key]

        Cov_glob = np.zeros((size,size))

        ntot = 0

        for key in settings.KeyOrder:

            # Determine # of Bins
            nbins = np.size(self.ExpData[key]['dBFs'])
            min = minimum[key]
            max = nbins # No Option yet to have a different Maximum. Could be implemented similar to the minimum though

            # Construct individual matrices to store predictions and measurements
            full_pred = Theory.FullBsgPrediction(self,key, settings.BasisExpansion, cn, norm)
            pred =  full_pred[min:max+1]
            meas = self.ExpData[key]['dBFs'][min:max+1]
            Cov = self.ExpData[key]['Cov'][min:max+1,min:max+1]
            
            # Add to global prediction
            pred_glob = np.append(pred_glob, pred)
            meas_glob = np.append(meas_glob, meas)
            Tools.SetSub(Cov_glob, ntot, ntot, Cov)
            ntot += max-min
        
        Cov_glob = np.array(np.linalg.inv(Cov_glob)) 
        
        Cov_glob = Tools.MakePositiveDefinite(Cov_glob, eps = 0) # Forces Matrice to be positive definite


        #Calculate chi2
        Dummy = np.dot((meas_glob - pred_glob).transpose(), Cov_glob)
        Chisq = np.dot(Dummy, (meas_glob - pred_glob))

        self.chisq = Chisq

        self.mb = self.mb1SPrediction(cn) # Store current mb & moments
        self.M1 = self.Moment(cn,1)
        self.M2 = self.Moment(cn,2)
        self.M3 = self.Moment(cn,3)
        self.la = settings.BasisExpansion 
        
        return Chisq