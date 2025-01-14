import uproot
import numpy as np
from Meas import Meas
import matplotlib.pyplot as plt
from Tools import Tools


Fmn_moments = Tools.pickle_to_dict("theory/Fmn_moments")

print(np.size(Fmn_moments['expx3']['Moment']))