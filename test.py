import numpy as np
import matplotlib.pyplot as plt 
from Meas import Meas
from Meas import settings
from matplotlib import rc
from Tools import Tools


collected_fits = Tools.pickle_to_dict('fit/collected_fits')
#collected_fits['%d' % (number_meas)][sublead_or_not]['%d' % (j+4)]['an']
with open('table.csv', 'w') as file:
    for numb_meas in range(2,5,1):
        for with_without in ['subleading', 'leading']:
            for numb_para in range(4,8,1):
                an = collected_fits['%d' % (numb_meas)][with_without]['%d' % (numb_para)]['an']
                cn = Meas.ConvertPars(an)
                norm = an[0]
                chi2 = collected_fits['%d' % (numb_meas)][with_without]['%d' % (numb_para)]['Chisq']
                mass = collected_fits['%d' % (numb_meas)][with_without]['%d' % (numb_para)]['mb']

                file.write("%d&%s&%d&" % (numb_meas, with_without,numb_para))
                file.write(r"%.2f&%.4f&%.4f\\" % (chi2, mass, norm))
                file.write('\n')

with open('cn_pars.csv', 'w') as f:
    for numb_meas in range(2,5,1):
        for with_without in ['subleading', 'leading']:
            for numb_para in range(4,8,1):
                an = collected_fits['%d' % (numb_meas)][with_without]['%d' % (numb_para)]['an']
                cn = Meas.ConvertPars(an)
                f.write("%d&%.1s&%d" % (numb_meas, with_without,numb_para))
                for val in cn:
                    line = " ".join(["%.3f" % val])
                    f.write('&'+ line)
                
                f.write(r'\\' + '\n')