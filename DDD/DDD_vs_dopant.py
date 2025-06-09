import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import re


def fitting_function(x,C,Dx):
    remaining_factor = 1-C*np.log10(1+x/Dx)
    return remaining_factor

dictionary = {
    "CdSeTe:Cu" : (0.632,2.06e11),
    "CdSeTe:As" : (0.805,2.19e11),
    "GaAs/Ge" : (0.322,1.86e9)
}

colorList = ["black","crimson", "blue"]

fig,ax = plt.subplots()
# props = dict(boxstyle='round', facecolor='white', alpha=0.5)
# fig.text(0.175,0.25, r"$\frac{\eta_1}{\eta_0}=1-C*log(1+\frac{DDD}{D_x})$", fontsize = 15,bbox=props)
curve_fit_abscissa = np.logspace(7,13,100)
for index, key in enumerate(dictionary.keys()):
    color = colorList[index]
    popt = dictionary[key]
    ax.plot(curve_fit_abscissa,fitting_function(curve_fit_abscissa,*popt),color = color, label = key, linewidth=2)

ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
ax.set_xscale("log")
ax.set_ylim(0,1.1)
ax.xaxis.set_tick_params(labelsize = 15)
ax.yaxis.set_tick_params(labelsize = 15)
ax.legend()
ax.set_xscale('log')

plt.savefig(r"C:\Users\Zach\Downloads\Figure_1.png", dpi = 400)