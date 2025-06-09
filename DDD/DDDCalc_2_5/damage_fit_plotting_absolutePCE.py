import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def fitting_function(x,C,Dx,BoL):
    remaining_factor = (1-C*np.log10(1+x/Dx))*BoL
    return remaining_factor


# Hard coded ddd curve parameters from most recent data set
fit_coef_dict = {
    "CdSeTe:Cu" : (.656,2.76e11,1),
    "CdSeTe:As" : (.551,9.82e10, 1),
    "Azur 3G28" : (.249,1.73e9, 1),
    "Azur 4G32" : (.359,6.94e9, 1),
    "CESI CTJ-LC" : (.354,4.81e9,1)
}

fig,ax = plt.subplots()

for key in fit_coef_dict.keys():
    curve_fit_abscissa = np.logspace(7,13,100)
    C = fit_coef_dict[key][0]
    Dx = fit_coef_dict[key][1]
    BoL_PCE = fit_coef_dict[key][2]
    ax.plot(curve_fit_abscissa,fitting_function(curve_fit_abscissa,C,Dx,BoL_PCE), label = key, linewidth=2)

ax.set_ylabel("PCE Remaining Factor", weight = 'bold', fontsize = 15)
ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
ax.set_xscale("log")
ax.set_ylim(0,1.1)
ax.xaxis.set_tick_params(labelsize = 12)
ax.yaxis.set_tick_params(labelsize = 12)
fig.tight_layout()
#ax.legend()
#plt.savefig(r"C:\Users\Zach\Documents\RJE\DDD_analysis\Reports\20241008_paper_draft"+"\ddd_comparison.png",dpi=400)
plt.show()