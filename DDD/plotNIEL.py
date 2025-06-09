import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

energy = 1 #MeV
fluence = 1e11 #ions/cm^2

C = .633
Dx = 2.32e11
# df = pd.read_csv(r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\dataNIEL.csv")
df2 = pd.read_csv(r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\dataNIEL2.csv")

def fitting_function(x,C,Dx):
    remaining_factor = 1-C*np.log10(1+x/Dx)
    return remaining_factor

fig, ax = plt.subplots()
# ax.plot(df2["Energy (MeV)"], df2["NIEL (MeV cm2 g-1)"], linewidth=2, label = 2)
NIEL = np.interp(energy,df2["Energy (MeV)"],df2["NIEL (MeV cm2 g-1)"])
DDD = NIEL*fluence
print(f"Energy: {energy} MeV, Fluence: {fluence:.2e} ions/cm^2, DDD: {DDD:.3e} MeV/g, EOL: {fitting_function(DDD,C,Dx):.2%}")
label = f"({energy}, {fluence:.2e}), EOL: {fitting_function(DDD,C,Dx):.2%}"
# color = iter(plt.cm.plasma(np.linspace(0,1,n_energy+1)))
# for energy in energyList:
#     c = next(color)
#     NIEL = np.interp(energy,df["Energy (MeV)"],df["NIEL (MeV cm2 g-1)"])
#     ax.plot(energy,NIEL,color=c,marker='o', label=f'{energy} MeV, {NIEL} '+r'$MeV cm^2 g^{-1}$')
#     print(f"{energy}: {NIEL*1e14:.3E}")

curve_fit_abscissa = np.logspace(9,13,100)
ax.plot(curve_fit_abscissa,fitting_function(curve_fit_abscissa,C,Dx))
ax.scatter(DDD,fitting_function(DDD,.633,2.32e11), label=label)

ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
ax.set_xscale("log")
ax.set_ylim(0,1.1)
ax.xaxis.set_tick_params(labelsize = 15)
ax.yaxis.set_tick_params(labelsize = 15)
ax.set_xscale('log')
ax.legend()


# # plt.savefig(r"C:\Users\Zach\Downloads\NIEL_Curve.png",dpi = 400)
plt.show()