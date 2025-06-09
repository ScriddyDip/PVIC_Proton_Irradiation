import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def fitting_function(x,C,Dx,BoL):
    remaining_factor = (1-C*np.log10(1+x/Dx))*BoL
    return remaining_factor

def nonlinear_electron_dependence(differential_spectrum_df, absorber_niel_df, quality_factor):
    # Extract energy values
    x = differential_spectrum_df["Energy\n(MeV)"].copy()
    
    # Interpolate NIEL values for the corresponding energies
    niel = np.interp(x, absorber_niel_df["Energy (MeV)"], absorber_niel_df["NIEL (MeV cm2 g-1)"])
    niel = niel**quality_factor
    niel_1mev_electron = np.interp(1, absorber_niel_df["Energy (MeV)"], absorber_niel_df["NIEL (MeV cm2 g-1)"])
    # Extract the fluence values
    fluence = differential_spectrum_df["Total mission\nfluence\n(/cm2/MeV)"]
    
    # Calculate ddd = fluence * NIEL
    ddd = fluence * niel
    nonlinear_coefficient = 1/(niel_1mev_electron**(quality_factor-1))

    # Perform trapezoidal integration over energy (x)
    cumulative_ddd = np.trapz(ddd, x)
    cumulative_ddd = cumulative_ddd*nonlinear_coefficient

    return cumulative_ddd

 #for protons or any particle with linear NIEL dependence 
def spectrum_ddd(differential_spectrum_df, absorber_niel_df):
    # Extract energy values
    x = differential_spectrum_df["Energy\n(MeV)"].copy()
    
    # Interpolate NIEL values for the corresponding energies
    niel = np.interp(x, absorber_niel_df["Energy (MeV)"], absorber_niel_df["NIEL (MeV cm2 g-1)"])
    
    # Extract the fluence values
    fluence = differential_spectrum_df["Total mission\nfluence\n(/cm2/MeV)"]
    
    # Calculate ddd = fluence * NIEL
    ddd = fluence * niel

    # Perform trapezoidal integration over energy (x)
    cumulative_ddd = np.trapz(ddd, x)

    return cumulative_ddd

def sum_ddd_binned(bin_fluence_df,bin_energy_upperbound_df,bin_energy_lowerbound_df,bin_mean_energy,absorber_niel_df):
    ddd = 0

    for i,fluence in enumerate(bin_fluence_df):
        #calculate bin width and convert from keV to MeV
        bin_width = 10**(-3)*float(bin_energy_upperbound_df.iloc[i]) - 10**(-3)*float(bin_energy_lowerbound_df.iloc[i])
        #convert from keV to MeV
        e_avg = 10**(-3)*float(bin_mean_energy.iloc[i])
        niel = np.interp(e_avg, absorber_niel_df["Energy (MeV)"], absorber_niel_df["NIEL (MeV cm2 g-1)"])
        ddd_step = fluence*niel*bin_width
        ddd = ddd + ddd_step
    
    return(ddd)

CdTe_p_niel_df = pd.read_csv(r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\cdte_proton_niel.csv")
CdTe_e_niel_df = pd.read_csv(r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\cdte_electron_niel.csv")
GaAs_p_niel_df = pd.read_csv(r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\gaas_niel.csv")
proton_differential_spectrum_df = pd.read_csv(r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\differential_proton_leo - Sheet1.csv")
electron_differential_spectrum_df = pd.read_csv(r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\differential_electron_leo - Sheet1.csv")

fit_coef_dict = {
    "CdSeTe:Cu" : (.632,2.45e11,1),
    "CdSeTe:As" : (.521,7.72e10,1),
    "Azur 3G28" : (.249,1.73e9,1),
    "Azur 4G32" : (.359,6.94e9,1),
    "CESI CTJ-LC" : (.354,4.81e9,1)
}

fig,ax = plt.subplots()

for key in fit_coef_dict.keys():
    C = fit_coef_dict[key][0]
    Dx = fit_coef_dict[key][1]
    BoL = fit_coef_dict[key][2]
    if key.startswith("CdSeTe"):
        p_niel_df = CdTe_p_niel_df
        e_niel_df = CdTe_e_niel_df
        linestyle = "solid"
    else:
        p_niel_df = GaAs_p_niel_df
        linestyle = "dashed"
        
    ddd_p = spectrum_ddd(proton_differential_spectrum_df,p_niel_df)
    ddd_e = nonlinear_electron_dependence(electron_differential_spectrum_df,e_niel_df, 1)
    print(f"Material: {key} EOL: {fitting_function(ddd_p,C,Dx)}\nMission Proton DDD: {ddd_p:.2e}\nMission Electron DDD: {ddd_e:.2e}")
    orbit_step_array = np.linspace(0,400,101)
    mission_length_array = 10*orbit_step_array #duration of orbit in days * number of orbit segments plotted
    ax.plot(mission_length_array,fitting_function(ddd_p*orbit_step_array,C,Dx,BoL), label = key, linewidth=2, linestyle = linestyle)

ax.set_ylabel("PCE Remaining Factor", weight = 'bold', fontsize = 15)
ax.set_xlabel("Days in Orbit", weight = 'bold', fontsize = 15)

ax.xaxis.set_tick_params(labelsize = 12)
ax.yaxis.set_tick_params(labelsize = 12)
fig.tight_layout()
ax.legend()
# plt.savefig(r"C:\Users\Zach\Documents\RJE\DDD_analysis\Reports\20241008_paper_draft"+"\days_in_orbit.png",dpi=400)
plt.show()

# Read in Shielded spectra

shielded_proton_spectrum_path_list = [
r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\30day5000km_0um.csv",
r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\30day5000km_25um.csv",
r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\30day5000km_50um.csv",
r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\30day5000km_100um.csv",
r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\30day5000km_150um.csv",
r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\30day5000km_300um.csv",
]

# [
# r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\DDDCalc_2_5\30day500km_0um.csv",
# r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\DDDCalc_2_5\30day500km_25um.csv",
# r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\DDDCalc_2_5\30day500km_50um.csv",
# r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\DDDCalc_2_5\30day500km_100um.csv",
# r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\DDDCalc_2_5\30day500km_150um.csv",
# r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\DDDCalc_2_5\30day500km_300um.csv",
# ]






fig,ax = plt.subplots()
fig1,ax1 = plt.subplots()
for key in fit_coef_dict.keys():
    print(key)
    C = fit_coef_dict[key][0]
    Dx = fit_coef_dict[key][1]
    BoL = fit_coef_dict[key][2]
    if key.startswith("CdSeTe"):
        p_niel_df = CdTe_p_niel_df
        linestyle = "solid"
    else:
        p_niel_df = GaAs_p_niel_df
        linestyle = "dashed"
        
    shielding_eol_list = []
    for path in shielded_proton_spectrum_path_list:
        shielded_proton_spectrum_df = pd.read_csv(path, skiprows=8)
        e_lb_df = shielded_proton_spectrum_df["Lower edge of energy bin [keV]"]
        e_ub_df = shielded_proton_spectrum_df[" Upper edge of energy bin [keV]"]
        e_meanb_df = shielded_proton_spectrum_df[" Mean energy of bin [keV]"]
        fluenceb_df = shielded_proton_spectrum_df[" Fluence [particles/cm2/bin]"]
        fluenceb_error_df = shielded_proton_spectrum_df[" Error in fluence/flux [particles/cm2/bin]"]
        shielding = path.split("_")[-1]
        shielding = shielding.split(".")[0]
        ddd_p = sum_ddd_binned(bin_energy_lowerbound_df=e_lb_df,
                        bin_energy_upperbound_df=e_ub_df,
                        bin_fluence_df=fluenceb_df,
                        bin_mean_energy=e_meanb_df,
                        absorber_niel_df=p_niel_df
                        )
        annual_ddd = (365/30)*ddd_p
        eol = fitting_function(annual_ddd,C,Dx)
        shielding_val = eval(shielding.split("u")[0])
        shielding_eol_list.append((shielding_val,eol))
        
        if shielding_val == 0:
            orbit_step_array = np.linspace(0,75,101)
            mission_length_array = 30*orbit_step_array #duration of orbit in days * number of orbit segments plotted
            ax1.plot(mission_length_array,fitting_function(ddd_p*orbit_step_array,C,Dx), label = key, linewidth=2.5, linestyle = linestyle)


    shielding_eol_list = np.array(shielding_eol_list)
    print(shielding_eol_list[:,0])
    ax.plot(shielding_eol_list[:,0],shielding_eol_list[:,1], label = key, linestyle = linestyle, linewidth = 2.5, marker = "o")
ax1.set_ylabel("PCE Remaining Factor", weight = 'bold', fontsize = 15)
ax1.set_xlabel("Days in Orbit", weight = 'bold', fontsize = 15)
ax1.set_ylim(-0.1,1.1)

ax1.xaxis.set_tick_params(labelsize = 12)
ax1.yaxis.set_tick_params(labelsize = 12)
fig1.tight_layout()
ax1.legend()
ax.set_ylabel("PCE Remaining Factor", weight = 'bold', fontsize = 15)
ax.set_xlabel("Shielding Thickness [um]", weight = 'bold', fontsize = 15)
ax.xaxis.set_tick_params(labelsize = 12)
ax.yaxis.set_tick_params(labelsize = 12)
fig.tight_layout()
ax.legend()
fig1.savefig(r"C:\Users\slamb\OneDrive\Documents\Publications\2025\EOL Predictions for Industry-grade CdSeTe PV Devices in LEO\Images\eol_vs_days_in_orbit_5000km_no_shielding.png",dpi=400)
fig.savefig(r"C:\Users\slamb\OneDrive\Documents\Publications\2025\EOL Predictions for Industry-grade CdSeTe PV Devices in LEO\Images\eol_vs_shielding_one_year_5000km.png",dpi=400)

plt.show()
