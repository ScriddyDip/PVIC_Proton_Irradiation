import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_average(numbers):
    return sum(numbers) / len(numbers)

def calculate_arc_length(c,r):
    central_angle = 2*np.arcsin(c/(2*r))
    arc_length = r*central_angle
    return arc_length

def calculate_circumference(r):
    return 2*np.pi*r

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
    "CdSeTe:Cu" : (.656,2.76e11,16.8),
    "CdSeTe:As" : (.551,9.82e10, 16.8),
    "Azur 3G28" : (.249,1.73e9, 28),
    "Azur 4G32" : (.359,6.94e9, 32),
    "CESI CTJ-LC" : (.354,4.81e9,28)
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
    print(f"Material: {key} EOL: {fitting_function(ddd_p,C,Dx,BoL)}\nMission Proton DDD: {ddd_p:.2e}\nMission Electron DDD: {ddd_e:.2e}")
    orbit_step_array = np.linspace(0,400,101)
    mission_length_array = 10*orbit_step_array #duration of orbit in days * number of orbit segments plotted
    ax.plot(mission_length_array,fitting_function(ddd_p*orbit_step_array,C,Dx,BoL), label = key, linewidth=2, linestyle = linestyle)

ax.set_ylabel("Absolute PCE (%)", weight = 'bold', fontsize = 15)
ax.set_xlabel("Days in Orbit", weight = 'bold', fontsize = 15)

ax.xaxis.set_tick_params(labelsize = 12)
ax.yaxis.set_tick_params(labelsize = 12)
ax.set_ylim(-0.1,33)
fig.tight_layout()
ax.legend()
#fig.savefig(r"C:\Users\slamb\OneDrive\Documents\Publications\2025\EOL Predictions for Industry-grade CdSeTe PV Devices in LEO\Images\Abs-PCE_vs_days_in_orbit_5000km_150um_shielding.png",dpi=400)
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
fig2,ax2 = plt.subplots()
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
        eol = fitting_function(annual_ddd,C,Dx,BoL)
        shielding_val = eval(shielding.split("u")[0])
        shielding_eol_list.append((shielding_val,eol))
        
        if shielding_val == 0:
            orbit_step_array = np.linspace(0,75,76)
            time_step_days = 30
            mission_length_array = time_step_days*orbit_step_array #duration of orbit in days * number of orbit segments plotted
            eol_by_time_in_orbit = fitting_function(ddd_p*orbit_step_array,C,Dx,BoL)
            ax1.plot(mission_length_array,eol_by_time_in_orbit, label = key, linewidth=2.5, linestyle = linestyle)
            
            e_yield_by_time_in_orbit = np.zeros(len(orbit_step_array), dtype=np.float64)
            nameplate = 1
            earth_radius = 6378
            orbit_height = 5000
            orbit_radius = earth_radius + orbit_height
            shadow_arc_length = calculate_arc_length(2*earth_radius,orbit_radius)
            orbit_circumference = calculate_circumference(orbit_radius)
            current_e_yield = 0
            for time_interval in range(len(eol_by_time_in_orbit)):
                if time_interval == 0:
                    current_e_yield = 0
                    e_yield_by_time_in_orbit[time_interval] = 0
                else:
                    current_e_yield = current_e_yield + (calculate_average(eol_by_time_in_orbit[time_interval - 1 : time_interval + 1]))
                    e_yield_by_time_in_orbit[time_interval] = ((current_e_yield)*(24)*(time_step_days)*(nameplate)*((orbit_circumference - shadow_arc_length)/orbit_circumference))/1000
                print((time_step_days*orbit_step_array[time_interval]),e_yield_by_time_in_orbit[time_interval])
            
            ax2.plot(mission_length_array,e_yield_by_time_in_orbit, label = key, linewidth=2.5, linestyle = linestyle)

    shielding_eol_list = np.array(shielding_eol_list)
    print(shielding_eol_list[:,0])
    ax.plot(shielding_eol_list[:,0],shielding_eol_list[:,1], label = key, linestyle = linestyle, linewidth = 2.5, marker = "o")

ax.set_ylabel("Absolute PCE (%)", weight = 'bold', fontsize = 15)
ax.set_xlabel("Shielding Thickness [um]", weight = 'bold', fontsize = 15)
ax.xaxis.set_tick_params(labelsize = 12)
ax.yaxis.set_tick_params(labelsize = 12)
ax.set_ylim(-0.1,33)
fig.tight_layout()
# ax.legend()

ax1.set_ylabel("Absolute PCE (%)", weight = 'bold', fontsize = 15)
ax1.set_xlabel("Days in Orbit", weight = 'bold', fontsize = 15)
ax1.xaxis.set_tick_params(labelsize = 12)
ax1.yaxis.set_tick_params(labelsize = 12)
ax1.set_ylim(-0.1,33)
fig1.tight_layout()
# ax1.legend()

ax2.set_ylabel("Total Energy Yield (MW*h)", weight = 'bold', fontsize = 15)
ax2.set_xlabel("Days in Orbit", weight = 'bold', fontsize = 15)
ax2.xaxis.set_tick_params(labelsize = 12)
ax2.yaxis.set_tick_params(labelsize = 12)
ax2.set_ylim(-0.1,45)
fig2.tight_layout()

#fig.savefig(r"C:\Users\slamb\OneDrive\Documents\Publications\2025\EOL Predictions for Industry-grade CdSeTe PV Devices in LEO\Images\Abs_PCE_vs_shielding_one_year_5000km.png",dpi=600)
fig1.savefig(r"C:\Users\slamb\OneDrive\Documents\Publications\2025\EOL Predictions for Industry-grade CdSeTe PV Devices in LEO\Images\Abs_PCE_vs_days_in_orbit_5000km_0um_shielding.png",dpi=600)
#fig2.savefig(r"C:\Users\slamb\OneDrive\Documents\Publications\2025\EOL Predictions for Industry-grade CdSeTe PV Devices in LEO\Images\total_e_yield_vs_days_in_orbit_5000km_150um_shielding.png",dpi=600)


plt.show()
