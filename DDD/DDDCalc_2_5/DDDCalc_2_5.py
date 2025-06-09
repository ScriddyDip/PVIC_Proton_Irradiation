import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import re


def damage_curve_correlator(energy, dopant, uncorrelated_dataset_df, curve_storage_dict=None, excitation_incidence="sunnyside"):
    # Filter the uncorrelated_dataset_df based on the conditions
    correlated_dataset_df = uncorrelated_dataset_df[
        (uncorrelated_dataset_df['Proton Acceleration Voltage (kV)'] == energy) &
        (uncorrelated_dataset_df['CdTe p-Type Dopant'] == dopant) &
        (uncorrelated_dataset_df['Excitation Incidence'] == excitation_incidence)
    ].copy()  # Using .copy() to avoid modifying the original dataframe
    
    #Sorting degradation curve
    correlated_dataset_df = correlated_dataset_df.sort_values(by='Proton Fluence (cm^-2)')

    if curve_storage_dict != None:
        # Create a key for the dictionary
        curve_key = f"{dopant}_{energy}"
        
        # Store the filtered dataframe in the dictionary
        curve_storage_dict[curve_key] = correlated_dataset_df
    
    # Return the filtered dataframe
    return correlated_dataset_df

def niel_scaling(energy, correlated_dataset_df, niel_curve_df):
    # Convert energy to MeV
    MeV_energy = (10**(-3)) * float(energy)
    
    # Interpolate NIEL for the given energy
    niel = np.interp(MeV_energy, niel_curve_df["Energy (MeV)"], niel_curve_df["NIEL (MeV cm2 g-1)"])
    
    # Create a copy of the input dataframe to avoid modifying the original
    niel_scaled_df = correlated_dataset_df.copy()
    
    # Scale the fluence column by the interpolated NIEL
    niel_scaled_df['Proton Fluence (cm^-2)'] = niel_scaled_df['Proton Fluence (cm^-2)'] * niel
    
    # Return the modified dataframe
    return niel_scaled_df

def fitting_function(x,C,Dx):
    remaining_factor = 1-C*np.log10(1+x/Dx)
    return remaining_factor


sample_data_df = pd.read_csv(r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\JVByControls_forDDD.csv")
absorber_niel_df = pd.read_csv(r"C:\Users\slamb\OneDrive\Desktop\Displacement Damage Dose\DDD\DDDCalc_2_5\cdte_proton_niel.csv")
dopant_types = ["Copper", "Arsenic"]

parameter_list = ["Median(PCE/(Control PCE))","Median(Voc/(Control Voc))","Median(Jsc/(Control Jsc))","Median(FF/(Control FF))"]

#energies to loop over
energy_list = [650,1000]
n_energy = len(energy_list)

#setting a "switch" in order to change set behavior later
switch = 0
if n_energy<3:
    n_energy=n_energy+1
    switch=1
    
#initialize figures
fig,axs = plt.subplots(2,2,figsize=(12,10))

#Looping over all the parameters in order to plot them in a 2x2 figure
for n, p in enumerate(parameter_list):
    #Setting variable for assignment of subplot location
    if n<=1:
        figure_index = 0
    elif n<=3:
        figure_index = 1
        n=n-2
        
    #initializing list to store data values to be fitted
    ddd_cumulative_fit_list = []
    
    #set up "iterator" to step through colors
    color = iter(plt.cm.plasma(np.linspace(0,1,n_energy+1)))
    
    #skip a color if switch is equal to 1
    if switch==1:
        dummy = next(color)
        
    #regular expression to extract label
    label_p = re.search(r'\((.*?)\)', p).group(1).split('/')[0]

# for loop to loop over energies for fitting degradation curves
    for energy in energy_list:
        correlated_data = damage_curve_correlator(
            energy=energy,
            dopant=dopant,
            uncorrelated_dataset_df=sample_data_df
        )

        scaled_data = niel_scaling(
            energy=energy,
            correlated_dataset_df=correlated_data,
            niel_curve_df=absorber_niel_df
            )
        
        ddd_list = scaled_data["Proton Fluence (cm^-2)"].values
        remaining_factor_list = scaled_data[p].values

        c = next(color)
        axs[n,figure_index].scatter(ddd_list,remaining_factor_list,color=c,label=f"{energy*(1e-3)} MeV")
        for i,ddd in enumerate(ddd_list):
            ddd_cumulative_fit_list.append((ddd_list[i],remaining_factor_list[i]))

#Sort list of tuples according to first position of tuple
    ddd_cumulative_fit_list.sort(key=lambda x: x[0])
    ddd_cumulative_fit_list = np.array(ddd_cumulative_fit_list)
    ddd_cumulative_list = ddd_cumulative_fit_list[:,0]
    rf_cumulative_list = ddd_cumulative_fit_list[:,1]

#fitting of data points and calculation of r-squared
    popt, pcov = curve_fit(fitting_function, ddd_cumulative_list, rf_cumulative_list)
    curve_fit_abscissa = np.logspace(np.log10(ddd_cumulative_list[0]),13,100)
    fit_y_values = fitting_function(curve_fit_abscissa,*popt)
    y_pred = fitting_function(ddd_cumulative_list,*popt)
    sst = np.sum((rf_cumulative_list - np.mean(rf_cumulative_list))**2)
    ssr = np.sum((rf_cumulative_list - y_pred)**2)
    r_squared = 1 - (ssr / sst)
    
    axs[n,figure_index].plot(curve_fit_abscissa,fit_y_values, color = "black", label = "\n"+f"C = {popt[0]:.3f}, "+r"$D_x$ = "+f"{popt[1]:.2e} MeV/g: " +r"$R^2=$" +f"{r_squared:.4f}", linewidth=2)
    axs[n,figure_index].set_ylabel(f"Relative Change in {label_p}", weight = 'bold', fontsize = 15)
    axs[n,figure_index].set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
    axs[n,figure_index].set_xscale("log")
    axs[n,figure_index].set_ylim(0,1.1)
    axs[n,figure_index].xaxis.set_tick_params(labelsize = 15)
    axs[n,figure_index].yaxis.set_tick_params(labelsize = 15)
    handles, labels = axs[n,figure_index].get_legend_handles_labels()
    axs[n,figure_index].legend(handles[::-1], labels[::-1], fontsize = 11)
    axs[n,figure_index].legend()
# plt.savefig(r"C:\Users\Zach\Documents\RJE\DDD_analysis\Reports\20241008_paper_draft"+f"\{dopant}_fit.png",dpi=400)
plt.show()

