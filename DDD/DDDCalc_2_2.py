import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import re

saveFig = True
titleBool = True
path = r"C:\Users\Zach\Documents\RJE\DDD_analysis\20240708_CdSeTe_Auburn_lambright\Plots\DDDCalc_2\DDDCalc_2_2"

# energylist = ['150','650','1000']
energylist = ['650','1000']

dopant_list = ["Cu","As"]

excitation_incidence_list = ["sunnyside","filmside"]

solar_cell_parameter_list = ["Mean(PCE (%))","Mean(Voc (V))","Mean(Jsc (mA/cm^2))","Mean(Fill Factor (%))"]


scp_uncertainty_list = ["Std Dev(PCE (%))","Std Dev(Voc (V))","Std Dev(Jsc (mA/cm^2))","Std Dev(Fill Factor (%))"]
    

df = pd.read_excel(r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\Task 24_AuburnSamples_JVSummary_SampleLevel__meansAndStdDev_forDDD.xlsx")
NIEL_df = pd.read_csv(r"C:\Users\Zach\Documents\VSCodeFolders\Python Scripts\DDD\dataNIEL.csv")

def uncertaintyFormulaDivMult(transformed_data, p, q, p_std, q_std):
    delta_transformed_data = transformed_data*np.sqrt((p_std/p)**2+(q_std/q)**2)
    return(delta_transformed_data)

def DamageCurveGrabber(df, energy, dopant, excitation_incidence, solar_cell_parameter, scp_uncertainty):

    buffer_df = df.copy()
    buffer_df = buffer_df.loc[(df["Excitation Incidence"] == excitation_incidence) & (df["Proton Acceleration Voltage (kV)"] == energy) & (df["Lightsoak (hours)"] == 0) & (df["CdTe p-Type Dopant"] == dopant)]

    sampleIDList = buffer_df["ShortID"].unique()
    sampleIDList = sampleIDList.tolist()
    damage_curve = []
    for ID in sampleIDList:
        try:
            temp_df = buffer_df.loc[buffer_df["ShortID"] == ID]
            efficiency_before_df = temp_df.loc[temp_df["Proton Exposure Toggle"] == "before"]
            efficiency_after_df = temp_df.loc[temp_df["Proton Exposure Toggle"] == "after"]
            efficiency_before = efficiency_before_df[solar_cell_parameter].iloc[0]
            efficiency_after = efficiency_after_df[solar_cell_parameter].iloc[0]
            stdDev_before = efficiency_before_df[scp_uncertainty].iloc[0]
            stdDev_after = efficiency_after_df[scp_uncertainty].iloc[0]
            remaining_factor = efficiency_after/efficiency_before
            uncertainty_Reff = uncertaintyFormulaDivMult(remaining_factor,p=efficiency_after,q=efficiency_before,p_std=stdDev_after,q_std=stdDev_before)
            fluence = temp_df["Proton Fluence (# cm^-2)"].unique()
            uncertainty_fluence = (0.1)*fluence
            data_point = ([fluence[0]],[remaining_factor],[uncertainty_Reff],[uncertainty_fluence[0]],[ID])
            damage_curve.append(data_point)                

        except:
            print(f"    Error grabbing data for sample {ID}")
            continue
    damage_curve.sort(key=lambda x: x[0])
    return(np.array(damage_curve))

def nielScaling(data_array, energy):
    
    fluence_array = data_array[:,0][:,0]
    fluence_uncertainty = data_array[:,3][:,0]
    MeV_energy = (10**(-3))*eval(energy)
    NIEL = np.interp(MeV_energy,NIEL_df["Energy (MeV)"],NIEL_df["NIEL (MeV cm2 g-1)"])
    NIEL_uncertainty = NIEL*0.1
    NIEL_scaled_array = np.copy(data_array)
    NIEL_scaled_array[:,0][:,0] = NIEL*fluence_array
    for i in range(len(fluence_uncertainty)):
        DDD_uncertainty = uncertaintyFormulaDivMult(transformed_data=NIEL_scaled_array[:,0][i,0],p=data_array[:,0][i,0],q=NIEL,p_std=fluence_uncertainty[i],q_std=NIEL_uncertainty)
        NIEL_scaled_array[:,3][i,0]=DDD_uncertainty

    return(NIEL_scaled_array)

def fitting_function(x,C,Dx):
    remaining_factor = 1-C*np.log10(1+x/Dx)
    return remaining_factor

def reduced_chi_squared(x_arraylike,y_arraylike,y_uncertainty_arraylike,fit_func, fit_params):
    Sigma = 0
    Ndata = len(x_arraylike)
    Nparams = len(fit_params)
    deg_free = Ndata-Nparams
    for index, x in enumerate(x_arraylike):
        numerator = (y_arraylike[index] - fit_func(x,*fit_params))**2
        denominator = (y_uncertainty_arraylike[index])**2
        Sigma = Sigma+(numerator/denominator)
    return(Sigma/deg_free)


#generating all plots for the dataset

for Dopant in dopant_list:
    for Excitation_Incidence in excitation_incidence_list:
        for supreme_index,Solar_cell_parameter in enumerate(solar_cell_parameter_list):
            Scp_uncertainty = scp_uncertainty_list[supreme_index]
            fig,ax = plt.subplots()
            n_energy = len(energylist)

            if n_energy < 3:
                n_energy = n_energy+1
            color = iter(plt.cm.plasma(np.linspace(0,1,n_energy+1)))
            fluence_all_energy_fitlist = []
            if n_energy < 4:
                next(color)
            for Energy in energylist: 
                Energy = eval(Energy)
                data = DamageCurveGrabber(df,energy=Energy,dopant=Dopant,excitation_incidence=Excitation_Incidence,solar_cell_parameter=Solar_cell_parameter,scp_uncertainty=Scp_uncertainty)
                data = np.delete(data,obj = 4,axis=1)
                data = np.array(data,float)
                col = next(color)
                ax.errorbar(x=data[:,0][:,0],y=data[:,1][:,0],yerr=data[:,2][:,0],xerr=data[:,3][:,0], marker = 's',color = col,ecolor="black",capsize=3, label=f"{Energy*(1e-3)} MeV",linestyle='')

                yavgData = []
                xavgData = []
                i=0
                while True:
                    try:
                        if data[:,0][i,0] == data[:,0][i+1,0]:
                            yavgData.append((data[:,1][i,0]+data[:,1][i+1,0])/2)
                            xavgData.append(data[:,0][i,0])
                            i=i+2
                        elif data[:,0][i,0] <= data[:,0][i+1,0]:
                            yavgData.append(data[:,1][i,0])
                            xavgData.append(data[:,0][i,0])
                            i=i+1
                    except:
                        break
                
                ax.plot(xavgData,yavgData,color = col,linewidth=2)
                for i in range(len(data[:,0][:,0])):
                    fluence_all_energy_fitlist.append((data[:,0][i,0],data[:,1][i,0],data[:,2][i,0]))
            fluence_all_energy_fitlist.sort(key=lambda x: x[0])
            fluence_all_energy_fitlist = np.array(fluence_all_energy_fitlist)
            fluence_all_energy_fluence_data = fluence_all_energy_fitlist[:,0]
            fluence_all_energy_Rf_data = fluence_all_energy_fitlist[:,1]
            fluence_all_energy_Rf_uncertainty_data = fluence_all_energy_fitlist[:,2]
            # popt, pcov = curve_fit(fitting_function, fluence_all_energy_fluence_data, fluence_all_energy_Rf_data, sigma=fluence_all_energy_Rf_uncertainty_data)
            # curve_fit_abscissa = np.logspace(np.log10(fluence_all_energy_fluence_data[0]),np.log10(fluence_all_energy_fluence_data[-1]),100)
            # redchisq = reduced_chi_squared(fluence_all_energy_fluence_data, fluence_all_energy_Rf_data, fluence_all_energy_Rf_uncertainty_data,fit_func=fitting_function,fit_params=popt)
            # ax.plot(curve_fit_abscissa,fitting_function(curve_fit_abscissa,*popt), color = "black", label = "\n"+f"C = {popt[0]:.3f}, "+r"$D_x$ = "+f"{popt[1]:.2e} MeV/g" +"\n"+r"$\chi_{red}^2=$"+f"{redchisq}", linewidth=2)

            ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
            ax.set_xlabel(r"Fluence [$\bf{H^+ \text{ions}/cm^2}$]", weight = 'bold', fontsize = 15)
            ax.set_xscale("log")
            ax.set_ylim(0,1.1)
            ax.xaxis.set_tick_params(labelsize = 15)
            ax.yaxis.set_tick_params(labelsize = 15)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], fontsize = 11)
            ax.set_xscale('log')
            if titleBool:
                fig.suptitle(f"{Dopant} {Excitation_Incidence} "+re.sub('[^A-Za-z0-9]+', '', Solar_cell_parameter))
            if saveFig:
                plt.savefig(path+f"\{Dopant}_{Excitation_Incidence}_{energylist[0]}to{energylist[-1]}_"+re.sub('[^A-Za-z0-9]+', '', Solar_cell_parameter)+"_fluence.png",dpi = 400)
                plt.close()


            fig,ax = plt.subplots()
            color = iter(plt.cm.plasma(np.linspace(0,1,n_energy+1)))
            DDD_all_energy_fitlist = []
            markerList = ["o","v","s"]
            if n_energy < 4:
                next(color)
            for index, Energy in enumerate(energylist):
                Energy = eval(Energy) 
                data = DamageCurveGrabber(df,energy=Energy,dopant=Dopant,excitation_incidence=Excitation_Incidence,solar_cell_parameter=Solar_cell_parameter,scp_uncertainty=Scp_uncertainty)
                data = np.delete(data,obj = 4,axis=1)
                data = np.array(data,float)   
                NIEL_data = nielScaling(data,energy=str(Energy))
                col = next(color)
                ddd_data = NIEL_data[:,0][:,0]
                Rf_data = NIEL_data[:,1][:,0]
                Rf_uncertainty_data = NIEL_data[:,2][:,0]
                label = f"{Energy*(1e-3)} MeV"
                ax.errorbar(x=ddd_data,y=Rf_data,yerr=Rf_uncertainty_data,xerr=NIEL_data[:,3][:,0], marker = "s", color = col,ecolor="black",capsize=3, linestyle='', label=label)
                for i in range(len(ddd_data)):
                    DDD_all_energy_fitlist.append((ddd_data[i],Rf_data[i],Rf_uncertainty_data[i]))

            DDD_all_energy_fitlist.sort(key=lambda x: x[0])
            DDD_all_energy_fitlist = np.array(DDD_all_energy_fitlist)
            DDD_all_energy_DDD_data = DDD_all_energy_fitlist[:,0]
            DDD_all_energy_Rf_data = DDD_all_energy_fitlist[:,1]
            DDD_all_energy_Rf_uncertainty_data = DDD_all_energy_fitlist[:,2]
            popt, pcov = curve_fit(fitting_function, DDD_all_energy_DDD_data, DDD_all_energy_Rf_data, sigma=DDD_all_energy_Rf_uncertainty_data)
            curve_fit_abscissa = np.logspace(np.log10(DDD_all_energy_DDD_data[0]),13,100)
            redchisq = reduced_chi_squared(DDD_all_energy_DDD_data, DDD_all_energy_Rf_data, DDD_all_energy_Rf_uncertainty_data,fit_func=fitting_function,fit_params=popt)
            ax.plot(curve_fit_abscissa,fitting_function(curve_fit_abscissa,*popt), color = "black", label = "\n"+f"C = {popt[0]:.3f}, "+r"$D_x$ = "+f"{popt[1]:.2e} MeV/g" +"\n"+r"$\chi_{red}^2=$"+f"{redchisq:.5f}", linewidth=2)

            ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
            ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
            ax.set_xscale("log")
            ax.set_ylim(0,1.1)
            ax.xaxis.set_tick_params(labelsize = 15)
            ax.yaxis.set_tick_params(labelsize = 15)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], fontsize = 11)
            ax.set_xscale('log')

            if titleBool:
                fig.suptitle(f"{Dopant} {Excitation_Incidence} "+re.sub('[^A-Za-z0-9]+', '', Solar_cell_parameter))
                # props = dict(boxstyle='round', facecolor='white', alpha=0.5)
                # fig.text(0.175,0.425, r"$\frac{\eta_{post}}{\eta_{pre}}=1-C*log(1+\frac{DDD}{D_x})$", fontsize = 15,bbox=props)
            if saveFig:
                    plt.savefig(path+f"\{Dopant}_{Excitation_Incidence}_{energylist[0]}to{energylist[-1]}_"+re.sub('[^A-Za-z0-9]+', '', Solar_cell_parameter)+"_DDD.png",dpi = 400)
                    plt.close()
            print(f"\{Dopant}_{Excitation_Incidence}_{energylist[0]}to{energylist[-1]}_"+re.sub('[^A-Za-z0-9]+', '', Solar_cell_parameter))
