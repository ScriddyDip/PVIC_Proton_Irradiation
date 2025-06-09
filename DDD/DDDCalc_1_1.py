import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt



plotSavingBool = False
plotTitlesEqsBool = True
no100Toggle = True


Excitation_Incidence = "sunnyside"
DopantList = ["Cu","As"]


def fitting_function(x,C,Dx):
    remaining_factor = 1-C*np.log(1+x/Dx)
    return remaining_factor

fitDict = {}

for Dopant in DopantList:

    NIEL_df = df = pd.read_csv(r"C:\Users\Zach\Documents\Python Scripts\DDD\dataNIEL.csv")
    df = pd.read_csv(r"C:\Users\Zach\Documents\Python Scripts\DDD\Task 24_AuburnSamples_JVSummary_SampleLevel_forDDD.csv")
    df = df.loc[df["Exposure Toggle"] == "expose"]
    df = df.loc[df["Excitation Incidence"] == Excitation_Incidence]
    df = df.loc[df["CdTe p-Type Dopant"] == Dopant]

    sampleIDList = df["ShortID"].unique()
    sampleIDList = sampleIDList.tolist()

    remaining_factor_Dict = {}

    for ID in sampleIDList:
        temp_df = df.loc[df["ShortID"] == ID]
        efficiency_before = temp_df.loc[temp_df["Proton Exposure Toggle"] == "before"]
        efficiency_after = temp_df.loc[temp_df["Proton Exposure Toggle"] == "after"]
        
        try:
            efficiency_before = efficiency_before["Median(Efficiency (%))"].iloc[0]
            efficiency_after = efficiency_after["Median(Efficiency (%))"].iloc[0]
        except:
            print(f"Error grabbing data for sample {ID}")
            continue

        remaining_factor = efficiency_after/efficiency_before
        remaining_factor_Dict[ID] = remaining_factor

    energyList = df["Proton Acceleration Voltage (kV)"].unique()
    energyList = energyList.tolist()
    energyList = [eval(i) for i in energyList]
    energyList.sort()
    energyList = [str(i) for i in energyList]

    fluenceList = df["Fluence (# cm^-2)"].unique()
    fluenceList = fluenceList.tolist()
    n_fluences = len(fluenceList)
    energyDict = {}

    for energy in energyList:
        temp_df = df.loc[df["Proton Acceleration Voltage (kV)"] == energy]
        temp_list = []

        for fluence in fluenceList:
            temp_df_2 = temp_df.loc[temp_df["Fluence (# cm^-2)"] == fluence]
            ID_list = temp_df_2["ShortID"].unique()
            ID_list = ID_list.tolist()

            temp_list.append([fluence,ID_list])

        energyDict[energy] = temp_list

    #Plotting

    n_energy = len(energyList)
    color = iter(plt.cm.plasma(np.linspace(0,1,n_energy+1)))



    fig, ax = plt.subplots()

    damage_Vs_fluence_curve_dict = {}
    damage_Vs_fluence_curve_dict_individualdata = {}

    for energy in energyList:
        avgdataList = []
        individualdataList = []
        c = next(color)
        label = f"{energy} keV"

        for i in range(len(energyDict[energy])):
            
            fluence = energyDict[energy][i][0]
            

            avg = 0

            n_samples = len(energyDict[energy][i][1])

            for j in range(n_samples):
                sampleID = energyDict[energy][i][1][j]
                try:
                    remaining_factor = remaining_factor_Dict[sampleID]
                    avg = avg + remaining_factor
                    individualdataList.append((eval(fluence),remaining_factor))
                except:
                    print(f"There was an error getting the data for {sampleID}")
                    if n_samples == 1:
                        avg = 0
                        continue
                    else:
                        n_samples = n_samples - 1
                        continue
            
            avg = avg/n_samples

            avgdataList.append((eval(fluence),avg)) 
        

        avgdataList.sort()
        avgdataList = np.array(avgdataList)
        damage_Vs_fluence_curve_dict[eval(energy)] = avgdataList
        damage_Vs_fluence_curve_dict_individualdata[eval(energy)] = individualdataList
        ax.plot(avgdataList[:,0],avgdataList[:,1], linewidth=2,marker='o', label = label, c = c)
        
        for k in individualdataList:
            ax.plot(*k,marker = ("."),c=c)
                
                
    ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
    ax.set_xlabel(r"Fluence [$\bf{Protons/cm^2}$]", weight = 'bold', fontsize = 15)
    ax.set_xscale("log")
    # ax.set_ylim(0,1.1)
    ax.xaxis.set_tick_params(labelsize = 15)
    ax.yaxis.set_tick_params(labelsize = 15)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])
    if plotTitlesEqsBool:
        fig.suptitle(r"$\bf{H^+}$ Irradiated "+f"{Dopant} Doped CdSeTe {Excitation_Incidence}", fontsize = 18,weight = 'bold')
    fig.tight_layout()

    if plotSavingBool:
        plt.savefig(r"D:\RJE\protRad\Plots"+f"\{Dopant}_{Excitation_Incidence}_CdSeTe_remainingEfficiencyVSfluence.png",dpi = 400)






        
    # Scaling to Damage Displacement Dose

    energyList = [eval(i) for i in energyList]

    fig,ax = plt.subplots()
    color = iter(plt.cm.plasma(np.linspace(0,1,n_energy+1)))

    for energy in energyList:
        MeV_energy = (1/1000)*energy
        c = next(color)

        NIEL = np.interp(MeV_energy,NIEL_df["Energy (MeV)"],NIEL_df["NIEL (MeV cm2 g-1)"])
        # OPERATES THE DDD CALCULATION ON THE DICTIONARY FOR THE FUTURE    
        buffer_df = damage_Vs_fluence_curve_dict[energy]
        buffer_df[:,0] = NIEL*buffer_df[:,0]

        ax.scatter(buffer_df[:,0],buffer_df[:,1], marker = "o", label = "_label_", color = c)

        # for k in damage_Vs_fluence_curve_dict_individualdata[energy]:
        #     x = NIEL*k[0]
        #     y = k[1]
        #     ax.plot(x,y,marker = ("."),c=c)

        popt, pcov = curve_fit(fitting_function, buffer_df[:,0], buffer_df[:,1])
        label = f"{MeV_energy} MeV: C = {popt[0]:.3f}, "+r"$D_x$ = "+f"{popt[1]:.2e} MeV/g"

        # xarray = np.linspace(buffer_df[0,0],buffer_df[-1,0],100)
        xarray = np.logspace(7,15,1000)

        ax.plot(xarray,fitting_function(xarray,*popt), linewidth=2, label = label, color = c)


    ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
    ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
    ax.set_xscale("log")
    ax.set_ylim(0,1.1)
    ax.xaxis.set_tick_params(labelsize = 15)
    ax.yaxis.set_tick_params(labelsize = 15)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], fontsize = 11)

    if plotTitlesEqsBool:    
        fig.suptitle(r"$\bf{H^+}$ Irradiated "+f"{Dopant} Doped CdSeTe {Excitation_Incidence}", fontsize = 18,weight = 'bold')
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        fig.text(0.2,0.45, r"$y=1-C*log(1+\frac{x}{D_x})$", fontsize = 15,bbox=props)
    fig.tight_layout()

    if plotSavingBool:
        plt.savefig(r"D:\RJE\protRad\Plots"+f"\{Dopant}_{Excitation_Incidence}_CdSeTe_damageRatio_vs_DDD.png",dpi = 400)


    fig,ax = plt.subplots()



    if no100Toggle:
        no100List = []
        for energy in energyList[1:]:
            MeV_energy = (1/1000)*energy
            buffer_no100 = damage_Vs_fluence_curve_dict[energy]
            for dataPoint in buffer_no100:
                no100List.append(dataPoint)
                
        dataArray = np.zeros((2,len(no100List)))
        for i in range(len(no100List)):
            dataArray[0,i] = no100List[i][0]
            dataArray[1,i] = no100List[i][1]

        dataArray = dataArray[:,dataArray[0,:].argsort()]
        dataArray = np.transpose(dataArray)
        ax.scatter(dataArray[:,0],dataArray[:,1], marker = "o", label = "_label_", color = "black")

        popt, pcov = curve_fit(fitting_function,dataArray[:,0],dataArray[:,1])
        label = f"C = {popt[0]:.3f}, "+r"$D_x$ = "+f"{popt[1]:.2e} MeV/g"

        # xarray = np.linspace(buffer_df[0,0],buffer_df[-1,0],100)
        xarray = np.logspace(7,15,1000)

        fitDict[Dopant] = popt
        ax.plot(xarray,fitting_function(xarray,popt[0],popt[1]), linewidth=2, label = label, color = "black")


        ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
        ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
        ax.set_xscale("log")
        ax.set_ylim(0,1.1)
        ax.xaxis.set_tick_params(labelsize = 15)
        ax.yaxis.set_tick_params(labelsize = 15)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], fontsize = 11)

        if plotTitlesEqsBool:    
            fig.suptitle(r"$\bf{H^+}$ Irradiated "+f"{Dopant} Doped CdSeTe {Excitation_Incidence}", fontsize = 18,weight = 'bold')
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            fig.text(0.2,0.45, r"$y=1-C*log(1+\frac{x}{D_x})$", fontsize = 15,bbox=props)
        fig.tight_layout()

        if plotSavingBool:
            plt.savefig(r"D:\RJE\protRad\Plots"+f"\{Dopant}_{Excitation_Incidence}_CdSeTe_damageRatio_vs_DDD_no100kev.png",dpi = 400)


fig, ax = plt.subplots()

for Dopant in DopantList:

    xarray = np.logspace(7,15,1000)
    ax.plot(xarray,fitting_function(xarray,fitDict[Dopant][0],fitDict[Dopant][1]), linewidth=2, label = Dopant)
    ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
    ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
    ax.set_xscale("log")
    ax.set_ylim(0,1.1)
    ax.xaxis.set_tick_params(labelsize = 15)
    ax.yaxis.set_tick_params(labelsize = 15)
    ax.legend(fontsize = 11)

    if plotTitlesEqsBool:    
        fig.suptitle(r"$\bf{H^+}$ Irradiated "+f"{Excitation_Incidence}_CdSeTe_DDDvsdopant.png", fontsize = 18,weight = 'bold')
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        fig.text(0.2,0.45, r"$y=1-C*log(1+\frac{x}{D_x})$", fontsize = 15,bbox=props)
    fig.tight_layout()

    if plotSavingBool:
        plt.savefig(r"D:\RJE\protRad\Plots"+f"\{Excitation_Incidence}_CdSeTe_DDDvsdopant.png",dpi = 400)

plt.show()