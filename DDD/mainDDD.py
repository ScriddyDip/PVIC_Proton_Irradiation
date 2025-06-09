import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt





Excitation_Incidence = "sunnyside"
Dopant = "Cu"

NIEL_df = df = pd.read_csv(r"C:\Users\Zach\Documents\Python Scripts\DDD\dataNIEL.csv")
df = pd.read_csv(r"C:\Users\Zach\Documents\Python Scripts\DDD\Task 24_AuburnSamples_JVSummary_SampleLevel_forDDD.csv")
df = df.loc[df["Exposure Toggle"] == "expose"]
df = df.loc[df["Excitation Incidence"] == Excitation_Incidence]
df = df.loc[df["CdTe p-Type Dopant"] == Dopant]

sampleIDList = df["ShortID"].unique()
sampleIDList = sampleIDList.tolist()
print(sampleIDList)
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
color = iter(plt.cm.viridis(np.linspace(0,1,n_energy)))



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
    ax.plot(avgdataList[:,0],avgdataList[:,1], linewidth=2, label = label, c = c)
    
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
fig.suptitle(f"Proton Irradiated {Dopant} Doped CdSeTe", fontsize = 18,weight = 'bold')
fig.tight_layout()

# plt.savefig(r"D:\RJE\protRad\Plots\",dpi = 400)
plt.show()





    
# Scaling to Damage Displacement Dose

energyList = [eval(i) for i in energyList]

fig,ax = plt.subplots()
color = iter(plt.cm.viridis(np.linspace(0,1,n_energy)))

for energy in energyList:
    MeV_energy = (1/1000)*energy
    label = f"{MeV_energy} MeV"
    c = next(color)

    NIEL = np.interp(MeV_energy,NIEL_df["Energy (MeV)"],NIEL_df["NIEL (MeV cm2 g-1)"])

    buffer_df = damage_Vs_fluence_curve_dict[energy]
    buffer_df[:,0] = NIEL*buffer_df[:,0]

    ax.plot(buffer_df[:,0],buffer_df[:,1], linewidth=2, label = label, c = c)

    for k in damage_Vs_fluence_curve_dict_individualdata[energy]:
        x = NIEL*k[0]
        y = k[1]
        ax.plot(x,y,marker = ("."),c=c)


ax.set_ylabel("Remaining Efficiency Factor", weight = 'bold', fontsize = 15)
ax.set_xlabel(r"Displacement Damage Dose [$\bf{MeV/g}$]", weight = 'bold', fontsize = 15)
ax.set_xscale("log")
ax.set_ylim(0,1.1)
ax.xaxis.set_tick_params(labelsize = 15)
ax.yaxis.set_tick_params(labelsize = 15)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1])
fig.suptitle(f"Proton Irradiated {Dopant} Doped CdSeTe", fontsize = 18,weight = 'bold')
fig.tight_layout()

# plt.savefig(r"D:\RJE\protRad\Plots\",dpi = 400)
plt.show()