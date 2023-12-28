import numpy as np
from matplotlib import pyplot as plt
import statistics

# Biontech_______________________________________________________________________________________________________________________
# placebo arm
biontech_placebo_path = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Biontech.txt"
file = open(biontech_placebo_path, 'r')
Lines = file.readlines()
k_inf_placebo_biontech = []

for line in Lines:
    k_inf_placebo_biontech.append(float(line.replace('\n', '')))

# vaccination arm
biontech_vacc_path = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Biontech_Vacc.txt"
file = open(biontech_vacc_path, 'r')
Lines = file.readlines()
k_inf_vacc_biontech = []

for line in Lines:
    k_inf_vacc_biontech.append(float(line.replace('\n', '')))

efficiacy_biontech = []

for placebo, vacc in zip(k_inf_placebo_biontech, k_inf_vacc_biontech):
    efficiacy_biontech.append(1 - (vacc/placebo))




mean_efficiacy_biotech = round(statistics.mean(efficiacy_biontech),4)
stdev_efficiacy_biontech = round(statistics.stdev(efficiacy_biontech),4)

fig, ax = plt.subplots()
plt.hist(efficiacy_biontech, bins=30)
plt.grid()
plt.title("Efficiacy Biontech")
plt.xlabel("Efficiacy")
plt.ylabel("Frequency")
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = "Uncertainity (Mean +- standard deviation):\n{} +- {} "
ax.text(0.05, 0.95, textstr.format(mean_efficiacy_biotech, stdev_efficiacy_biontech), transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
plt.show()



# Moderna _________________________________________________________________________________________________________________________

# placebo arm
moderna_placebo_path = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Moderna.txt"
file = open(moderna_placebo_path, 'r')
Lines = file.readlines()
k_inf_placebo_moderna = []

for line in Lines:
    k_inf_placebo_moderna.append(float(line.replace('\n', '')))

# vaccination arm
moderna_vacc_path = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Moderna_Vacc.txt"
file = open(moderna_vacc_path, 'r')
Lines = file.readlines()
k_inf_vacc_moderna = []

for line in Lines:
    k_inf_vacc_moderna.append(float(line.replace('\n', '')))

efficiacy_moderna = []

for placebo, vacc in zip(k_inf_placebo_moderna, k_inf_vacc_moderna):
    efficiacy_moderna.append(1 - (vacc/placebo))




mean_efficiacy_moderna = round(statistics.mean(efficiacy_moderna),4)
stdev_efficiacy_moderna = round(statistics.stdev(efficiacy_moderna),4)

fig, ax = plt.subplots()
plt.hist(efficiacy_moderna, bins=30)
plt.grid()
plt.title("Efficiacy Moderna")
plt.xlabel("Efficiacy")
plt.ylabel("Frequency")
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = "Uncertainity (Mean +- standard deviation):\n{} +- {} "
ax.text(0.05, 0.95, textstr.format(mean_efficiacy_moderna, stdev_efficiacy_moderna), transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
plt.show()



# Astra Zeneca _______________________________________________________________________________________________________________________

# placebo arm
astra_placebo_path = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Astra.txt"
file = open(astra_placebo_path, 'r')
Lines = file.readlines()
k_inf_placebo_astra = []

for line in Lines:
    k_inf_placebo_astra.append(float(line.replace('\n', '')))

# vaccination arm
astra_vacc_path = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Astra_Vacc.txt"
file = open(astra_vacc_path, 'r')
Lines = file.readlines()
k_inf_vacc_astra = []

for line in Lines:
    k_inf_vacc_astra.append(float(line.replace('\n', '')))

efficiacy_astra = []

for placebo, vacc in zip(k_inf_placebo_astra, k_inf_vacc_astra):
    efficiacy_astra.append(1 - (vacc/placebo))




mean_efficiacy_astra = round(statistics.mean(efficiacy_astra),4)
stdev_efficiacy_astra = round(statistics.stdev(efficiacy_astra),4)

fig, ax = plt.subplots()
plt.hist(efficiacy_astra, bins=30)
plt.grid()
plt.title("Efficiacy Astra Zeneca")
plt.xlabel("Efficiacy")
plt.ylabel("Frequency")
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = "Uncertainity (Mean +- standard deviation):\n{} +- {} "
ax.text(0.05, 0.95, textstr.format(mean_efficiacy_astra, stdev_efficiacy_astra), transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
plt.show()



vaccinations = ["Biontech", "Moderna", "Astra Zeneca"]
uncertainities = [stdev_efficiacy_biontech, stdev_efficiacy_moderna, stdev_efficiacy_astra]
plt.bar(vaccinations, uncertainities)
plt.grid()
plt.xlabel("Vaccination")
plt.ylabel("Standard deviation of efficiacy")
plt.show()