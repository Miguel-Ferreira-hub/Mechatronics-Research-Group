import os
from scipy.integrate import simpson
import pandas as pd

# ---------------------------------------------------------------------------
# Data Paths
# ---------------------------------------------------------------------------

# Cell 12 - Cell 10 - Cell 1 - Old Cell (1)
directory = [r'C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Cell 12',
    r'C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Cell 10',
    r'C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Cell 1',
    r'C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Old Data'] 

file_name = ['New Cell 12 Full Discharge.txt',
    'Full Discharge.txt',
    'Cell 1 Full Discharge.txt',
    'FULL DISCHARGE FAST FILTERED.txt']

labels = ['Battery 12','Battery 10','Battery 01','Old Battery 01']

cell12_charge = [['Cell 12 Charge to 5%.txt','Cell 12 Charge to 10%.txt','Cell 12 Charge to 15%.txt','Cell 12 Charge to 20%.txt',
    'Cell 12 Charge to 25%.txt','Cell 12 Charge to 30%.txt','Cell 12 Charge to 35%.txt','Cell 12 Charge to 40%.txt','Cell 12 Charge to 45%.txt',
    'Cell 12 Charge to 50%.txt','Cell 12 Charge to 55%.txt','Cell 12 Charge to 60%.txt','Cell 12 Charge to 65%.txt','Cell 12 Charge to 70%.txt',
    'Cell 12 Charge to 75%.txt','Cell 12 Charge to 80%.txt','Cell 12 Charge to 85%.txt','Cell 12 Charge to 90%.txt','Cell 12 Charge to 95%.txt',
    'Cell 12 Charge to 100%.txt']]

cell10_charge = [['Cell 10 Charge to 5%.txt','Cell 10 Charge to 10%.txt','Cell 10 Charge to 15%.txt','Cell 10 Charge to 20%.txt',
    'Cell 10 Charge to 25%.txt','Cell 10 Charge to 30%.txt','Cell 10 Charge to 35%.txt','Cell 10 Charge to 40%.txt','Cell 10 Charge to 45%.txt',
    'Cell 10 Charge to 50%.txt','Cell 10 Charge to 55%.txt','Cell 10 Charge to 60%.txt','Cell 10 Charge to 65%.txt','Cell 10 Charge to 70%.txt',
    'Cell 10 Charge to 75%.txt','Cell 10 Charge to 80%.txt','Cell 10 Charge to 85%.txt','Cell 10 Charge to 90%.txt','Cell 10 Charge to 95%.txt',
    'Cell 10 Charge to 100%.txt']]

cell1_charge = [['Cell 1 Charge to 5%.txt','Cell 1 Charge to 10%.txt','Cell 1 Charge to 15%.txt','Cell 1 Charge to 20%.txt',
    'Cell 1 Charge to 25%.txt','Cell 1 Charge to 30%.txt','Cell 1 Charge to 35%.txt','Cell 1 Charge to 40%.txt','Cell 1 Charge to 45%.txt',
    'Cell 1 Charge to 50%.txt','Cell 1 Charge to 55%.txt','Cell 1 Charge to 60%.txt','Cell 1 Charge to 65%.txt','Cell 1 Charge to 70%.txt',
    'Cell 1 Charge to 75%.txt','Cell 1 Charge to 80%.txt','Cell 1 Charge to 85%.txt','Cell 1 Charge to 90%.txt','Cell 1 Charge to 95%.txt',
    'Cell 1 Charge to 100%.txt']]

old_charge = [[None]]

old_soc = [0.00,25.4,34.8,44.1,53.4,62.7,72.0,81.4,90.7,100.0]

charge_file_names = cell12_charge + cell10_charge + cell1_charge + old_charge

# ---------------------------------------------------------------------------
# Compute Capacity, SoC and SoH 
# ---------------------------------------------------------------------------

# Calculate capacity and SoH
def compute(directory,file_name):
    path = os.path.join(directory,file_name)

    data = pd.read_csv(path,sep="\t")

    time = data.iloc[:,0]
    current = data.iloc[:,1]

    capacity = abs(simpson(current,time) / 3600)
    SoH = (capacity / 2100) * 100

    return capacity, SoH

# Calculate SoCs
def compute_soc(capacity,directory,charge_file_names):
    states = []
    soc = 0
    states.append("0.00%")
    
    for file in charge_file_names:
        if file is not None:
            path = os.path.join(directory,file)
            data = pd.read_csv(path,sep="\t")
            time = data.iloc[:,0]
            current = data.iloc[:,1]
            soc += (simpson(current,time) / (capacity*3600)) * 100
            result = f"{soc:.2f}%"
            states.append(result)

    return states

# Print values
count = 0
for d, f, label in zip(directory,file_name,labels):
    capacity, SoH = compute(d,f)

    print(f"Capacity ({label}): {capacity}")
    print(f"State of Health ({label} - SoH): {SoH:.2f}%")

    states = compute_soc(capacity,d,charge_file_names[count])
    count += 1
    print(f"{label} SoC: {states}")

print(f"{label[-1]} SoC: {old_soc}")