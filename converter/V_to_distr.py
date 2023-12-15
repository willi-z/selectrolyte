import csv
import numpy as np


datas= []

diameters = []
volumes = []
vol_scale = 1.0 # 1.27

with open('./data/pore_V_data.csv', 'r') as csvfile:
    reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
    header = reader.__next__()
    print(header)
    for data in reader:
        diameter = np.float64(data[0]) # [A]
        vol = np.float64(data[1]) # [cm^3/g] 
        vol *= 100000000**3 # -> [A**3/g]

        # account for volume shrinking
        diameter *= 1.0/np.power(vol_scale, 1.0/3)
        vol *= 1.0/vol_scale

        diameters.append(diameter)
        volumes.append(vol)

for i in range(len(diameters)):
    diameter = diameters[i]
    vol = volumes[i]
    if i > 0:
        vol -= volumes[i-1]
    vol_sphere = np.pi * diameter**2 / 3.0
    number = vol / vol_sphere
    if i > 0:
        number += datas[i-1][1]
    datas.append([diameter, number])

with open('./data/pore_distr_data.csv', 'w') as csvfile:
    writer = csv.writer(csvfile,quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["pore width [A]", "comulative number of pores []"])
    writer.writerows(datas)