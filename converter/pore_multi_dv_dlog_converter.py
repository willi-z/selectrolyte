import csv
import numpy as np
from pathlib import Path




folder = Path.cwd() / "data/poresizes/"
specimens = ["O1_50", "O2_40", "O4_40"]

volume_shrinkage = {
   "O1_50": 0.274, 
   "O2_40": 0.25, 
   "O4_40": 0.17, 
}

for specimen in specimens:
    datas_dV = [] # pore size
    with open(str(folder / "dv_dlog" / (specimen + '.csv')), 'r') as csvfile:
        reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
        header = reader.__next__()
        print(header[0], ", ", header[1])
        for data in reader:
            x,y = data[0], data[1]
            if x == '' or y == '':
                break
            dV_dw = np.float64(y) * np.float64(x) # dV/dlog(w) * w = dV/dw [cm^3/g/nm]
            datas_dV.append([np.float64(x)*10.0 / np.power(1.0 - volume_shrinkage[specimen], 1./3.), np.float64(dV_dw)*10]) # [nm to A]

    datas_dV.sort(key=lambda e: e[0])

    vol = 0.0
    datas_V = []
    for i in range(len(datas_dV)-1):
        p0 = datas_dV[i]
        p1 = datas_dV[i+1]
        width0, dV0 = p0[0], p0[1]
        width1, dV1 = p1[0], p1[1]
        dV = (dV0 + dV1) / 2.0
        width = (width0 + width1) / 2.0
        vol += dV * (width1 - width0)
        datas_V.append([width, vol])

    file = folder / "comulative_pore_volume" / (specimen + '.csv')
    with open(str(file), 'w') as csvfile:
        writer = csv.writer(csvfile,quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["pore width [A]", "comulative pore volume [cm^3/g]"])
        writer.writerows(datas_V)