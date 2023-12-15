import csv
import numpy as np


datas_dV = [] # pore size
entry = 0

with open('./data/pore_dv_dlog_data.csv', 'r') as csvfile:
    reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
    header = reader.__next__()
    print(header[2*entry])
    header = reader.__next__()
    print(header[0+2*entry], ", ", header[1+2*entry])
    for data in reader:
        x,y = data[0+2*entry], data[1+2*entry]
        if x == '' or y == '':
            break
        datas_dV.append([np.float64(x)*10.0, np.float64(y)]) # [nm to A]

datas_dV.sort(key=lambda e: e[0])

vol = 0.0
datas_V = []
for entry in datas_dV:
    width, dV = entry[0], entry[1]
    vol += dV
    datas_V.append([width, vol])


with open('./data/pore_V_data.csv', 'w') as csvfile:
    writer = csv.writer(csvfile,quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["pore width [A]", "comulative pore volume [cm^3/g]"])
    writer.writerows(datas_V)