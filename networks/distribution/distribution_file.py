from .idistribution import IDiameterDistribution
from pathlib import Path
import csv
import numpy as np

class CSVFileDiameterDistribution(IDiameterDistribution):
    def __init__(self, path: Path):
        self.xs = []
        self.ys = []
        with open(path, 'r') as csvfile:
            reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
            _header = reader.__next__()
            # print(header)
            for data in reader:
                self.xs.append(np.float64(data[0]))
                self.ys.append(np.float64(data[1]))

        self.ys = np.array(self.ys, dtype=np.float64) - self.ys[0]
        self.ys = self.ys / self.ys[-1] # normalized
        self.xs = np.array(self.xs, dtype=np.float64) * 1e-10 # A to m

    
    def generate_n_diameters(self,  num: int):
        seeds = np.random.rand(num)
        Dpores = np.interp(seeds, self.ys, self.xs)
        return Dpores