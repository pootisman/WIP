import numpy as np
from pywt import dwt
from matplotlib import pyplot as plt
from random import gauss

class ClusterDetector:
    # Winsize = window size in nanoseconds
    # sampstep = step between individual samples in nanoseconds
    def __init__(self, sampstep: float = -1.0):
        self.sampstep = sampstep

    # Makes discrete CIR, return Moving Average Ratio of it
    def __MAR_calc__(self, M: int, CIRp: list = [], CIRd: list = []):
        res = []

        if self.sampstep == -1.0:
            self.sampstep = np.min(np.diff(CIRd))/2.0

        discreteCIRd = np.linspace(start=np.min(CIRd), num=((np.max(CIRd))/self.sampstep),
                                   stop=(np.max(CIRd)))
        print(discreteCIRd)
        discreteCIRp = np.interp(discreteCIRd, CIRd, CIRp)
        print(discreteCIRp)

        for i in range(M, (len(discreteCIRp) - 1 - M)):
            res.append(10.0 * np.log(np.sum((discreteCIRp[i:(i+M-1)])**2.0) / np.sum((discreteCIRp[i-M:i-1])**2.0)))

        return res

    def Detect(self, CIRp: list = [0]*40+[np.abs(gauss(30/(i+1), 1)) for i in range(1000)], CIRd: list = [i for i in range(1040)], M: int = 250, eps: float = 0.1):
        (MARa, MARc) = dwt(self.__MAR_calc__(M=M, CIRp=CIRp, CIRd=CIRd), 'db2', mode='zero')

        print(MARc)

        MARcf = [i if i > 0 else 0 for i in MARc]

        print(MARcf)

        MARcf = MARcf / max(MARcf)

        print(MARcf)

        MARcfd = [i for i in MARcf if i > eps]

        print(MARcfd)

        threshold = (1 - np.std(MARcfd))/np.mean(MARcfd)

        print(threshold)

        print(np.max(MARcf))

        discreteCIRd = np.linspace(start=0.0, num=((np.max(CIRd))/self.sampstep),
                                   stop=(np.max(CIRd)))
        discreteCIRp = np.interp(discreteCIRd, CIRd, CIRp)

        plt.plot(MARcf)
        plt.plot(discreteCIRp/max(discreteCIRp))
        plt.show()











