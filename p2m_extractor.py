import scipy.io as sio
import os
import glob

data_fmt = {'cir': ['pn', 'phase', 'delay', 'power'], 'doa': ['pn', 'AoA', 'EoA', 'power'],
            'dod': ['pn', 'AoD', 'EoD', 'power']}


class P2MConnector:
    def __init__(self, target: str = ''):
        assert os.path.exists(os.path.abspath(target)) and os.path.isdir(os.path.abspath(target)), 'Requested target does not exists!'
        self.path = os.path.abspath(target)
        self.files = glob.glob(target+'*.p2m')
        self.rxs = dict()
        self.txs = dict()

    def extract_params(self, datas:dict = {'cir': 2, 'dod': 1, 'doa': 1}):
        data_dict = dict()
        for i in self.files:
            data_desc = os.path.basename(i).split(sep='.')[1]
            TX_desc = os.path.basename(i).split(sep='.')[2]
            RX_desc = os.path.basename(i).split(sep='.')[3]
            if data_desc in datas.keys():
                file = open(i)
                for j in range(datas[data_desc]):
                    file.readline()
                numpoints = int(file.readline())

                for j in range(numpoints):
                    infovec = file.readline().split(' ')
                    ckey = 'TX[{}]->RX[{}]'.format(TX_desc, j)

                    if ckey not in data_dict.keys():
                        data_dict[ckey] = dict()

                    for k in range(int(infovec[1])):
                        datavec = file.readline().split(' ')
                        lowdict = dict()
                        for l in datavec:
                            lowkey = data_fmt[data_desc][datavec.index(l)]
                            lowdict[lowkey] = float(l)

                            #print(l, end=' ')

                        if k not in data_dict[ckey]:
                            data_dict[ckey][k] = lowdict
                        else:
                            data_dict[ckey][k].update(lowdict)

        for i in data_dict.items():
            datadd = dict()
            for j in i[1].items():
                for k in j[1].keys():
                    if k not in datadd.keys():
                        datadd[k] = list()
                        datadd[k].append(j[1][k])
                    else:
                        datadd[k].append(j[1][k])
            sio.savemat(file_name=i[0] + '.mat', mdict=datadd)

        #for i in data_dict.keys():
        #    sio.savemat(file_name=i+'.mat', mdict=data_dict[i])


if __name__ == '__main__':
    pc = P2MConnector(target='/home/alexey/sss/')
    pc.extract_params()