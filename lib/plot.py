#!/usr/bin/env python3
from base import LobsterParse
import matplotlib.pyplot as plt


class LobsterCohpcarPlot(LobsterParse):

    def plot(self, cohp_bond_str:str, energy_in_y=False):
        data = self.cohp_get_bond_from_str(cohp_bond_str)
        x = data.values[:, 0]
        x_label = 'E-Ef(eV)'
        y = data.values[:, 1] * (-1)
        y_label = '-COHP'
        if energy_in_y:
            x,y              = y,x
            x_label, y_label = y_label, x_label
            plt.figure(figsize=(5,7))
        title = '-COHP('+data.columns[1].split('(')[-1].split(')')[0]+')'
        
        plt.plot(x, y)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        return plt

    def test(self):
        #test1 = self.cohp_get_bond_from_str('C21->O22')
        #print(test1)
        #test3 = self.cohp_get_EMICOHP(test1)
        #print(test3, type(test3)) 
        aaa = self.plot('Cu17->C21', energy_in_y=True)
        aaa.savefig('test.png')


if __name__ == "__main__":
    aaa = LobsterCohpcarPlot('../example/spin-2')
    aaa.test()   