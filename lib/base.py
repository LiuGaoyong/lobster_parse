#!/usr/bin/env python3
import os, sys, re
import pandas as pd 
from scipy.integrate import simps


class LobsterParse():
    '''
    用于综合解析LOBSTER结果文件的类
    '''
    def __init__(self, lobster_path:str):
        if os.path.exists(os.path.join(lobster_path, 'lobsterout')):
            self.file_list = os.listdir(lobster_path)
        elif os.path.exists(os.path.join(lobster_path, 'lobsterin')):
            print(lobster_path, "is a LOBSTER project.\n" + "    but output files don't exist !!!")
            exit(101)
        else:
            print(lobster_path, "isn't a LOBSTER project !!!")
            exit(102)
        self.root_path  =   lobster_path

    def cohp_as_str(self, is_coop=False):
        if is_coop:
            target = 'COOPCAR.lobster'
        else:
            target = 'COHPCAR.lobster'
        if target in self.file_list:
            with open(os.path.join(self.root_path, target)) as f:
                data = f.read()
            return data
        else:
            print(os.path.join(self.root_path, target), "doesn't exist !!!")
            exit(103)

    def cohp_as_DataFrame(self, is_coop=False):
        data = self.cohp_as_str(is_coop)
        data = data.splitlines()
        second_line = data[1].split()
        num_bonds   = int(second_line[0])
        bonds       = data[2:num_bonds+2]
        bonds_data  = data[num_bonds+2:]
        columns     = ['Energy(eV)']
        columns.append( 'COHP('+bonds[0]+')') 
        columns.append('ICOHP('+bonds[0]+')')
        for i in range(1, num_bonds):
            bond_name = bonds[i].split('(')[0].split(':')[-1]
            columns.append( 'COHP('+bond_name+')') 
            columns.append('ICOHP('+bond_name+')')
        result = []
        for i in bonds_data:
            result.append(i.split())
        result = pd.DataFrame(result, columns=columns, dtype=float)
        return result

    def cohp_get_bond(self, bond='', is_coop=False):
        df = self.cohp_as_DataFrame(is_coop)
        columns = [ df.columns[0] ]
        if len(bond) == 2:
            print('sdfasdfas')
            for name in df.columns[1:]:
                if re.match('.*'+bond[0] + '[0-9]*->' + bond[1] + '[0-9]*', name):
                    print(name)
                    columns.append(name)
        elif len(bond) == 1:
            for name in df.columns[1:]:
                if re.match(bond[0] + '[0-9]*', name):
                    columns.append(name)
        else:
            columns = df.columns[:3]
        return df.get(columns)


    def cohp_get_EMICOHP(self, bond_pd):
        index_x, index_y = bond_pd.columns[:2]
        #获取Fermi能级以下的数据
        data =  bond_pd.loc[bond_pd[index_x] < 0.0]         
        #计算EMICOHP
        x = data.loc[:, index_x]
        y = data.loc[:, index_y]
        yE = x * y
        return round(simps(yE,x) / simps(y,x),  5)


    def test(self):
        test1 = self.cohp_get_bond(bond=['Co1', 'Cu13'])
        print(test1)
        #test2 = self.cohp_as_DataFrame()
        test3 = self.cohp_get_EMICOHP(test1)
        print(test3)



if __name__ == "__main__":
    aaa = LobsterParse('../example')
    aaa.test()