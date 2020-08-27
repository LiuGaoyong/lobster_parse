#!/usr/bin/env python3
import os, sys, re
import numpy as np
import pandas as pd 
from scipy.integrate import simps


class LobsterParse():
    '''
    用于综合解析LOBSTER结果文件的类
    '''
    def __init__(self, lobster_path:str):
        if os.path.exists(os.path.join(lobster_path, 'lobsterout')):
            self.file_list = os.listdir(lobster_path)
            with open(os.path.join(lobster_path, 'lobsterout')) as f:
                data = f.read().split('\n')
            spin_line = [line for line in data if 'spin' in line]
            if len(spin_line) > 0:
                self.is_spin_polarized = True
            else:
                self.is_spin_polarized = False
        elif os.path.exists(os.path.join(lobster_path, 'lobsterin')):
            print(lobster_path, "is a LOBSTER project.\n" + "    but output files don't exist !!!")
            exit(101)
        else:
            print(lobster_path, "isn't a LOBSTER project !!!")
            exit(102)
        self.root_path  =   lobster_path


    def _cohp_as_str(self, is_coop=False):
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
        data = self._cohp_as_str(is_coop)
        data = data.splitlines()
        second_line = data[1].split()
        num_bonds   = int(second_line[0])
        bonds       = data[2:num_bonds+2]
        bonds_data  = data[num_bonds+2:]
        columns     = ['Energy(eV)']
        columns.append( 'COHP('+bonds[0]+')') 
        columns.append('ICOHP('+bonds[0]+')')
        for i in range(1, num_bonds):
            bond_name   = bonds[i].split('(')[0].split(':')[-1]
            bond_length = bonds[i].split('(')[-1].split(')')[0]
            columns.append( 'COHP('+bond_name+')_'+bond_length+'') 
            columns.append('ICOHP('+bond_name+')_'+bond_length+'')
        result = []
        for i in bonds_data:
            result.append(i.split())
        result = pd.DataFrame(result, dtype=float)
        #对自旋极化和非极化体系的处理
        if self.is_spin_polarized:
            df1 = result.loc[:,0]
            df2 = pd.DataFrame(result.values[:, 1:num_bonds*2+1])
            df3 = pd.DataFrame(result.values[:, num_bonds*2+1:] )
            result_up = pd.concat([df1,df2], axis=1, ignore_index=True)
            result_up.columns = columns
            result_down = pd.concat([df1,df3], axis=1, ignore_index=True)
            result_down.columns = columns
            return (result_up + result_down)/2
        else:
            result.columns = columns
            return result

    def cohp_get_bond_from_str(self, bond_str:str, is_coop=False):
        df = self.cohp_as_DataFrame(is_coop)
        columns = [ df.columns[0] ]
        for name in df.columns[1:]:
            if bond_str in name and 'ICOHP' not in name:
                columns.append(name)
        return df.get(columns)

    def cohp_get_EMICOHP(self, bond_pd):
        #获取Fermi能级以下的数据
        data =  bond_pd.loc[bond_pd[bond_pd.columns[0]] <= 0.0]        
        #计算EMICOHP
        x = data[data.columns[0]]
        delta_x = round(np.mean([x[i+1]-x[i] for i in range(len(x)-1)]), 5)
        y = data[data.columns[1:]]
        ICOHP = np.sum(y * delta_x)
        yE = y.mul(x, axis=0)
        EMICOHP = np.sum(yE * delta_x) / ICOHP
        #print('  parse ICOHP:   \n', ICOHP)
        #print('  parse EMICOHP: \n', EMICOHP)
        return EMICOHP

    def cohp_get_EMICOHP_all_as_df(self):
        data = self.cohp_get_EMICOHP(self.cohp_as_DataFrame())
        result = [ [data.index[0].split('(')[-1].replace(')', ''), 
                    0, 
                    data.values[0], 
                    data.values[1] ] ]
        for i in range(2,len(data),2):
            bond_id = data.index[i].split('(')[-1].replace(')', '').split('_')
            bond_name, bond_length = bond_id[0], float(bond_id[1])
            cohp    = data.values[i]
            emicohp = data.values[i+1]
            result.append([bond_name, bond_length, cohp, emicohp])
        result = pd.DataFrame(result)        
        return result

    def test(self):
        #test1 = self.cohp_get_bond(bond=['Co1', 'Cu13'])
        #test1 = self.cohp_get_bond_from_str('C21->O22')
        #print(test1)
        #print(test1.info())
        #print(self.is_spin_polarized)
        #test2 = self.cohp_as_DataFrame()
        #test3 = self.cohp_get_EMICOHP(test1)
        #print(test2)
        #print(test3, type(test3))
        self.cohp_get_EMICOHP_all_as_df().to_csv('EMICOHPLIST.parse', sep='\t')
        print(self.cohp_get_EMICOHP_all_as_df())



if __name__ == "__main__":
    aaa = LobsterParse('../example/spin-2')
    aaa.test()