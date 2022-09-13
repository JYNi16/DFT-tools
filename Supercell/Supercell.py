#This script is coded for making supercell
#Author: njy 

import math 
import re 
import numpy as np 

class make_sc():
    
    def __init__(self, da, db, dc):
        
        self.da = da
        self.db = db 
        self.dc = dc 
        self.sup = da*db*dc
        self.sup_l = [da, db, dc]
        self.f1 = open('POSCAR','r')
        self.f2 = open('supercell.vasp', 'a')
    
    def ele_name(self, ele_or):       
        ele_cn = ''
        ele_s = [str(st) for st in ele_or.split()]
        for i in range(len(ele_s)):
            ele_cn += "    "
            ele_cn += str(ele_s[i])    
        return ele_cn
    
    def ele_lat(self, ele_or):       
        ele_cn = ''
        ele_s = [float(st) for st in ele_or.split()]
        for i in range(len(ele_s)):
            ele_cn += "    "
            ele_cn += str(np.round(ele_s[i] * self.sup_l[i], 7))   
        return ele_cn
    
    def ele_num(self, ele_or):
        ele_cn = ""        
        ele_s = [int(st) for st in ele_or.split()]
        for i in range(len(ele_s)):
            ele_cn += "    "
            ele_cn += str(ele_s[i] * self.sup) 
        nums = sum(ele_s) 
        return ele_cn, nums
    
    def ele_coord(self, ele_or, j, k, l, delta_a, delta_b, delta_c):
        ele_cn = ""
        sup_idx = [j, k, l]
        sup_delta = [delta_a, delta_b, delta_c]
        ele_s = [float(st) for st in ele_or.split()]
        for i in range(len(ele_s)):
            ele_cn += "    "
            ele_cn += str(np.round(((ele_s[i]/self.sup_l[i])+(sup_idx[i]*sup_delta[i]))%1,7))
            print
        
        return ele_cn 
    
      
    def fill(self):
        
        #f2 = open('supercell.vasp', 'a')
        
        self.f1.readline()
        self.f1.readline()
        print('supercell', file=self.f2)
        print('1.0', file=self.f2)
        a = self.f1.readline()
        lat_a = self.ele_lat(a)
        print(lat_a, file=self.f2)
        b = self.f1.readline()
        lat_b = self.ele_lat(b)
        print(lat_b, file=self.f2)
        c = self.f1.readline()
        lat_c = self.ele_lat(c)
        print(lat_c, file=self.f2)
        elename = self.f1.readline()
        ele_name = self.ele_name(elename)
        print(ele_name, file=self.f2)
        
        elenum = self.f1.readline()
        ele_num, total_elenum = self.ele_num(elenum)
        print("total_elenum is:", total_elenum)
        print(ele_num, file=self.f2)
        self.f1.readline()
        print('Direct', file=self.f2)
         
        for i in range(total_elenum):
            idx = self.f1.readline()
            print("idx is:", idx)
            #idx = [float(s) for s in idx.split()]
            delta_a = round(1/self.da, 6)
            delta_b = round(1/self.db, 6)
            delta_c = round(1/self.dc, 6)
            for j in range(self.da): 
                for k in range(self.db):
                    for l in range(self.dc):
                        coord = self.ele_coord(idx, j, k, l, delta_a, delta_b, delta_c)
                        print(coord, file=self.f2)
                        print("coord is:", coord)
    		
    def close(self):
        self.f1.close()
        self.f2.close()

if __name__=="__main__":
    sup = make_sc(1,2,2)
    sup.fill()
    sup.close()
    