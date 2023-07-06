# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:57:29 2023

@author: Ni Jinyang 
"""
import numpy as np

def not_empty(s):
    return s and s.strip()

class diff_POSCAR():
    
    def __init__(self, diff_nums):
        super(diff_POSCAR, self)
        self.num = diff_nums
        self.ion_PE_x, self.ion_PE_y, self.ion_PE_z = [], [], []
        self.ion_FE_x, self.ion_FE_y, self.ion_FE_z = [], [], []
        self.PE_lat_x, self.PE_lat_y, self.PE_lat_z = 0, 0, 0 
        self.FE_lat_x, self.FE_lat_y, self.FE_lat_z = 0, 0, 0
        self.delta_x, self.delta_y, self.delta_z = [], [], []
        self.ele_num = 0
        self.elename, self.elenum = []
    
    def write_PE_pos(self, filename):
       
        with open(filename,'r') as f1:
            f1.readline()
            f1.readline()
            a = f1.readline() 
            b = f1.readline() 
            c = f1.readline() 
            
            if "PE" in filename:
            	self.PE_lat_x = [float(s) for s in a.split()]
            	self.PE_lat_y = [float(s) for s in b.split()]
            	self.PE_lat_z = [float(s) for s in c.split()]
            else:
            	self.FE_lat_x = [float(s) for s in a.split()]
            	self.FE_lat_y = [float(s) for s in b.split()]
            	self.FE_lat_z = [float(s) for s in c.split()]
            
            #elename = f1.readline()
            self.elename = f1.readline().split()
            elenum = f1.readline() 
            elenum = [int(s) for s in elenum.split()]
            self.ele_num = elenum[0] + elenum[1] + elenum[2]
            f1.readline()
            lines = f1.readlines()
            lines = filter(not_empty, lines)
            
            if "PE" in filename:
                for line in lines:
                    value = [float(s) for s in line.split()]
                    self.ion_PE_x.append(value[0])
                    self.ion_PE_y.append(value[1])
                    self.ion_PE_z.append(value[2])
            else:
                for line in lines:
                    value = [float(s) for s in line.split()]
                    self.ion_FE_x.append(value[0])
                    self.ion_FE_y.append(value[1])
                    self.ion_FE_z.append(value[2])            
        
        f1.close() 
        
    def diff_PE_FE(self):
        
        for i in range(len(self.ion_PE_x)):
        	a = self.ion_FE_x[i] - self.ion_PE_x[i]
        	self.delta_x.append((a - round(a))/(self.num-1))
        	b = self.ion_FE_y[i] - self.ion_PE_y[i]
        	self.delta_y.append((b - round(b))/(self.num-1))
        	c = self.ion_FE_z[i] - self.ion_PE_z[i]
        	self.delta_z.append((c - round(c))/(self.N-1))
    
    def write_pos(self):
        
        for i in range(self.num): 
        	s = 'POSCAR' + '_'+str(i)
        	with open(s,'w') as f3: 
        		print(self.elename[0],'  ', self.elename[1],'   ', self.elename[2], file = f3)
        		print('1.0',file=f3)
        		print('   ',((self.FE_lat_x[0]-self.PE_lat_x[0])/self.num)*i + self.PE_lat_x[0],'   ', ((self.FE_lat_x[1]-self.PE_lat_x[1])/self.num)*i + self.PE_lat_x[1],'   ', ((self.FE_lat_x[2]-self.PE_lat_x[2])/self.num)*i + self.PE_lat_x[2], file=f3)
        		print('   ',((self.FE_lat_y[0]-self.PE_lat_y[0])/self.num)*i + self.PE_lat_y[0],'   ', ((self.FE_lat_y[1]-self.PE_lat_y[1])/self.num)*i + self.PE_lat_y[1],'   ', ((self.FE_lat_y[2]-self.PE_lat_y[2])/self.num)*i + self.PE_lat_y[2], file=f3)
        		print('   ',((self.FE_lat_z[0]-self.PE_lat_z[0])/self.num)*i + self.PE_lat_z[0],'   ', ((self.PE_lat_z[1]-self.PE_lat_z[1])/self.num)*i + self.PE_lat_z[1],'   ', ((self.PE_lat_z[2]-self.PE_lat_z[2])/self.num)*i + self.PE_lat_z[2], file=f3)
        		print('  ',self.elename[0],'  ',self.elename[1],'   ', self.elename[2], file = f3)
        		print('  ',self.elenum[0],'   ',self.elenum[1],'   ', self.elenum[2], file = f3)
        		print('Direct',file = f3)
        		for j in range(0,1):
        			print('  ',(i*self.delta_x[j]+self.ion_PE_x[j])%1.0,'   ', (i*self.delta_y[j]+self.ion_PE_y[j])%1.0,'   ', (i*self.delta_z[j]+self.ion_PE_z[j])%1.0, file=f3)
        		for k in range(1,2):
        			print('  ',(i*self.delta_x[k]+self.ion_PE_x[k])%1.0,'   ', (i*self.delta_y[k]+self.ion_PE_y[k])%1.0,'   ', (i*self.delta_z[k]+self.ion_PE_z[k])%1.0, file=f3)
        		for l in range(2,5):
        			print('  ',(i*self.delta_x[l]+self.ion_PE_x[l])%1.0,'   ', (i*self.delta_y[l]+self.ion_PE_y[l])%1.0,'   ', (i*self.delta_z[l]+self.ion_PE_z[l])%1.0, file=f3)
        		#for l in range(elenum[2],elenum[3]):
        			#print((delta_x[l]+ion_PE_x[l])%1.0,'   ', (delta_y[l]+ion_PE_y[l])%1.0,'   ', (delta_z[l]+ion_PE_y[l])%1.0, file=f3)