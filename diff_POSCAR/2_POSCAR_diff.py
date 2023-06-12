# Author Jinyang Ni
# This script is coded for generating the diff files between two POSCARs such as FE and PE
import numpy as np 
import math 


def read_files(filename):
    '''
    str:filename
    ion_x, ion_y, ion_z : the
    '''
    ion_x, ion_y, ion_z = [], [], [] 
    
    with open(filename,'r') as f1: 
    	f1.readline()
    	f1.readline()
    	a1 = f1.readline() 
    	a1 = [float(s) for s in a1.split()]
    	b1 = f1.readline()
    	b1 = [float(s) for s in b1.split()]
    	c1 = f1.readline()
    	c1 = [float(s) for s in c1.split()]
    	elename = f1.readline()
    	elename = elename.split()
    	elenum = f1.readline() 
    	elenum = [int(s) for s in elenum.split()]
    	num = elenum[0] + elenum[1] + elenum[2]
    	f1.readline()
    	lines = f1.readlines()
    	for line in lines:
    		value = [float(s) for s in line.split()]
    		ion_x.append(value[0])
    		ion_y.append(value[1])
    		ion_z.append(value[2])
    f1.close() 
    
    return ion_x, ion_y, ion_z, elenum, elename, a1, b1, c1

def diff_pos_FE_PE(N, ion_PE_x, ion_PE_y, ion_PE_z, ion_FE_x, ion_FE_y, ion_FE_z):
    
    delta_x = [] 
    delta_y = []
    delta_z = [] 
    
    print("ion_PE_x is:", len(ion_PE_x), len(ion_FE_x))
    
    for i in range(len(ion_PE_x)):
    	a = ion_FE_x[i]-ion_PE_x[i]
    	delta_x.append((a - round(a))/(N-1))
    	b = ion_FE_y[i]-ion_PE_y[i]
    	delta_y.append((b - round(b))/(N-1))
    	c = ion_FE_z[i]-ion_PE_z[i]
    	delta_z.append((c - round(c))/(N-1))
    
    return delta_x, delta_y, delta_z  

def write_pos(num, ion_PE_x, ion_PE_y, ion_PE_z, elenum, elename, a1, b1, c1, a2, b2, c2, delta_x, delta_y, delta_z):
    print("elenum is:", elenum, "elename is:", len(elename))
    
    for i in range(num): 
    	s = 'POSCAR' + '_'+str(i)
    	with open(s,'w') as f3: 
    		print(elename[0],'  ', elename[1],'   ',elename[2], file = f3)
    		print('1.0',file=f3)
    		print('   ',((a2[0]-a1[0])/num)*i + a1[0],'   ', ((a2[1]-a1[1])/num)*i + a1[1],'   ', ((a2[2]-a1[2])/num)*i + a1[2], file=f3)
    		print('   ',((b2[0]-b1[0])/num)*i + b1[0],'   ', ((b2[1]-b1[1])/num)*i + b1[1],'   ', ((b2[2]-b1[2])/num)*i + b1[2], file=f3)
    		print('   ',((c2[0]-c1[0])/num)*i + c1[0],'   ', ((c2[1]-c1[1])/num)*i + c1[1],'   ', ((c2[2]-c1[2])/num)*i + c1[2], file=f3)
    		print(' ',elename[0],'  ',elename[1],'   ', elename[2], file = f3)
    		print(' ',elenum[0],'   ',elenum[1],'   ', elenum[2], file = f3)
    		print('Direct',file = f3)
    		for j in range(0,1):
    			print('  ',(i*delta_x[j]+ion_PE_x[j])%1.0,'   ', (i*delta_y[j]+ion_PE_y[j])%1.0,'   ', (i*delta_z[j]+ion_PE_z[j])%1.0, file=f3)
    		for k in range(1,2):
    			print('  ',(i*delta_x[k]+ion_PE_x[k])%1.0,'   ', (i*delta_y[k]+ion_PE_y[k])%1.0,'   ', (i*delta_z[k]+ion_PE_z[k])%1.0, file=f3)
    		for l in range(2,5):
    			print('  ',(i*delta_x[l]+ion_PE_x[l])%1.0,'   ', (i*delta_y[l]+ion_PE_y[l])%1.0,'   ', (i*delta_z[l]+ion_PE_z[l])%1.0, file=f3)
    		#for l in range(elenum[2],elenum[3]):
    			#print((delta_x[l]+ion_PE_x[l])%1.0,'   ', (delta_y[l]+ion_PE_y[l])%1.0,'   ', (delta_z[l]+ion_PE_y[l])%1.0, file=f3)


def main(diff_num):
    
    ion_PE_x, ion_PE_y, ion_PE_z, elenum, elename, a1, b1, c1 = read_files("POSCAR-PE.vasp")
    ion_FE_x, ion_FE_y, ion_FE_z, elenum, elename, a2, b2, c2 = read_files("POSCAR-FE.vasp")
    delta_x, delta_y, delta_z = diff_pos_FE_PE(diff_num, ion_PE_x, ion_PE_y, ion_PE_z, ion_FE_x, ion_FE_y, ion_FE_z)
    write_pos(diff_num, ion_PE_x, ion_PE_y, ion_PE_z, elenum, elename, a1, b1, c1, a2, b2, c2, delta_x, delta_y, delta_z)


if __name__=="__main__":
    
    diff_num = 18
    main(diff_num)