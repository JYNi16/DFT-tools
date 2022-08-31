# This code is used to calculate the differenece between two POSCARs 
import numpy as np 
import math 

s = input("the number of POSCARs between FE and PE structure")
N = int(s)

ion_PE_x = []
ion_PE_y = []
ion_PE_z = []
ion_FE_x = []
ion_FE_y = []
ion_FE_z = []
delta_x = []
delta_y = []
delta_z = []
 

with open('POSCAR-PE','r') as f1: 
	f1.readline()
	f1.readline()
	a1 = f1.readline() 
	a1 = [float(s) for s in a1.split()]
	b1 = f1.readline()
	b1 = [float(s) for s in b1.split()]
	c1 = f1.readline()
	c1 = [float(s) for s in c1.split()]
	elename = f1.readline()
	#elename = [char(s) for s in elename.split()]
	elenum = f1.readline() 
	elenum = [int(s) for s in elenum.split()]
	num = elenum[0] + elenum[1] + elenum[2]
	f1.readline()
	lines = f1.readlines()
	for line in lines:
		value = [float(s) for s in line.split()]
		ion_PE_x.append(value[0])
		ion_PE_y.append(value[1])
		ion_PE_z.append(value[2])

with open('POSCAR-FE','r') as f2: 
	f2.readline()
	f2.readline()
	a2 = f2.readline() 
	a2 = [float(s) for s in a2.split()]
	b2 = f2.readline()
	b2 = [float(s) for s in b2.split()]
	c2 = f2.readline()
	c2 = [float(s) for s in c2.split()]
	elename = f2.readline()
	elename = [str(s) for s in elename.split()]
	elenum = f2.readline() 
	elenum = [int(s) for s in elenum.split()]
	f2.readline()
	lines = f2.readlines()
	for line in lines:
		value = [float(s) for s in line.split()]
		ion_FE_x.append(value[0])
		ion_FE_y.append(value[1])
		ion_FE_z.append(value[2])
        

for i in range(num):
	a = ion_FE_x[i]-ion_PE_x[i]
	delta_x.append((a - round(a))/(N-1))
	b = ion_FE_y[i]-ion_PE_y[i]
	delta_y.append((b - round(b))/(N-1))
	c = ion_FE_z[i]-ion_PE_z[i]
	delta_z.append((c - round(c))/(N-1))
	
	
for i in range(N): 
	s = 'POSCAR' + '_'+str(i)
	with open(s,'w') as f3: 
		print(elename[0],'  ', elename[1],'   ',elename[2], file = f3)
		print('1.0',file=f3)
		print('   ',((a2[0]-a1[0])/N)*i + a1[0],'   ', ((a2[1]-a1[1])/N)*i + a1[1],'   ', ((a2[2]-a1[2])/N)*i + a1[2], file=f3)
		print('   ',((b2[0]-b1[0])/N)*i + b1[0],'   ', ((b2[1]-b1[1])/N)*i + b1[1],'   ', ((b2[2]-b1[2])/N)*i + b1[2], file=f3)
		print('   ',((c2[0]-c1[0])/N)*i + c1[0],'   ', ((c2[1]-c1[1])/N)*i + c1[1],'   ', ((c2[2]-c1[2])/N)*i + c1[2], file=f3)
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
        


		
