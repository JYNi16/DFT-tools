import numpy as np 
import math 

Sa = '  '
Sb = '  '
Sc = '  '

with open('231.vasp','r') as f: 
	f.readline()
	f.readline()
	a = f.readline()
	b = f.readline()
	c = f.readline()
	a = [float(s) for s in a.split()]
	b = [float(s) for s in b.split()]
	c = [float(s) for s in c.split()]
	elename = f.readline().split()
	elenum = f.readline().split()
	ion_pos = f.read()
	for i in range(0,10,2):
		Sa = '  ' + str(round(a[0]*(1 - i/100),6)) + '  ' + str(round(a[1]*(1 - i/100),6)) + '  ' + str(round(a[2]*(1 - i/100),6))
		Sb = '  ' + str(round(b[0]*(1 - i/100),6)) + '  ' + str(round(b[1]*(1 - i/100),6)) + '  ' + str(round(b[2]*(1 - i/100),6))
		Sc = '  ' + str(round(c[0]*(1 + 0/100),6)) + '  ' + str(round(c[1]*(1 + 0/100),6)) + '  ' + str(round(c[2]*(1 + 0/100),6))
		namefile = 'POSCAR_'+str(i)
		with open (namefile, 'w') as f2:
			print('POSCAR', file = f2)
			print(' 1.0', file = f2)
			print(Sa, file = f2)
			print(Sb, file = f2)
			print(Sc, file = f2)
			print(" ".join(elename),  file = f2)
			print(" ".join(elenum), file = f2)
			print('Direct', file = f2)
			print(ion_pos, file = f2)