# This code is used to make a supercell 
import numpy as np 
import math 
import re 

s = input("please input the da db dc:")
s = [float(s) for s in s.split()]
da = int(s[0])
db = int(s[1])
dc = int(s[2])
sup = da*db*dc

with open('supercell.vasp','a') as f2:
	with open('POSCAR','r') as f1: 
		f1.readline()
		f1.readline()
		print('supercell',file=f2)
		print('1.0',file=f2)
		a = f1.readline()
		a = [float(s) for s in a.split()]
		print('  ',round(a[0]*da,7),'  ',round(a[1]*da,7),'  ',round(a[2]*da,7),file=f2) 
		b = f1.readline()
		b = [float(s) for s in b.split()]
		print('  ',round(b[0]*db,7),'  ',round(b[1]*db,7),'  ',round(b[2]*db,7),file=f2)
		c = f1.readline()
		c = [float(s) for s in c.split()]
		print('  ',round(c[0]*dc,7),'  ',round(c[1]*dc,7),'  ',round(c[2]*dc,7),file=f2) 
		elename = f1.readline()
		elename = [str(s) for s in elename.split()]
		elenum = f1.readline()
		elenum = [int(s) for s in elenum.split()]
		f1.readline()
		print('   ',elename[0],'   ',elename[1],file=f2)
		print('   ',elenum[0]*sup,'   ',elenum[1]*sup,'  ',file=f2)
		print('Direct',file=f2)
		total_elenum = elenum[0] + elenum[1]
		for i in range(total_elenum):
			i = f1.readline()
			i = [float(s) for s in i.split()]
			delta_a = round(1/da, 6)
			delta_b = round(1/db, 6)
			delta_c = round(1/dc, 6)
			for j in range(da): 
				for k in range(db):
					for l in range(dc):
						print('    ',round(((i[0]/da)+(j*delta_a))%1,7),'   ', round(((i[1]/db)+(k*delta_b))%1,7),'   ', round(((i[2]/dc)+(l*delta_c))%1,7),file=f2)




