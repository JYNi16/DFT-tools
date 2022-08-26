
import re
import numpy as np

with open('band.yaml','r') as f:
	while True:
		l=f.readline()
		s=re.match('(\w+):\s+(\d+)',l)

		if s is not None:
			if s.group(1)=='nqpoint':
				nqpoint=int(s.group(2))
			if s.group(1)=='npath':
				npath=int(s.group(2))
			if s.group(1)=='natom':
				natom=int(s.group(2))
				break


	nband=natom*3
	frequency=np.zeros((nqpoint,nband),dtype=float)
	distance=np.zeros(nqpoint,dtype=float)

	i=0
	for line in f:
		s1=re.match('\s+(\w+)\:\s+(\-?\d+\.\d+)',line)
		if s1 is not None:
			if s1.group(1)=='distance':
				distance[i]=float(s1.group(2))
				j=0
				while j<nband:
					s1=re.match('\s+(\w+)\:\s+(\-?\d+\.\d+)',f.readline())
					if s1 is not None:
						if s1.group(1)=='frequency':
							frequency[i][j]=float(s1.group(2))
							j=j+1
				i=i+1

#print(frequency[3])
with open('1.dat','w') as f1:
	i1=0
	while i1<nband:
		j1=0
		while j1<nqpoint:
			f1.write('%.7f    %.10f\n'%(distance[j1],frequency[j1][i1]))
			j1=j1+1
		f1.write('\n')
		i1=i1+1