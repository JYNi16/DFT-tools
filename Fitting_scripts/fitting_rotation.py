# -*- coding: utf-8 -*-
"""
@author: Jinyang Ni
# This script is coded for fitting the magnetic exchange parameters (J and K)
# In the input file E.dat, one site's spin is fixed while the other change from 0 to Pi.
"""

import numpy as np 
import matplotlib.pyplot as plt
import math  
from scipy.optimize import curve_fit

p = math.pi

def Forth(x,a,b,c):
	return a*(np.cos(x)) + b*(np.cos(x)*np.cos(x)) + c

def Second(x,a,b):
    return a*(np.cos(x)) + b

def read_data():
    x = [] 
    y = []
    with open ("E.dat",'r') as f:
        data = f.readlines() 
        for i in range(len(data)): 
            value = [float(s) for s in data[i].split()]
            x.append(value[0]*p/18)
            y.append(value[1]*1000)
    f.close()
    
    const = y[0]
    for i in range(len(y)):
        y[i] -= const       
    return x, y


x, y = read_data()
popt, pcov = curve_fit(Second,x,y)
J = popt[0]
B = popt[1]
yvals = Second(x,J,B)
print('Now the Fitting process only include J :')
print('J is : {} meV'.format(round(J, 6)))
print('constant: {}'.format(round(B, 6)))

popt, pcov = curve_fit(Forth,x,y)
J_2 = popt[0]
K_2 = popt[1]
B_2 = popt[2]
y_K = Forth(x,J_2,K_2,B_2)
print("==========================================")
print('Now the Fitting process include J and K :')
#print('popt:' ,popt)
print('J is : {} meV'.format(round(J_2, 6)))
print('K is : {} meV'.format(round(K_2, 6)))
print('constant: {}'.format(round(B_2, 6)))

font = {'family': "Times New Roman", "weight":"normal", "size":20,}
plt.figure(1, figsize=(10,8))
plot1 = plt.scatter(x,y, s = 100, c="black", label='Tight Binding results')
plot2 = plt.plot(x,yvals,'r',label=r'Fitting $A+JCos(\theta)$', linewidth=4)
plot2 = plt.plot(x,y_K,'green',label=r'Fitting $A+JCos(\theta) + KCos^{2}(\theta)$', linewidth=4)
plt.xlabel(r'Rotation angle($\theta$)', font)
plt.ylabel('E($meV$)', font)
plt.xticks(fontproperties = "Times New Roman", fontsize=20)
plt.yticks(fontproperties = "Times New Roman", fontsize=20)
plt.legend(loc=4, prop = {'family': "Times New Roman", "weight":"normal", "size":18,})
plt.savefig("JK.png")
#plt.ylim(-0.5, 70)
plt.show()