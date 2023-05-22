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
            y.append(value[1])
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
print('J is : {} meV'.format(round(J*1000, 6)))
print('constant: {}'.format(round(B, 6)))

popt, pcov = curve_fit(Forth,x,y)
J_2 = popt[0]
K_2 = popt[1]
B_2 = popt[2]
y_K = Forth(x,J_2,K_2,B_2)
print("==========================================")
print('Now the Fitting process include J and K :')
#print('popt:' ,popt)
print('J is : {} meV'.format(round(J_2*1000, 6)))
print('K is : {} meV'.format(round(K_2*1000, 6)))
print('constant: {}'.format(round(B_2, 6)))

plot1 = plt.plot(x,y,'s',label='original values')
plot2 = plt.plot(x,yvals,'r',label='only J values')
plot2 = plt.plot(x,y_K,'b',label='add K values')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc=4)
plt.title('curve_fit')
plt.savefig("JK.png")
plt.show()