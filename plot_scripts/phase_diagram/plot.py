# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:59:29 2023

@author: 26526
"""

import numpy as np
import matplotlib.pyplot as plt

J = [] 
A = []
C = []

with open("phase.txt", "r") as f:
    value = f.readlines()
    for i in range(len(value)):
        v = value[i].split()
        J.append(float(v[0]))
        A.append(float(v[1]))
        C.append(int(v[2]))
    
    f.close() 

plt.scatter(J, A, cmap=C)
    

    