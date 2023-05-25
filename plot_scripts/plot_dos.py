# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:46:12 2023

@author: 26526
"""
import numpy as np
import matplotlib.pyplot as plt 

def read_data(filename):
    print("filename is:", filename)
    E, spin_up, spin_dn = [], [], [] 
    with open(filename+".dat", "r") as f:
        values = f.readlines()
        for i in range(1, len(values)):
            v = [float(s) for s in values[i].split()]
            E.append(v[0])
            spin_up.append(v[1])
            spin_dn.append(v[2])
    f.close() 
    return E, spin_up, spin_dn

def input_s(files):
    Energy = []
    updos = []
    dndos = []
    for f in files:
        E, sup, sdn = read_data(f)
        Energy.append(E)
        updos.append(sup)
        dndos.append(sdn)
    
    return Energy, updos, dndos

def plot(file_path): 
    E, sup, sdn = input_s(file_path) 
    for i in range(len(sup)):
        plt.plot(E[i],sup[i]) 
        plt.plot(E[i],sdn[i], label = file_path[i])
        plt.xlim(-10, 10)
        plt.legend(loc=4)
    plt.show()

if __name__=="__main__":
    file_path = ["dx2"]
    plot(file_path)