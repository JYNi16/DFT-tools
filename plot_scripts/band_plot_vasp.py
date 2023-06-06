# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:58:14 2023
@author: Ni Jinyang 

1. This script is coded for plot the band structure
2. The data is preprocessed by the Vaspkit 

"""
import numpy as np
import matplotlib.pyplot as plt 

def write_band_file(file_path):
    with open(file_path, "r") as f:
        vs = f.readlines()
    
    num_bands = len(vs[1].split()) - 1
    k_v = []
    E_v = np.zeros((len(vs)-1, num_bands))
    for i in range(1, len(vs)):
        value = [float(s) for s in vs[i].split()]
        k_v.append(value[0])
        E_v_test = value[1:]
        E_v[i-1,:] = value[1:]    
        
    #print("kpoints is:", k_v)
    #print("E_v shape is:", E_v.shape)
    f.close() 
    
    return k_v, E_v 

def plot_band():
    k_v_dn, E_v_dn = write_band_file("REFORMATTED_BAND_DW.dat")
    k_v_up, E_v_up = write_band_file("REFORMATTED_BAND_UP.dat")
    num_ks, num_bands = E_v_dn.shape[0], E_v_dn.shape[1] 
    print("num_ks is:", num_ks, "num_bands is:", num_bands)
    
    for i in range(num_bands):
        print(len(k_v_dn), len(E_v_dn[:,i]))
        plt.plot(k_v_dn, E_v_dn[:,i], "black")
        plt.plot(k_v_up, E_v_up[:,i], "red")
    
    k_sym_points = [0.0, 1.033, 1.629, 2.821]
    k_sym_label = ["G", "M", "K", "G"]
    plt.ylim(-2, 6)
    plt.xlim(0, 2.82)
    plt.ylabel("Energy(eV)")
    plt.xlabel("Kpoints")
    plt.xticks(k_sym_points, k_sym_label)
    plt.show()
    
if __name__=="__main__":
    plot_band()