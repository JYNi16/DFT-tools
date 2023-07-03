# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:58:14 2023
@author: Ni Jinyang 

1. This script is coded for plot the band structure
2. Input datas include two files -- xx_UP.dat and xx__DN.dat, which are preprocessed by the Vaspkit 

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

def write_wannier_band():
    num_bands = 1
    K_v, E_v = [], []
    with open("wannier90_band.dat", "r") as f:
        wv = f.readlines() 
    for i in range(len(wv)):
        if wv[i].strip() == "":
            num_bands += 1
            pass
        else:
            value = [float(s) for s in wv[i].split()]
            K_v.append(value[0])
            E_v.append(value[1])
    
    k_points = int(len(K_v)/num_bands)
    E_v_w = np.zeros((k_points, num_bands))
    
    for j in range(num_bands):
        E_v_w[:,j] = np.array(E_v[j*k_points:(j+1)*k_points])
    
    #print("num_bands is:", E_v_w.shape, "K_points is:", k_points)
    
    return np.array(K_v[:274]), E_v_w, num_bands
            

def plot_band():
    k_v, E_v = write_band_file("REFORMATTED_BAND.dat")
    num_ks, num_bands = E_v.shape[0], E_v.shape[1]
    k_v_w, E_v_w, num_bands_w = write_wannier_band()
    print(k_v_w.shape)
    
    #plot the band structure calculated by DFT
    fig = plt.figure(1,figsize=(10,8))
    for i in range(num_bands):
        print(len(k_v), len(E_v[:,i]))
        plt.plot(k_v, E_v[:,i], "black")
    
    #plot the band structure fitting by Wannier90
    for j in range(num_bands_w):
        plt.plot(k_v_w, E_v_w[:,j], "red")
    
    plt.hlines(0,0,4)
    
    #plt.plot(k_v_w, E_v_w)
    k_sym_points = [0.0, 0.505, 0.797, 1.380]
    k_sym_label = ["G", "M", "K", "G"]
    plt.ylim(-10, 10)
    #plt.xlim(0, 1.38)
    plt.ylabel("Energy(eV)")
    plt.xlabel("Kpoints")
    #plt.xticks(k_sym_points, k_sym_label)
    plt.savefig("FeI3_soc.png", dpi=300)
    plt.show()
    
    
if __name__=="__main__":
    plot_band()
