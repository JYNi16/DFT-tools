# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 15:06:51 2023

@author: 26526
"""
#This script is codding for plot the projection bands 

import numpy as np 
import matplotlib.pyplot as plt

def read_band_data(file_path):
    K, Eng, Orb_weight, p_band = [], [], [], []
    
    with open(file_path, "r") as f:
        datas = f.readlines() 
        k_infos = datas[1].split()
        num_kpoints, num_bands = int(k_infos[-2]), int(k_infos[-1])
        print("kpoints is:", num_kpoints, "num_bands is:", num_bands)
        for i in range(2, len(datas)):
            if ("#" in datas[i] or datas[i].strip() == ""):
                pass 
            else:
                value = [float(s) for s in datas[i].split()]
                K.append(value[0])
                Eng.append(value[1])
                Orb_weight.append(value[2])
                #p_band.append(value[2])
    
    #print("len(K)",len(K), "len=(Eng)",len(Eng), "len(Orb)",len(Orb_weight))
    
    E_bands = np.zeros((num_kpoints, num_bands))
    Orb = np.zeros((num_kpoints, num_bands))
    
    for j in range(num_bands):
        E_bands[:,j] = np.array(Eng[j*num_kpoints:(j+1)*num_kpoints])
        Orb[:,j] = np.array(Orb_weight[j*num_kpoints:(j+1)*num_kpoints])
    
    print(E_bands.shape, Orb.shape)
    
    return np.array(K[:2*num_kpoints]), E_bands, Orb

def plot_band():
    k_points, E_bands, weights = read_band_data("PBAND_SUM_SOC_Te.dat")
    num_bands = E_bands.shape[1]
    
    #plot band structure 
    num_kpoints = int(len(k_points)/2)
    print("num_kpoints is:", num_kpoints)
    fig = plt.figure(1, figsize=(10, 14))
    k_sym_points = [0.0, 0.837, 1.804]
    
    #plot total bands
    plt.subplot(221)
    for i in range(num_bands):
        if i % 2 == 0:
            plt.plot(k_points[:num_kpoints], E_bands[:,i], "Red")
        else:
            plt.plot(k_points[num_kpoints:], E_bands[:,i], "Red")
    plt.xlim(0, 0.527)
    plt.ylim(-0.5, 0.5)
    k_sym_label = ["M", r"$\Gamma$", "K"]
    plt.ylabel("Energy(eV)")
    plt.xlabel("Kpoints")
    plt.xticks(k_sym_points, k_sym_label)
    
    #plot pbands 
    plt.subplot(222)
    for j in range(num_bands):
        if j % 2 == 0:
            plt.scatter(k_points[:num_kpoints], E_bands[:,j], c = weights[:,j], cmap="Blues")
        else:
            plt.scatter(k_points[num_kpoints:], E_bands[:,j], c = weights[:,j], cmap="Blues")
    
    plt.colorbar()  
    plt.xlim(0,0.527)
    plt.ylim(-0.5, 0.5)
    k_sym_label = ["M", r"$\Gamma$", "K"]
    #plt.ylabel("Energy(eV)")
    plt.xlabel("Kpoints")
    plt.xticks(k_sym_points, k_sym_label)
    plt.savefig("MBT_pband.png", dpi=300)
    plt.show()

def plot_sub_band(file_name):
    k_points, E_bands, weights = read_band_data(file_name)
    num_bands = E_bands.shape[1]
    
    num_kpoints = int(len(k_points)/2)
    print("num_kpoints is:", num_kpoints)
    fig = plt.figure(1, figsize=(6, 7))
    k_sym_points = [0.0, 0.837, 1.804]
    
    for j in range(num_bands):
        if j % 2 == 0:
            plt.scatter(k_points[:num_kpoints], E_bands[:,j], c = weights[:,j], cmap="Blues")
        else:
            plt.scatter(k_points[num_kpoints:], E_bands[:,j], c = weights[:,j], cmap="Blues")
    
    plt.colorbar()  
    plt.xlim(0,1.804)
    plt.ylim(-1, 1)
    k_sym_label = ["M", r"$\Gamma$", "K"]
    #plt.ylabel("Energy(eV)")
    band_name = file_name.split(".")[0]
    plt.title(band_name)
    plt.xlabel("Kpoints")
    plt.xticks(k_sym_points, k_sym_label)
    plt.savefig(band_name+".png", dpi=300)
    plt.show()
    
def plot_Bi_Te_band():
    k_points, E_bands, weights = read_band_data("PBAND_SUM_SOC_Bi.dat")
    k_points_Te, E_bands_Te, weights_Te = read_band_data("PBAND_SUM_SOC_Te.dat")
    num_bands = E_bands.shape[1]
    
    #plot band structure 
    num_kpoints = int(len(k_points)/2)
    print("num_kpoints is:", num_kpoints)
    fig = plt.figure(1, figsize=(6, 7))
    k_sym_points = [0.0, 0.316, 0.527]
    
    #plot Bi pbands
    for i in range(num_bands):
        if i % 2 == 0:
            plt.scatter(k_points[:num_kpoints], E_bands[:,i], c = weights[:,i], cmap="Reds")
        else:
            plt.scatter(k_points[num_kpoints:], E_bands[:,i], c = weights[:,i], cmap="Reds")
    plt.colorbar() 
    #plot Te pbands 
    for j in range(num_bands):
        if j % 2 == 0:
            plt.scatter(k_points_Te[:num_kpoints], E_bands_Te[:,j], c = weights_Te[:,j], cmap="Blues")
        else:
            plt.scatter(k_points_Te[num_kpoints:], E_bands_Te[:,j], c = weights_Te[:,j], cmap="Blues")
    
    plt.colorbar()  
    plt.xlim(0, 0.527)
    plt.ylim(-0.5, 0.5)
    k_sym_label = ["M", r"$\Gamma$", "K"]
    #plt.ylabel("Energy(eV)")
    plt.xlabel("Kpoints")
    plt.xticks(k_sym_points, k_sym_label)
    plt.savefig("MBT_pband_Bi_Te.png", dpi=300)
    plt.show()

    

if __name__=="__main__":
    #file_name = ["PBAND_SUM_DW_Bi.dat", "PBAND_SUM_DW_Bi_Px.dat", "PBAND_SUM_DW_Bi_Py.dat", "PBAND_SUM_DW_Bi_Pz.dat"]
    #for file in file_name:
    #    plot_sub_band(file)  
    plot_Bi_Te_band()