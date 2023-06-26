# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:39:15 2023

@author: 26526
"""
import numpy as np 
import matplotlib.pyplot as plt 
import os, glob

#Recipicoal Lattice 
G = np.array([0,0])
K = np.array([1/3, 1/3])
M = np.array([0.5, 0])

npoints = 200

#define K-path
kmg = np.linspace(M,G,npoints)
kgk = np.linspace(G,K,npoints)
kkm = np.linspace(K,M,npoints)

##K点相对距离
def Dist(r1,r2):
    return np.linalg.norm(r1-r2)  

lmg=Dist(M,G)
lgk=Dist(G,K)
lkm=Dist(K,M)

lk = np.linspace(0, 1, npoints)
xmg = lmg * lk 
xgk = lgk * lk + xmg[-1]
xkm = lkm * lk + xgk[-1]

kpath = np.concatenate((xmg, xgk, xkm), axis=0)

def read_file(file_path):
    with open(file_path, "r") as f:
        E_v = np.ones((npoints, 2))
        Eng = f.readlines()
        l = len(Eng)
        for i in range(l):
            e = Eng[i].split()
            E_v[i][0] = e[1]
            E_v[i][1] = e[2]
    
    f.close()
    
    return E_v 
    

def read_Ev():
    file_path = "E:\\magnon\\2D_honeycomb_FM\\lattice_model\\Kitaev"
    path_1 = "M_G.dat"
    path_2 = "G_K.dat"
    path_3 = "K_M.dat"
    band_path = [path_1, path_2, path_3]
    
    E_v_k = []
    
    for i in range(len(band_path)):
        print("ith band_path is:", band_path[i])
        filepath = os.path.join(file_path, band_path[i])
        print("file_path is:", filepath)
        E_v = read_file(filepath)
        E_v_k.append(E_v)
    
    return E_v_k 

E = read_Ev() 
E_mg, E_gk, E_km = E[0], E[1], E[2]
E_up = np.hstack((E_mg[:,0], E_gk[:,0], E_km[:,0]))
E_dn = np.hstack((E_mg[:,1], E_gk[:,1], E_km[:,1]))

plt.plot(kpath, E_up)
plt.plot(kpath, E_dn)
plt.ylim(0, 20)
plt.show()
        