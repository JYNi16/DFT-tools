# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:39:15 2023

@author: 26526
"""
import numpy as np 
import matplotlib.pyplot as plt 
import os, glob

#Recipicoal Lattice 
G = np.array([0,0,0])
X = np.array([0.5, 0.0, 0.5])
W = np.array([0.5, 0.25, 0.75])

npoints = 200

#define K-path
kgx = np.linspace(G,X,npoints)
kxw = np.linspace(X,W,npoints)
kwg = np.linspace(W,G,npoints)

##K点相对距离
def Dist(r1,r2):
    return np.linalg.norm(r1-r2)  

lgx=Dist(G,X)
lxw=Dist(X,W)
lwg=Dist(W,G)

lk = np.linspace(0, 1, npoints)
xgx = lgx * lk 
xxw = lxw * lk + xgx[-1]
xwg = lwg * lk + xxw[-1]

kpath = np.concatenate((xgx, xxw, xwg), axis=0)

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
    file_path = "D:\\Magnetic Theory\\Y2V2O7\\magnon\\J1"
    path_1 = "G_X.dat"
    path_2 = "X_W.dat"
    path_3 = "W_G.dat"
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
#plt.ylim(0, 20)
plt.show()