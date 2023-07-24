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
L = np.array([0.5, 0.5, 0.5])

npoints = 200

#define K-path
kgx = np.linspace(G,X,npoints)
kxw = np.linspace(X,W,npoints)
kwg = np.linspace(W,G,npoints)
kgl = np.linspace(G,L,npoints)

##K点相对距离
def Dist(r1,r2):
    return np.linalg.norm(r1-r2)  

lgx=Dist(G,X)
lxw=Dist(X,W)
lwg=Dist(W,G)
lgl=Dist(G,L)

lk = np.linspace(0, 1, npoints)
xgx = lgx * lk 
xxw = lxw * lk + xgx[-1]
xwg = lwg * lk + xxw[-1]
xgl = lgl * lk + xwg[-1]

kpath = np.concatenate((xgx, xxw, xwg, xgl), axis=0)


def read_file(file_path):
    with open(file_path, "r") as f:
        E_v = np.ones((npoints, 4))
        Eng = f.readlines()
        l = len(Eng)
        for i in range(l):
            e = Eng[i].split()
            E_v[i][0] = e[1]
            E_v[i][1] = e[2]
            E_v[i][2] = e[3]
            E_v[i][3] = e[4]
    
    f.close()
    
    return E_v 
    

def read_Ev():
    file_path = "E:\\Magnon\\Y2V2O7\\J1"
    path_1 = "G_X.dat"
    path_2 = "X_W.dat"
    path_3 = "W_G.dat"
    path_4 = "G_L.dat"
    band_path = [path_1, path_2, path_3, path_4]
    
    E_v_k = []
    
    for i in range(len(band_path)):
        print("ith band_path is:", band_path[i])
        filepath = os.path.join(file_path, band_path[i])
        print("file_path is:", filepath)
        E_v = read_file(filepath)
        E_v_k.append(E_v)
    
    return E_v_k 

plt.figure(1, figsize = (7,6))
E = read_Ev() 
E_gx, E_xw, E_wg, E_gl = E[0], E[1], E[2], E[3]
E_1 = np.hstack((E_gx[:,0], E_xw[:,0], E_wg[:,0], E_gl[:,0]))
E_2 = np.hstack((E_gx[:,1], E_xw[:,1], E_wg[:,1], E_gl[:,1]))
E_3 = np.hstack((E_gx[:,2], E_xw[:,2], E_wg[:,2], E_gl[:,2]))
E_4 = np.hstack((E_gx[:,3], E_xw[:,3], E_wg[:,3], E_gl[:,3]))

plt.plot(kpath, E_1)
plt.plot(kpath, E_2)
plt.plot(kpath, E_3)
plt.plot(kpath, E_4)
#plt.ylim(0, 20)
plt.show()
