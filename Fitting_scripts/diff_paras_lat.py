# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 10:21:53 2023

@author: 26526
"""
import numpy as np
from math import pi
import os

def write_spin_fdf(path, file_name):
    with open(file_name, "r") as a:
        origin = a.read() 
    
    diff_agl = np.linspace(0, pi, 19, endpoint=True)
    for i in range(len(diff_agl)):        
        Sx = round(np.cos(diff_agl[i]), 6)
        Sy = round(np.sin(diff_agl[i]), 6)
        site1 = "1 1 0 0"
        site2 = "2 " + str(Sx) + " " + str(Sy) + " " + "0"
        new_file = path + "//lat_mode_" + str(i) + ".fdf"
        with open(new_file, "w") as a1:
            print(origin, file=a1)
            print("%block spin_directions_set", file=a1)
            print(site1, file = a1)
            print(site2, file = a1)
            print("%endblock spin_directions_set", file=a1)
        a1.close()      
    a.close()
    print("write spin files is completed !")


def write_lat_fdf(path_name, eg, t2g, po, pdpi, pdsig):
    with open("lattice_model.fdf", "r") as f1:
        origin = f1.read()
    
    file_name = path_name + "//lat_" + str(eg) + "_" + str(t2g) + "_" + str(po) + "_" + str(pdpi) + "_" + str(pdsig)+".fdf"
    with open(file_name, "w") as f:
        # write original values
        print(origin, file=f)
        
        #write d and p orbital informations
        print("%block Orbitals", file=f)
        print("1 5 #is,norb", file=f)
        print(str(eg) + " " + str(eg) + " " + str(t2g) + " " + str(t2g), " " + str(t2g), file=f)
        print("2 3", file=f)
        print(str(po) + " " + str(po) + " " + str(po), file=f)
        print("%endblock Orbitals" + "\n", file=f)
        
        #print hopping parameters
        print("%block Slater_Koster_Hop_Parameters", file=f)
        print("1 2  3.82165     # Cu-O", file=f)
        print("1", file=f)
        print("3 2 2 1", file=f)
        print(str(pdpi) + " " + str(pdsig), file=f)
        print("%endblock Slater_Koster_Hop_Parameters", file=f)
    
    f1.close()
    f.close()
    
    return file_name
        


def main():
    for eg in np.linspace(-4, -2.1, 2, endpoint=True):
        for tg in np.linspace(-6, -4.1, 2, endpoint=True):
            for po in np.linspace(-2.0,-1.1,2, endpoint=True):
                for pdpi in np.linspace(-2, -1.1, 2, endpoint=True):
                    for pdsig in np.linspace( 0.1, 1.0, 2, endpoint=True):
                        eg = round(eg,2)
                        tg = round(tg,2)
                        po = round(po,2)
                        pdpi = round(pdpi,2)
                        pdsig = round(pdsig, 2)
                        print("eg is:", eg, "t2g is:", tg, "po is:", po, "pdpi is:", pdpi, "pdsig is:", pdsig)
                        path_name = "EG_"+str(eg)+"//T2G_"+str(tg)+"//PO_"+str(po)+"//PDP_"+str(pdpi)+"//PDS_"+str(pdsig)    
                        if not os.path.exists(path_name):
                            os.makedirs(path_name)
                            print("the new directory is cerated !")                           
                        file_name = write_lat_fdf(path_name, eg, tg, po, pdpi, pdsig)
                        
                        write_spin_fdf(path_name, file_name)


if __name__=="__main__":
    main()