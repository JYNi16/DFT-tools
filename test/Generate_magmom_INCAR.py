# -*- coding: utf-8 -*-
"""
### Generating the INCAR with SOC calculations 

@author: Curry
"""

import numpy as np
from math import pi
import os

def write_INCAR(S_m, state_num):
    with open("INCAR", "r") as a:
        origin = a.read() 
    
    new_file = "INCAR_" + str(state_num)
    
    with open(new_file, "w") as a2:
        print(origin, file=a2)
        print("MAGMOM = " + S_m, file = a2)
        print("M_CONSTR = " + S_m, file = a2)
    
    print("write INCAR is completed !")


def float_str(s):
    
    t = ''
    for i in s:
        t += str(i) + " "
    
    return t

def four_state(state, s1, s2, so, n1, n2, N_m, N_nm):
    
    S_m = ""
    
    if state == 2:
        s1 = [-x if x != 0 else x for x in s1]
    elif state == 3:
        s2 = [-x if x != 0 else x for x in s2] 
    elif state == 4:
        s1 = [-x if x != 0 else x for x in s1]
        s2 = [-x if x != 0 else x for x in s2]
    else:
        s1 = s1 
        s2 = s2
    
    for i in range(1, N_m+1): 
        
        if i == n1: 
            S_m += float_str(s1)
        elif i == n2:
            S_m += float_str(s2)
        else:
            S_m += float_str(so)
    
    S_m += str(N_nm*3) + "*0"
    
    return S_m

def main(n1, n2, N_m, N_nm, S):
    
    dxyz1, dxyz2, dxyz3 = 0, 1, 2
    
    S_xyz = [[S, 0, 0], [0, S, 0], [0, 0, S]]
    #S_xyz = [str[i] for i in Sxyz]
    
    s1 = S_xyz[dxyz1]  
    s2 = S_xyz[dxyz2]
    so = S_xyz[dxyz3]
    
    for state_num in [1, 2, 3, 4]:
        S_m = four_state(state_num, s1, s2, so, n1, n2, N_m, N_nm)
        print("s1 is:", S_m)
        
        write_INCAR(S_m, state_num) 


if __name__=="__main__":
    
    main(1, 3, 6, 10, 5)