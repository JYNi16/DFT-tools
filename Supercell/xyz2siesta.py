#!/usr/bin/env python
"""
This program transform xyz file to siesta input coord file.

"""
import sys
class xyz:
    def __init__(self):
        pass
    
    def read_xyz(self,infile):
        #open infile and read info
        xyz_file=open(infile,'r')
        #all_lines=xyz_file.readlines()
        self.nions=int(xyz_file.readline())
        xyz_file.readline()
        self.labels=[]
        self.coord=[]
        for i in range(self.nions):
            line=xyz_file.readline()
            label,x,y,z=line.split()
            x,y,z=map(float,(x,y,z))
            self.labels.append(label)
            self.coord.append((x,y,z))
        #print self.labels
        #print self.coord
        xyz_file.close()        
        return None
        
class siesta:
    def __init__(self):
        pass

    def write_siesta(self,cxyz,posfile,outfile):
        #get info from cxyz(class xyz),and write siesta coord file
        #first parse cxyz
        self.label_diff=[]
        self.labels=cxyz.labels
        self.nions=cxyz.nions
        self.coord=cxyz.coord
        #init self.label_index
        self.label_index=[0]*self.nions
        #print self.label_index
        for i in range(self.nions):
            found=0
            for j in range(len(self.label_diff)):
                if cxyz.labels[i].upper() == self.label_diff[j].upper() :
                    found=1
                    self.label_index[i]=j+1
            #
            if not found:
                self.label_diff.append(cxyz.labels[i])
                self.label_index[i]=len(self.label_diff)
        #
        #print self.label_index
        #print self.label_diff

        
        s_file=open(outfile,'w')
        sys.stdout=s_file
        print "NumberOfAtoms  ",self.nions    
        print "number_of_species  ",len(self.label_diff)
        if posfile:
            pf=open(posfile,'r')
            line=pf.readline()
            a=float(pf.readline().split()[0])
            print "LatticeConstant       1 Ang"
            print "%block LatticeVectors"
            for i in range(3):
                vec=map(float,pf.readline().split())
                vec[0]=vec[0]*a
                vec[1]=vec[1]*a
                vec[2]=vec[2]*a
                print vec[0],vec[1],vec[2]
            print "%endblock LatticeVectors"
            pf.close()
        print "AtomicCoordinatesFormat  NotScaledCartesianAng"
        print "%block AtomicCoordinatesAndAtomicSpecies"
        for i in range(self.nions):
            print '%16.12f %16.12f %16.12f '%(self.coord[i][0], \
                  self.coord[i][1],self.coord[i][2]), \
                  self.label_index[i], \
                  self.labels[i],i+1
        print "%endblock AtomicCoordinatesAndAtomicSpecies"
        s_file.close()
        
def fhelp():
    print "Usage: python xyz2siesta.py xyzfile"            
    print "   or: python xyz2siesta.py xyzfile POSCAR"
def main():
    posfile=''
    if len(sys.argv) == 2:
        xyzfile=sys.argv[1]
    elif len(sys.argv) == 3:
        xyzfile=sys.argv[1]
        posfile=sys.argv[2]
    else:
        fhelp()
        sys.exit()

    #print xyzfile
    cxyz=xyz()
    cxyz.read_xyz(xyzfile)
    csiesta=siesta()
    csiesta.write_siesta(cxyz,posfile,'siesta.ANI')
    
if __name__=="__main__":
    main()
