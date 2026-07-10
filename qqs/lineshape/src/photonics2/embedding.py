# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 14:48:32 2022

@author: 87930
"""
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Poscar, Kpoints
import numpy as np
import sys

def conference_matrix(A,B,d):#小到大对应
    mat=[]
    #print(lenA)
    #print(B)
    for m in range(len(A)):
        dd = 1000000.0
        ref=-1
        
        for i in range(len(B)):
            #print(A[m]+d,B[i])
            dnow=np.dot(A[m]-B[i]+d,A[m]-B[i]+d)
            #print(dnow)
            if dnow<dd:
                dd = dnow
                ref=i
        mat.append(ref)
    return mat
    

class embedding:
    def __init__(self,defect_lattice,perfect_lattice):
        self.defect = Poscar.from_file("{}".format(defect_lattice))
        self.perfect = Poscar.from_file("{}".format(perfect_lattice))
        #self.perfect.write_file("POSCAR.out2")
        
        #print(self.defect.structure.distance_matrix)
    
    def embedd(self,**parameter):
        a1=parameter.get("atom1",-1)-1 
        a2=parameter.get("atom2",-1)-1
        if a1==-1 or  a2==-1:
            print("input error!")
            sys.exit()
        
        else:
            self.d=self.perfect.structure._sites[a2].coords-self.defect.structure._sites[a1].coords
            #print(self.perfect.structure._sites[a2].coords,self.defect.structure._sites[a1].coords,self.d)
            self.defect_coord=np.array([self.defect.structure._sites[i].coords for i in range(len(self.defect.structure._sites))])
            #print(self.defect_coord)
            self.perfect_coord=np.array([self.perfect.structure._sites[i].coords for i in range(len(self.perfect.structure._sites))])

            self.ref_mat=conference_matrix(self.defect_coord,self.perfect_coord,self.d)#得到对应矩阵
            #print(self.ref_mat)
            
            #for i in range(len(self.defect.structure._sites)):
                #print(self.defect_coord[i]+self.d-self.perfect_coord[self.ref_mat[i]] )
            
            delta=0.0
            count=0            
            for i in range(len(self.defect_coord)):

                #print("distance=",sum((self.defect_coord[i]-self.defect_coord[a1])**2))
                if sum((self.defect_coord[i]-self.defect_coord[a1])**2)>parameter.get("rang",4)**2:
                    delta-=self.defect_coord[i]-self.perfect_coord[self.ref_mat[i]]
                    count+=1
            #print("count=",count)
            #print(self.d)
            self.d=delta/float(count)
            #print(self.d)
            
           
            #print(re)
            remove_more=[]
            for i in range(len(self.perfect_coord)):

                #print("distance=",sum((self.defect_coord[i]-self.defect_coord[a1])**2))
                if sum((self.perfect_coord[i]-self.perfect_coord[a2])**2)<parameter.get("rang",4)**2:
                    if i not in self.ref_mat:
                        remove_more.append(i)
            print(remove_more)
            #self.perfect.structure.remove_sites(remove_more)
            #self.perfect.write_file("POSCAR.remove_close.vasp")
            
            self.perfect.structure.remove_sites(self.ref_mat+remove_more)
            self.perfect.write_file("POSCAR.remove_defect_cell.vasp")
                    

            for i in range(len(self.defect.structure._sites)):
                frac_coord=np.linalg.solve(self.perfect.structure.lattice.matrix,self.defect_coord[i]+self.d)
                #print(self.defect.structure.species[i],self.defect_coord[i]+self.d)
                self.perfect.structure.append(self.defect.structure.species[i] , frac_coord)
          
            self.perfect.write_file("POSCAR.embedded.vasp")

    def read_force_constants(self,file):
        if "0" == file:
            print("error")
        else:
            f=open(file,"r")
            natom=int(f.readline().split()[0])
            print(natom)
            force_constants=np.zeros(natom*natom*9).reshape((natom,natom,3,3))
            for n in range(natom*natom):
                line=f.readline().split()
                i,j=int(line[0])-1,int(line[1])-1
                for a in range(3):
                    line = f.readline().split()
                    for b in range(3):
                        force_constants[i][j][a][b] = float(line[b])
            f.close()
            return force_constants
        
    def write_force_constants(self,force_constants,file,natom):
        if "0" == file:
            print("error")
        else:
            f=open(file,"w")
            f.write(str(natom)+"\t"+str(natom)+"\n")
            for i in range(natom):
                for j in range(natom):
                    f.write(str(i+1)+"\t"+str(j+1)+"\n")
                    for a in range(3):
                        for b in range(3):
                            f.write(str(force_constants[i][j][a][b])+"\t")
                        f.write("\n")
            f.close()
            
        
    
    def create_force_constant(self,r1,r2,ndefect_defect,ndefect_large,**parameter):
        self.single_force=self.read_force_constants(parameter.get("single_force_constant","0"))
        self.defect_force=self.read_force_constants(parameter.get("defect_force_constant"",0)) 
        self.single = Poscar.from_file("{}".format(perfect_lattice))
        self.large = Poscar.from_file("{}".format(perfect_lattice))
        self.defect = Poscar.from_file("{}".format(perfect_lattice))
        force_constants=np.zeros(self.large.natoms**2*9).reshape((self.large.natoms,self.large.natoms,3,3))
            
        ref=conference_matrix()
        for i in range(self.large.natoms):
            for j in range(self.large.natoms):
                if distance(i,j)>r2:
                    continue
                elif distance(i,natoms) <r1 or distance(j,natoms) <r1:
                    force_constants[i][j]=self.defect_force[ref.find(i)][ref.find(j)]
                else:
                    #find force in singlecell
                
        
        
        
        self.write_force_constants(force_constants,parameter.get("large_force_constant","0"),511)
        
        #self.large_force = force_constants=np.zeros(self.large.natoms*self.large.natoms*3*3).reshape((self.large.natoms,self.large.natoms,3,3))
        

            
            
            

            
            
                    
            
            


        
        
        

    