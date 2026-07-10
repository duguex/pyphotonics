import yaml
import numpy as np
from scipy import constants as constant
from photonics2.configuration_coordinate import ConfigurationCoordinate
import copy
import matplotlib.pyplot as plt
import os

class lattice_inform:
    def __init__(self,lattice,natom,points):
        self.lattice_vector=np.array(lattice)
        self.natom=natom
        self.atom=[]
        self.mass=[]
        self.atom_frac_position=[]
        for a in points:
            self.atom.append(a['symbol'])
            self.mass.append(a['mass'])
            self.atom_frac_position.append(a['coordinates'])
        self.mass=np.array(self.mass)
        self.atom_frac_position=np.array(self.atom_frac_position)
        self.atom_position=[]
        for a in range(self.natom):
            pos=self.atom_frac_position[a][0]*self.lattice_vector[0]+self.atom_frac_position[a][1]*self.lattice_vector[1]+self.atom_frac_position[a][2]*self.lattice_vector[2]
            self.atom_position+=[pos]#*self.lattice_vector[0]
        c=self.atom[0]
        self.formula={}
        count=0
        for a in self.atom:
            if c == a:
                count+=1
            else:
                self.formula[c]=count
                c=a
                count=1
        self.formula[c]=count
        
            
class phonon_point:
    def __init__(self,data):
        self.frequency=data['frequency']
        self.eigenvector=[]
        self.eigenvector_real=[]
        for each in data['eigenvector']:
            self.eigenvector+=[[complex(a[0],a[1]) for a in each]]
            self.eigenvector_real+=[[a[0] for a in each]]
        self.eigenvector=np.array(self.eigenvector)
        self.eigenvector_real=np.array(self.eigenvector_real)
        #print(self.eigenvector)
            
        

class phonon:
    def __init__(self,position,band):
        self.position=position
        self.phonon=[]
        for p in band:
            self.phonon+=[phonon_point(p)] 
    
    
class struct:
    def __init__(self,**para):
        cc = ConfigurationCoordinate()
        if para.get("phonon",0)!=0:
            f=open(para.get("phonon",0),"r")
            self.data=yaml.load(f,Loader=yaml.CLoader)
            self.lattice=lattice_inform(self.data['lattice'], self.data['natom'], self.data['points'])
            self.phonon_all_point=[phonon(a['q-position'],a['band']) for a in self.data['phonon']]
            if para.get("POSCAR",0)!=0:
                cc = ConfigurationCoordinate()
                self.poscar=cc.read_poscar(para.get("POSCAR",0))
                self.lattice.lattice_vector=np.array(self.poscar.lattice.matrix)
                self.lattice.atom_frac_position=np.array(self.poscar.frac_coords)
                self.lattice.atom_position=np.array(self.poscar.lattice.get_cartesian_coords(self.poscar.frac_coords))
                
            f.close()#only first phonon mode for now
        if para.get("POSCAR1",0)!=0:
            self.poscar1=cc.read_poscar(para.get("POSCAR1",0))
        if para.get("POSCAR2",0)!=0:
            self.poscar2=cc.read_poscar(para.get("POSCAR2",0))
        if para.get("POSCAR3",0)!=0:
            self.poscar3=cc.read_poscar(para.get("POSCAR3",0))
        if para.get("POSCAR4",0)!=0:
            self.poscar4=cc.read_poscar(para.get("POSCAR4",0))
    def fold(self,d):
        nml= [ round(a) for a in np.linalg.solve(self.poscar1.lattice.matrix,d)]
        #print(d,nml)
        return d-np.matmul(self.poscar1.lattice.matrix,nml)
        
    def visual_mode(self,k_points,n,**para):
        f=open("./mode_eigenvector_data","w")
        f.write("system:"+str(self.lattice.formula))
        f.write("\nq-position:"+str(self.phonon_all_point[k_points].position)+"  "+str(n)+"th mode\n")
        f.write("atom\tx\ty\tz\tdx\tdy\tdz\n")
        mode=self.phonon_all_point[k_points].phonon[n].eigenvector
        for i in range(len(self.lattice.atom)):
            f.write(str(self.lattice.atom[i])+"\t")
            for j in range(3):
                f.write(str(self.lattice.atom_position[i][j])+"\t")
            for j in range(3):
                f.write(str(mode[i][j].real*0.1)+"\t")
            f.write("\n")
            
            
        
        f=open(para.get("outdir","")+"MODE"+str(n+1)+".xsf","w")
        f.write("CRYSTAL\nPRIMVEC\n")
        for a in self.lattice.lattice_vector:
            for b in a:
                f.write("\t"+str(b))
            f.write("\n")
        f.write('PRIMCOORD\n')
        f.write(str(self.lattice.natom)+"\n")


        for i in range(self.lattice.natom):
            f.write(str(self.lattice.atom[i])+"\t")
            for j in range(3):
                f.write(str(self.lattice.atom_position[i][j])+"\t")
            for j in range(3):
                f.write(str(self.phonon_all_point[k_points].phonon[n].eigenvector_real[i][j])+"\t")
            f.write("\n")
        f.close
        print("phonon No.",n+1," freq=",self.phonon_all_point[k_points].phonon[n].frequency)
            
        
    def visual_DR(self,**para):
        cc = ConfigurationCoordinate()
        k=para.get("k",1.0)
        poscar2=cc.read_poscar(para.get("POSCAR2",0))
        poscar1=cc.read_poscar(para.get("POSCAR1",0))
        d=(poscar2.cart_coords-poscar1.cart_coords)
        #print(d)
        D_R=[]
        for i in range(len(d)):
            D_R.append(self.fold(d[i]).tolist())
        d=k*np.array(D_R)
        
        f=open(para.get("outdir","")+"DR.xsf","w")
        f.write("CRYSTAL\nPRIMVEC\n")
        for a in self.lattice.lattice_vector:
            for b in a:
                f.write("\t"+str(b))
            f.write("\n")
        f.write('PRIMCOORD\n')
        f.write(str(self.lattice.natom)+"\n")


        for i in range(self.lattice.natom):
            f.write(str(self.lattice.atom[i])+"\t")
            for j in range(3):
                f.write(str(self.poscar1.cart_coords[i][j])+"\t")
            for j in range(3):
                f.write(str(d[i][j])+"\t")
            f.write("\n")
        f.close
            

        
    def move_thorough_phonon(self,k,k_of_mode,k2,k_of_mode2,**para):
        
        d=k*self.phonon_all_point[0].phonon[k_of_mode].eigenvector_real
        d+=k2*self.phonon_all_point[0].phonon[k_of_mode2].eigenvector_real
              
        moved=copy.deepcopy(self.poscar1)
        for i in range(len(moved.cart_coords)):
            moved.replace(i, moved.species[i],moved.cart_coords[i]+d[i],True)
            
        q=0
        for i in range(self.lattice.natom):
            
            q+=np.dot(d[i],d[i])*self.lattice.mass[i]
    
        print((q)**0.5,"A*amu^0.5")
         
        return moved,(q)**0.5
        """
        f=open(para.get("name","out"),"w")
        f.write("moved deltaQ="+str(q/constant.physical_constants["atomic mass constant"][0])+"\n1.0\n")
        for a in self.lattice.lattice_vector:
            for b in a:
                f.write("\t"+str(b))
            f.write("\n")
        for key in self.lattice.formula:
            f.write(key+"  ")
        f.write("\n")
        for key in self.lattice.formula:
            f.write(str(self.lattice.formula[key])+"  ")
        f.write("\nCart\n")
        
        for a in moved_lattice.tolist():
            for b in a:
                #print(b)
                f.write(str(b)+"\t")
            f.write("\n")
        f.close()
        """
       
    
def mix_two_poscar(sel,k,**para):
    d=(sel.poscar2.cart_coords-sel.poscar1.cart_coords)
    #print(d)
    D_R=[]
    for i in range(len(d)):
        D_R.append(sel.fold(d[i]).tolist())
    d=k*np.array(D_R)

    moved=copy.deepcopy(sel.poscar3)
    for i in range(len(moved.cart_coords)):
        moved.replace(i,moved.species[i],moved.cart_coords[i]+d[i],True)

        
    moved.to(filename=para.get("name","out"))
        
    q=0
    for i in range(len(sel.poscar1.cart_coords)):
        q+=np.dot(d[i],d[i])
    print((q)**0.5) 

    return (q)**0.5
    
    
def mix_two_poscar_for_round(sel,theta,mass,**para):
    d1=(sel.poscar2.cart_coords-sel.poscar1.cart_coords)
    #print(d)
    D_R=[]
    for i in range(len(d1)):
        D_R.append(sel.fold(d1[i]).tolist())
    d1=np.array(D_R)

    d2=(sel.poscar3.cart_coords-sel.poscar1.cart_coords)
    #print(d)
    D_R=[]
    for i in range(len(d2)):
        D_R.append(sel.fold(d2[i]).tolist())
    d2=np.array(D_R)
    
    d=d1*np.cos(theta)+d2*np.sin(theta)
    
    
    moved=copy.deepcopy(sel.poscar4)
    for i in range(len(moved.cart_coords)):
        moved.replace(i,moved.species[i],moved.cart_coords[i]+d[i],True)

        
    moved.to(filename=para.get("name","out"))
        
    q=0
    for i in range(len(sel.poscar1.cart_coords)):
        q+=np.dot(d[i],d[i])*mass[i]
    print((q)**0.5) 

    return (q)**0.5        
   
def allmode(a,outdire):
                
    
    print(a.lattice.formula)
    print(a.lattice.natom)
    for i in range(3*a.lattice.natom):
        a.visual_mode(0,i,outdir=outdire)
        
        
        
       


masscsgo=[9.12267699651288e-26, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26]
masskso=[9.12267699651288e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 6.492425458764678e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26]
masskso2=copy.deepcopy(masscsgo)
massbso=[9.12267699651288e-26, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 2.280368483989782e-25, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26]

massvsi=[]
for i in range(215):
    massvsi+=[28.085500]
for i in range(64):
    masskso[i+1]=5.324518517052899e-26



#=struct(phonon="./CsSO-pbesol/band.ex.nU.yaml")
#allmode(a,"./CsSO-pbesol/ex-noU-mode/")

"""  

a=struct(phonon="./CsSO-pbesol/band.ex.yaml",POSCAR1="./CsSO-pbesol/zx-CONTCAR")
a.visual_DR(outdir="./CsSO-pbesol/ZX-EXHS-",POSCAR1="./CsSO-pbesol/zx-CONTCAR",POSCAR2="./CsSO-pbesol/ex-hs")
a.visual_DR(outdir="./CsSO-pbesol/YZ-EXHS-",POSCAR1="./CsSO-pbesol/YZ-CONTCAR",POSCAR2="./CsSO-pbesol/ex-hs")
a.visual_DR(outdir="./CsSO-pbesol/x2-y2-EXHS-",POSCAR1="./CsSO-pbesol/x2-y2-CONTCAR",POSCAR2="./CsSO-pbesol/ex-hs")

a=struct(phonon="./CsSO-pbesol/band.yaml",POSCAR1="./CsSO-pbesol/zx-CONTCAR")
a.visual_DR(outdir="./CsSO-pbesol/ZX-ls-",POSCAR1="./CsSO-pbesol/zx-CONTCAR",POSCAR2="./CsSO-pbesol/ls")
a.visual_DR(outdir="./CsSO-pbesol/YZ-ls-",POSCAR1="./CsSO-pbesol/YZ-CONTCAR",POSCAR2="./CsSO-pbesol/ls")
a.visual_DR(outdir="./CsSO-pbesol/x2-y2-ls-",POSCAR1="./CsSO-pbesol/x2-y2-CONTCAR",POSCAR2="./CsSO-pbesol/ls")



#a=struct(phonon="./CsSO-pbesol/band.yaml")
#allmode(a,"./CsSO-pbesol/mode/")phonon="./CsSO-pbesol/band.yaml",


    


#mix_two_poscar(a,0.27,name="./CsSO-pbesol/POSCAR-LS3")

#mass=[9.12267699651288e-26, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 2.206946952101311e-25, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 5.324518517052899e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26, 2.656762874216004e-26]

index=0
for k in [0.25,0.5,0.75,1,1.25,1.5,1.75]:
    index+=1
    mix_two_poscar(a,k,mass,name="./CsSO-pbesol/vertical/POSCAR-"+str(index))
for k in [0.25,0.5,0.75,1,1.25,1.5,1.75]:
    index+=1
    mix_two_poscar(a,-k,mass,name="./CsSO-pbesol/vertical/POSCAR-"+str(index))
 


    

 
#a=struct(POSCAR1="C:/Users/权照勇/Desktop/fsdownload/POSCAR",\
         #POSCAR2="C:/Users/权照勇/Desktop/fsdownload/CONTCAR",POSCAR3="C:/Users/权照勇/Desktop/fsdownload/POSCAR") 
a=struct(phonon="./NMD/band.yaml",POSCAR1="./NMD/POSCAR")
allmode(a,"./NMD/mode/")

for i in range(3):
    poscar_out=a.move_thorough_phonon(-0.00, 48,0.01*(i+1), 49)   
    poscar_out.to(filename="./NMD/POSCAR-"+str(i+1)) 


#a=struct(POSCAR1="./NMD/POSCAR",\
         #POSCAR2="./NMD/CONTCAR",POSCAR3="./NMD/POSCAR") 
   
#mix_two_poscar(a,-1,masskso2,name="./NMD/POSCAR.ANT")   
    

kk=[0.2,0.4,0.6,0.8,0.9,1,1.1,1.2,1.4,1.6,1.8,2] 
index=0
for k in range(len(kk)):
    index+=1
    mix_two_poscar(a,-kk[-k-1],masskso2,name="./BaSO/ls_path/POSCAR-"+str(index))
for k in range(len(kk)):
    index+=1
    mix_two_poscar(a,kk[k],masskso2,name="./BaSO/ls_path/POSCAR-"+str(index))
 
  
#allmode("./BaSO/band.yaml","./BaSO/mode/")


#mix_two_poscar_for_round(a,5/4*np.pi,massbso,name="BaSO/POSCAR")


#mix_two_poscar(a,0.5,masskso2,name="./CsSO-pbesol/POSCAR-halfx2-y2") 










a=struct(phonon="./CsSO-pbesol/band.yaml",POSCAR1="./CsSO-pbesol/hs-noU") 
#allmode(a,"./KSO-pbesol/EXmode/")
kk=[0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2]
kk=[(i+1)*0.1 for i in range(20)]
#kk=[0.1,0.3,0.5,0.7,1.3,1.5,1.7,1.9]
#print(kk2,kk)

index=0
for k in range(len(kk)):
    index+=1
    #q=a.mix_two_poscar(k,name="./smallu/POSCAR-"+str(index))
    poscar_out=a.move_thorough_phonon(-kk[-k-1]*0.045, 384,0.0, 385)
    poscar_out.to(filename="./CsSO-pbesol/nUQx/POSCAR-"+str(index))
    

for k in range(len(kk)):
    index+=1
    #q=a.mix_two_poscar(k,name="./smallu/POSCAR-"+str(index))
    poscar_out=a.move_thorough_phonon(kk[k]*0.045, 384,0.0, 385)
    poscar_out.to(filename="./CsSO-pbesol/nUQx/POSCAR-"+str(index))

index=0
for k in range(len(kk)):
    index+=1
    #q=a.mix_two_poscar(k,name="./smallu/POSCAR-"+str(index))
    poscar_out=a.move_thorough_phonon(0.0, 384,-(kk[-k-1]-0.1)*0.045, 385)
    poscar_out.to(filename="./CsSO-pbesol/nUQy/POSCAR-"+str(index))
    

for k in range(len(kk)):
    index+=1
    #q=a.mix_two_poscar(k,name="./smallu/POSCAR-"+str(index))
    poscar_out=a.move_thorough_phonon(0.0, 384,(kk[k]-0.1)*0.045, 385)
    poscar_out.to(filename="./CsSO-pbesol/nUQy/POSCAR-"+str(index))

    #f.write("POSCARtest-"+str(index)+"\tq=\t"+str(q)+"\n")
#f.close()







kk=[1.2,1.4,1.6,1.8,2.0,2.5,3.0] 
index=100
for k in kk:
    index+=1
    a=struct(POSCAR1="./CsSO-pbesol/ex-hs-noU",POSCAR2="./CsSO-pbesol/yz-noU-CONTCAR",POSCAR3="./CsSO-pbesol/ex-hs-noU")
    print(k)
    mix_two_poscar(a,k,masscsgo,name="./CsSO-pbesol/yz-ex-hs-noU/POSCAR-"+str(index)) 
"""



"""
a=struct(phonon="./CsSO-pbesol/band.ex.nU.yaml",POSCAR1="./CsSO-pbesol/ex-hs-noU") 
#allmode(a,"./KSO-pbesol/EXmode/")
kk=np.array([0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2])*0.045*4


for kx in range(len(kk)):
    print(kx)
    poscar_out=a.move_thorough_phonon(kk[kx], 384, 0, 385)
    poscar_out.to(filename="./CsSO-pbesol/yz_DMS/parallel/POSCAR-"+str(kx+10))
    
  
for kx in range(len(kk)):
    print(kx)
    poscar_out=a.move_thorough_phonon(0, 384, kk[kx], 385)
    poscar_out.to(filename="./CsSO-pbesol/yz_DMS/vertical/POSCAR-"+str(kx+10))  


for kx in range(len(kk)):
    for ky in range(len(kk)):
        print(kx,ky)
        poscar_out=a.move_thorough_phonon(-kk[kx]-0.4*0.045*4, 384, kk[ky]-1.1*0.045*4, 385)
        poscar_out.to(filename="./CsSO-pbesol/x2-y2_DMS/POSCAR-"+str(kx)+"_"+str(ky))

   






""" 

def fmix(x,y,pm,hwx,hwy): 
    kxx=-0.624
    kyy=-0.624
    kxy=-0.533
    a3=-0.0125
    g=-0.037
    H=0.005
    Sx=0.005
    Sy=0.0
    
    kxx=-0.575
    kyy=-0.575
    kxy=-0.455
    a3=-0.0125
    g=-0.037
    H=0.005
    Sx=0.0092
    Sy=0.0
    a=hwx*kxx*x+Sx+hwx*g*(x**2-y**2)+hwx*H*(x**3+x*y**2)
    b=hwy*kxy*y+Sy-2*hwx*g*x*y+hwx*H*(y**3+y*x**2)

    A=np.array([(-a+pm*np.linalg.norm([a,b]))/-b,1])
    A=A/np.linalg.norm(A)
    return  A[0]**2


def solve(y,a,st=[-5,5]):
    
    if st[1]-st[0]<1e-7:
        return st[0]
    
    
    down=fmix(st[0],y,-1,0.04083,0.04088) -a
    up=fmix(st[1],y,-1,0.04083,0.04088) -a
    mid=fmix(st[0]/2+st[1]/2,y,-1,0.04083,0.04088) -a
    #print(st[0],down,st[0]/2+st[1]/2,mid,st[1],up)
    
    
    if down*up >0 :
        return 1e10
    elif down*mid <0:
        return solve(y,a,st=[st[0],st[0]/2+st[1]/2])
    else:
        return solve(y,a,st=[st[0]/2+st[1]/2,st[1]])
    


    
def f(x,y,pm,hwx,hwy): 
    kxx=-0.624
    kyy=-0.624
    kxy=-0.533
    a3=-0.0125
    g=-0.030
    H=0.005
    Sx=0.005
    Sy=0.0
    
    kxx=-0.575
    kyy=-0.575
    kxy=-0.455
    a3=-0.0125
    g=-0.030
    H=0.005
    Sx=0.0092
    Sy=0.0
    
    
    a=hwx*kxx*x+Sx+hwx*g*(x**2-y**2)+hwx*H*(x**3+x*y**2)
    b=hwy*kxy*y+Sy-2*hwx*g*x*y+hwx*H*(y**3+y*x**2)
    c=-a
    E=0.5*(a+c+pm*np.sqrt( (a+c)**2-4*(a*c-b**2) ))
    
    return E+hwx*0.5*x**2+hwy*0.5*y**2+hwx*a3*(x**3-3*x*y**2)

   
def mix_tutu(hwx,hwy,pp=np.array([])):
    fig=plt.figure(figsize=(10, 10),dpi=300)
       
    Q=4.07497E-24
    hbar=1.05457E-34
    e=1.60218E-19

    Qx=Q*(hwx*e)**0.5/hbar
    Qy=Q*(hwy*e)**0.5/hbar
    
    n = 200
    x = np.linspace(-0.4,0.4,n)
    y = np.linspace(-0.4,0.4,n)
    X,Y = np.meshgrid(x,y)
    #print(X,Y)
    
    Z=copy.deepcopy(X)
    
    for i in range(len(X)):
        for j in range(len(X[0])):
            Z[i][j]=fmix(X[i][j]*Qx,Y[i][j]*Qy,-1,hwx,hwy)
    
    poscar_x=[]
    poscar_y=[]
    dx=0.144906
    dy=0.144182
    dxe=0.1451144287559751/4
    dye=0.2883422597067737/8
            


    for i in range(9):
        for j in range(14):
            poscar_x+=[(j-7)*dxe]
            poscar_y+=[(i)*dye]
            a=0.5*x+3**0.5*0.5*y
            b=-0.5*y+3**0.5*0.5*x



    #print(poscar_x)
    
    ax=fig.add_subplot(1,1,1)
    #ax.contourf(X, Y, Z, 30, alpha=.75, cmap='jet')
    
    ax.contourf(X, Y, (f(X*Qx,Y*Qy,1,hwx,hwy)-f(X*Qx,Y*Qy,-1,hwx,hwy))**2, 30, alpha=.75, cmap='jet')
    #ax.contourf(X, Y, Z, 20, alpha=.75, cmap='jet')
    cs=ax.contour(X, Y, Z, 20, colors='black')
    ax.clabel(cs, inline=True, fontsize=10)
    #temp=np.array([[ -0.008005582898720753 , 0.04],[ -0.08397542892762837 , 0.08],[ -0.16717013754059332 , 0.12],[ -0.2595447091976949 , 0.16],[ -0.3632429250908367 , 0.2],[ -0.48237316589226403 , 0.24],[ -0.6229828143671975 , 0.28],[ -0.796206734019853 , 0.32],[ -1.0293861342305537 , 0.36],[ -1.4375492744753235 , 0.4]])
    #ax.scatter(poscar_x,poscar_y,color='y')
    ax.scatter(pp.T[0],pp.T[1],color='r')
    
    ax.set_xlabel("$Q_X \ (\sqrt{amu}$Å) ")
    ax.set_ylabel("$Q_Y \ (\sqrt{amu}$Å) ")
    ax.set_title("APS of $E \otimes e \ \epsilon  -$ (eV)")
    
    ax.set_xlim(-0.4,0.4)
    ax.set_ylim(-0.4,0.4)
    
    fig.show()
    return    

def dic():
    Q=4.07497E-24
    hbar=1.05457E-34
    e=1.60218E-19
    hwx=0.04083
    hwy=0.04083
    Qx=Q*(hwx*e)**0.5/hbar
    Qy=Q*(hwy*e)**0.5/hbar
    aa=[0.995,0.975,0.94,0.9,0.85,0.77,0.65,0.5,0.35,0.2,0.1,0.05,0.02,0.005]
    yy=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6])**0.7
    
    p=[]
    dic={}
    fac=[0.05, 0.11, 0.164, 0.208, 0.23 , 0.255, 0.265, 0.26, 0.245, 0.22, 0.182, 0.15, 0.105, 0.058]
    #fac=[0.055, 0.12, 0.19, 0.235, 0.26 , 0.27, 0.275, 0.275, 0.265, 0.245, 0.212, 0.16, 0.100, 0.050]
    for i in range(len(aa)):
    
        a=round(aa[i],3)
        dici=[]
        x=round(solve(0.2*Qy,a)/Qx,5)
        y=0.2
        for y in yy*fac[i]:
            x=round(solve(y*Qy,a)/Qx,5)
            y=round(y,5)
            #print("[",x,",",y,end="],")
            p+=[[x,y]]
            dici+=[[x,y]]
            dic[a]=dici
    
    p=np.array(p)
    print(fac)
    print(dic)
    #mix_tutu(0.04083,0.04088,pp=p)
    mix_tutu(hwx,hwy,pp=p)
    return dic

def E():
    al=dic()
    
    #a=struct(phonon="./KSO-pbesol/band.yaml",POSCAR1="./KSO-pbesol/hs") 
    #allmode(a,"./CsSO-pbesol/mix/")
    Qxtokx=-4.028173837143261
    Qytoky=4.0029534632853875
    
    #Qxtokx=-4.028173837143261/0.36923*0.3722800297898284 	 
    #Qytoky=4.0029534632853875/0.29459*0.29662248223503157
    #print(Qxtokx,Qytoky)
    #0.05
    
    
    for key in al:
        #os.mkdir("./KSO-pbesol/mix_DMS/"+str(key))
        print("./CsSO-pbesol/mix_DMS\t/"+str(key))
        for i in range(len(al[key])):
            print(al[key][i][0],"\t",al[key][i][1])
            #poscar_out=a.move_thorough_phonon(-al[key][i][0]/Qxtokx, 384, al[key][i][1]/Qytoky, 385)
            #poscar_out.to(filename="./KSO-pbesol/mix_DMS/"+str(key)+"/POSCAR-"+str(i))
            #print("-------------")
           
    return


#dic()
  
"""          
a=struct(phonon="./K2TiF6/band.yaml",POSCAR1="./K2TiF6/2Eg-CONTCAR") 
allmode(a,"./K2TiF6/mode/" )
freq_eV=0.004135665538536#ev=freq*freq_eV
Q=4.0749709938E-24#m*kg^0.5
hbar=1.05E-34
e=1.60E-19
for i in range(216): #range(len(a.phonon_all_point[0].phonon)):
    index=0
    ii=i
    print("phonon",i+1)
    for k in np.array([1.2,3.5,7])/(a.phonon_all_point[0].phonon[ii].frequency/freq_eV*1e-3):
        index+=1
        poscar_out,q=a.move_thorough_phonon(k/40, ii,0,1)
        Ehw=0.5*(q*Q)**2*(a.phonon_all_point[0].phonon[ii].frequency*freq_eV/hbar*1.60E-19)**2/1.60E-19
        print("E=",Ehw/0.0023061461027758497*0.5092188321318761*85.107,"meV")
        poscar_out.to(filename="./K2TiF6/DMS/POSCAR-"+str(ii+1)+"_"+str(index))
        
 
a=struct(POSCAR1="./CsSO-pbesol/ex-hs",POSCAR2="./CsSO-pbesol/yz-CONTCAR",POSCAR3="./CsSO-pbesol/ex-hs") 
for i in range(11):
    mix_two_poscar(a,i/10,name="./CsSO-pbesol/new/exhs_to_ex/POSCAR-"+str(i))
"""        

"""
a=struct(phonon="./KSO-pbesol/band.yaml",POSCAR1="./KSO-pbesol/ex-hs") 
#allmode(a,"./KSO-pbesol/EXmode/")

kk=np.array([0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2])*0.045*2

kkX=np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])*0.045*2*4
kkY=np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])*0.045*2*4
for kx in range(len(kkX)):
    for ky in range(len(kkY)):
        print(kkX[kx],kkY[ky])
        poscar_out,q=a.move_thorough_phonon(-kkX[kx], 384, -kkY[ky], 385)
        print(q)
        poscar_out.to(filename="./KSO-pbesol/new/test-ex/POSCAR-"+str(kx)+"_-"+str(ky))

for kx in range(len(kkX)):
    for ky in range(len(kkY)):
        print(kkX[kx],kkY[ky])
        poscar_out,q=a.move_thorough_phonon(kkX[kx], 384, -kkY[ky], 385)
        print(q)
        poscar_out.to(filename="./KSO-pbesol/new/test-ex/POSCAR--"+str(kx)+"_-"+str(ky))


for kx in range(len(kkX)):
    for ky in range(len(kkY)):
        print(kkX[kx],kkY[ky])
        poscar_out,q=a.move_thorough_phonon(-kkX[kx], 384, kkY[ky], 385)
        print(q)
        poscar_out.to(filename="./KSO-pbesol/new/test-ex/POSCAR-"+str(kx)+"_"+str(ky))

for kx in range(len(kkX)):
    for ky in range(len(kkY)):
        print(kkX[kx],kkY[ky])
        poscar_out,q=a.move_thorough_phonon(kkX[kx], 384, kkY[ky], 385)
        print(q)
        poscar_out.to(filename="./KSO-pbesol/new/test-ex/POSCAR--"+str(kx)+"_"+str(ky))

"""
    
    
    
a=struct(phonon="./zlq/band.ex2.yaml") 

allmode(a,"./zlq/modes_ex2")