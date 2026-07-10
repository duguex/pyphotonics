import numpy as np
from numpy import fft
from photonics2.xyz import XYZ
import sys
from photonics2.configuration_coordinate import ConfigurationCoordinate
from scipy import constants as constant
import matplotlib.pyplot as plt
import yaml
import multiprocessing as mp
import math
import copy

def realfft(func,emax):
    return fft.fft(func)/float(len(func))*emax

def realifft(func,emax):
    return fft.ifft(func)*emax

def read_force(file):
    #print("????")
    f=open(file,"r")
    line=f.readline()
    while  "FORCE"  in line:
        f.readline()
        line=f.readline()
        #print(line)
    position=[]
    force=[]
    while "--" not in line:
        position+=[line.split()[0:2]]
        force+=[[float(line.split()[i])*constant.eV*1e10 for i in range(3,6)]]
        line=f.readline()
    return position,np.array(force)

def occupy(freq,temperature):
    a=constant.physical_constants["Planck constant over 2 pi"][0] * freq \
                        / temperature /constant.physical_constants["Boltzmann constant"][0]\
                            / constant.physical_constants["Planck constant over 2 pi in eV s"][0] #frequency换算
    return 1/(np.exp(a)-1)


def gaussian(omega, omega_k, sigma):#归一化高斯函数
    return 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(omega - omega_k)**2 / sigma**2 / 2)


def get_S_omega(omega, sigma,frequencies,S,jtmodes,skipmodes):
    sum = 0
    for k in range(skipmodes+1,len(S)):
        if k not in jtmodes:
            if frequencies[k]<0.005:
                sum += S[k] * gaussian(omega, frequencies[k], sigma)
            else:
                sum += S[k] * gaussian(omega, frequencies[k], 0.3*sigma)
    return sum

def get_C_omega(omega, sigma,frequencies,S,jtmodes,skipmodes,temperature):#omega实际是输入的能量 sigma为高斯函数展宽
    sum = 0
    for k in range(skipmodes+1,len(S)):
        if k not in jtmodes:
            sum += occupy(frequencies[k], temperature)*S[k] * gaussian(omega, frequencies[k], sigma)
    return sum

def el_ph_S(delta_width,omega_set,frequencies,S,jtmodes,skipmodes,index):
    #print("here is el_ph_S! len(omega_set)",len(omega_set))
    result=[get_S_omega(o, delta_width,frequencies,S,jtmodes,skipmodes) for o in omega_set]
    #print("task",index,"finished!")
    return result,index
    
def el_ph_C(delta_width,omega_set,frequencies,S,jtmodes,skipmodes,temperature,index):
    #print("here is el_ph_C! len(omega_set)",len(omega_set))
    result=[get_C_omega(o, delta_width,frequencies,S,jtmodes,skipmodes,temperature) for o in omega_set]
    #print("task",index,"finished!")
    return result,index

class lattice_inform:
    def __init__(self,lattice,natom,points):
        self.lattice_vector=lattice
        self.natom=natom
        self.atom=[]
        self.mass=[]
        self.atom_position=[]
        for a in points:
            self.atom.append(a['symbol'])
            self.mass.append(a['mass'])
            self.atom_position.append(a['coordinates'])
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
        
        
class phonon_mode:
    def __init__(self,data,mass):
        self.frequency=data['frequency']
        self.eigenvector=[]
        
        for i in range(len(data['eigenvector'])):
            aa=data['eigenvector'][i]
            self.eigenvector+=[[a[0]/np.sqrt(mass[i]) for a in aa]]
        
        self.eigenvector=np.array(self.eigenvector)/np.linalg.norm(np.array(self.eigenvector).reshape(-1,1))
        #print("normalized!!!")
        return

class phonon_band_in_one_kpoint:
    def __init__(self,position,band,mass):
        self.position=position
        self.phonon=[]
        self.mass=mass
        for p in band:
            self.phonon.append(phonon_mode(p,mass))
        


class Photoluminescence:
    def phonopy_read_modes(self):
        modes = np.zeros((self.numModes, self.numAtoms, 3))

        try:
            band = open(self.path, 'r')
        except OSError:
            print("Could not open/read file: band.yaml")
            sys.exit()

        for line in band:
            if "  band:" in line:
                break

        for i in range(self.numModes):
            band.readline()
            band.readline()
            band.readline()

            for a in range(self.numAtoms):
                line = band.readline()

                line = band.readline().replace(",", "")

                parts = line.strip().split()
                modes[i][a][0] = float(parts[2])

                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][1] = float(parts[2])

                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][2] = float(parts[2])

        band.close()

        return modes

    def phonopy_read_frequencies(self):
        frequencies = np.zeros(self.numModes)
        try:
            band = open(self.path, 'r')
        except OSError:
            print("Could not open/read file: band.yaml")
            sys.exit()

        for line in band:
            if "  band:" in line:
                break

        for i in range(self.numModes):
            band.readline()
            line = band.readline()

            parts = line.strip().split()
            frequencies[i] = float(parts[1])

            line = band.readline()

            for a in range(self.numAtoms):
                band.readline()
                band.readline()
                band.readline()
                band.readline()

        band.close()
        return frequencies

    def fold(self,d):
        nml= [ round(a) for a in np.linalg.solve(self.lattice_vector.T,d)]
        #if abs(sum(nml))>0.5:
            #print(nml)
        return d-np.matmul(self.lattice_vector.T,nml)
        
        
    def plotdata(self,data,label):
        plt.figure(figsize=(10, 10))
        plt.xlabel(label)
        plt.plot(data)
        plt.show()

    def get_phonon(self,path):
        f=open(path,"r")
        data=yaml.load(f,Loader=yaml.CLoader)
        lattice=lattice_inform(data['lattice'], data['natom'], data['points'])
        m=np.array([a*constant.physical_constants["atomic mass constant"][0] for a in lattice.mass])
        ph=phonon_band_in_one_kpoint(data['phonon'][0]['q-position'],data['phonon'][0]['band'],m) 
        modes=np.array([each.eigenvector for each in ph.phonon])
        fre=np.array([each.frequency for each in ph.phonon])
        
        f.close()
        return m,modes,fre
    
    def read_grd_ex_pos(self,str_g,str_e,shiftVector,n_defect,defect_range):
        if '.xyz' in str_g:
            self.g = XYZ(str_g).coordinates
            self.e = XYZ(str_e).coordinates
        else:
            cc = ConfigurationCoordinate()
            self.g = cc.read_poscar(str_g)
            self.e = cc.read_poscar(str_e)
            self.lattice_vector=self.g.lattice.matrix
            
            self.g.translate_sites(
                range(len(self.g.frac_coords)), shiftVector, frac_coords=False)
            self.e.translate_sites(
                range(len(self.e.frac_coords)), shiftVector, frac_coords=False)
            
            
            lg = self.g.lattice
            le = self.e.lattice
            self.g = lg.get_cartesian_coords(self.g.frac_coords)
            self.e = le.get_cartesian_coords(self.e.frac_coords)
        
    
        self.numAtoms = len(self.g)
        self.D_R=[]
        for i in range(self.numAtoms):
            self.D_R.append(self.fold(self.e[i]-self.g[i]).tolist())
            #if abs(self.D_R[i][0])>0.1 or abs(self.D_R[i][1])>0.1 or abs(self.D_R[i][2])>0.1:
                #print(i+1,self.D_R[i])
        self.D_R=np.array(self.D_R)
        
        #solve the period lattice problem
    
        if n_defect>0:
            n_defect=n_defect-1
            sumd=np.array([0.0,0.0,0.0])
            count=0
            for i in range(self.numAtoms):
                if np.dot(self.e[i]-self.e[n_defect],self.e[i]-self.e[n_defect]) > defect_range :
                    sumd+=self.D_R[i]
                    count+=1
            sumd/=float(count)
            for i in range(self.numAtoms):
                self.D_R[i]-=sumd

    def get_S_omega(self, omega, sigma):
        sum = 0
        for k in range(self.skipmodes+1,len(self.S)):
            if k not in self.jtmodes:
                if self.frequencies[k]<0.035:
                    sum += self.S[k] * gaussian(omega, self.frequencies[k], sigma)
                else:
                    sum += self.S[k] * gaussian(omega, self.frequencies[k], sigma)
        return sum
    
       
    def get_C_omega(self, omega, sigma):#omega实际是输入的能量 sigma为高斯函数展宽
        sum = 0
        for k in range(self.skipmodes+1,len(self.S)):
            if k not in self.jtmodes:
                sum += self.occupy(self.frequencies[k])*self.S[k] * gaussian(omega, self.frequencies[k], sigma)
        return sum
        
    def occupy(self,freq):
        a=constant.physical_constants["Planck constant over 2 pi"][0] * freq \
                         / self.temperature /constant.physical_constants["Boltzmann constant"][0]\
                             / constant.physical_constants["Planck constant over 2 pi in eV s"][0] #frequency换算
        return 1/(np.exp(a)-1)

    def write_S(self, file_name):
        f = open(file_name, 'w')
        for i in range(len(self.S_omega)):
            # f.write(str(self.omega_set[i]) + "\t" + str(self.S_omega[i])+'\n')
            f.write(str(self.S_omega[i])+'\n')
        f.close()


    def el_ph(self,**parameter):
        self.temperature = parameter.get("temperature",0)
        self.delta_width=parameter.get("delta_width",6e-3)
        self.jtmodes=parameter.get("jtmodes",[])
        self.omega_set = np.linspace(
            self.min_energy, self.max_energy, int((self.max_energy-self.min_energy)*self.resolution)) #in energy(maybe eV)
        self.S_omega = [self.get_S_omega(o, self.delta_width) for o in self.omega_set]
        if self.temperature>1e-2:
            self.C_omega = [self.get_C_omega(o, self.delta_width) for o in self.omega_set]
        else:
            self.C_omega=self.S_omega
        with open("C_omega.data","w") as f:
            for a in self.C_omega:
                f.write(str(a)+"\n")
        print("el-ph OK!")


    def el_ph_(self,**parameter):
        self.temperature = parameter.get("temperature",0)
        self.delta_width=parameter.get("delta_width",6e-3)
        self.jtmodes=parameter.get("jtmodes",[])
        self.omega_set = np.linspace(self.min_energy, self.max_energy, int((self.max_energy-self.min_energy)*self.resolution)) #in energy(maybe eV)
        

        num_thread=8
        lenth=math.ceil(len(self.omega_set)/num_thread)
        #print("lenth",len(self.omega_set)/num_thread,lenth)
        pool=mp.Pool(processes=num_thread)
        result=[]
        
        for i in range(num_thread):
            pool.apply_async( el_ph_S ,[ self.delta_width,self.omega_set[i*lenth:(i+1)*lenth] ,self.frequencies,self.S,self.jtmodes,self.skipmodes,i],callback=result.append)
        pool.close()
        pool.join()
        self.S_omega=[]
        for j in range(num_thread):
            for data in result:
                if j==data[1]:
                    self.S_omega+=data[0]

        if self.temperature>1e-2:
            result=[]
            for i in range(num_thread):
                pool.apply_async( el_ph_C ,[ self.delta_width,self.omega_set[i*lenth:(i+1)*lenth] ,self.frequencies,self.S,self.jtmodes,self.skipmodes,self.temperature,i],callback=result.append)
            pool.close()
            pool.join()
            self.C_omega=[]
            for j in range(num_thread):
                for data in result:
                    if j==data[1]:
                        self.C_omega+=data[0]
        else:
            self.C_omega=self.S_omega

        print("el-ph OK!")
        return
        
    
    def PL(self, **parameter): 
        self.EZPL=parameter.get("EZPL",0)
        self.gamma=parameter.get("gamma",1)
        self.process=parameter.get("process","emission")
        SHR=parameter.get("SHR",0)

        Gt = []
        
        r = 1/self.resolution
        if "emission" in self.process:
            St = realfft(self.S_omega,self.max_energy-self.min_energy)
            Ct = realfft(self.C_omega,self.max_energy-self.min_energy)
            St = fft.fftshift(St)
            Ct = fft.fftshift(Ct)
        elif "absorption" in self.process:
            St = realifft(self.S_omega,self.max_energy-self.min_energy)
            Ct = realifft(self.C_omega,self.max_energy-self.min_energy)
            St = fft.fftshift(St)
            Ct = fft.fftshift(Ct)
        else:
            print("process input error, input:",self.process)
            sys.exit()
        #self.plotdata(St,"st")
        
        if self.temperature>1e-2 :
            for i in range(len(St)):
                t = r*(i-len(St)/2)
                Gt += [np.exp(St[i]+Ct[i]+Ct[len(St)-i-1]-St[0]-2*Ct[0]-SHR-self.gamma*np.abs(t))]
        elif self.temperature<=1e-2 :
            for i in range(len(St)):
                t = r*(i-len(St)/2)
                Gt += [np.exp(St[i]-St[0]-SHR-self.gamma*np.abs(t))]
        #self.plotdata(Gt,"Gt")
        A = fft.ifft(Gt)*len(Gt)/(self.max_energy-self.min_energy)
        tA = A.copy()
        for i in range(len(A)):
            A[(int((self.EZPL-self.min_energy)*self.resolution)-i) % len(A)] = tA[i]
            
        self.A=A
        return A
        #print(A.tolist())

    def PLA(self):
        # Now, shift the ZPL peak to the EZPL energy value

        I = []
        r = 1/self.resolution
        if "emission" in self.process:
            for i in range(len(self.A)):
                #t = r*(i-len(self.A)/2)
                t = r*(i+self.min_energy*self.resolution)
                I += [self.A[i]*((t)*r)**3]
        elif "absorption" in self.process:
            for i in range(len(self.A)):
                #t = r*(i-len(self.A)/2)
                t = r*(i+self.min_energy*self.resolution)
                I += [self.A[i]*((t)*r)]
        else:
            print("process input error, input:",self.process)
            sys.exit()
            
        self.I=np.array(I)
        
        return np.array(I)

    def __init__(self, path, method ,**parameter ):
        self.path = path 
        self.method = method #phonopy or phonopy-siesta or vasp used for phonon calculation
        self.resolution = parameter.get("resolution",1000) #resolution means points in 1eV

        if "large" in method :
            str_g= parameter.get("POSCAR_GRD")
            pos,self.fg=read_force(parameter.get("force_grd"))
            #print(self.fg)
            pos,self.fe=read_force(parameter.get("force_ex"))
            shiftVector=parameter.get('shift_vector=', [0.0,0.0,0.00])
            if '.xyz' in str_g:
                self.g = XYZ(str_g).coordinates
            else:
                cc = ConfigurationCoordinate()
                self.g = cc.read_poscar(str_g)
                self.lattice_vector=self.g.lattice.matrix
                
                self.g.translate_sites(
                    range(len(self.g.frac_coords)), shiftVector, frac_coords=False)
                
                lg = self.g.lattice
                self.g = lg.get_cartesian_coords(self.g.frac_coords)
    
            self.numAtoms = len(self.g)
            
        else:
             self.read_grd_ex_pos( parameter.get("POSCAR_GRD"),parameter.get("POSCAR_EX"),parameter.get('shift_vector=', [0.0,0.0,0.00]),parameter.get("n_defect",0),parameter.get("defect_range",4.0))   
            
        #struct aligment
        
        #structure data in Angstrom (1e-10m)    
     

        if parameter.get('m',0) !=0: 
            self.numModes = parameter.get('nmodes', 3*self.numAtoms)
            self.Modes=self.phonopy_read_modes()
            self.frequencies=self.phonopy_read_frequencies()
            self.m=parameter.get('m')
        else:
            self.m,self.Modes,self.frequencies=self.get_phonon(path)
            self.numModes=len(self.Modes)#kg
            
        for i in range(self.numModes):
            if "vasp" in method:
                self.frequencies[i] = self.frequencies[i] / 1000
            elif "phonopy-siesta" in method:
                self.frequencies[i] = self.frequencies[i] * \
                    0.004135665538536 * 0.727445665  # THz
            elif  "phonopy" in method:
                self.frequencies[i] = self.frequencies[i] * \
                    0.004135665538536  # THz
            if self.frequencies[i] <= 0.005:
                self.skipmodes=i
        print("!!!!!!!!!!!!!Warning!!!!!!!!!!!!\n skipmode energy 0.005")    
        print("skipmodes",self.skipmodes)
        #print("Phonon modes readed successfully!")
        #read phonon modes
         #Photoluminescence can accept use's input of mass list (m=[m1,m2....] and  nomber of  modes (nmodes = int),then it reads eigenvector and frequences automatically.
         #or if Photoluminescence  accepts none of these imformation, it will read all these needed data automatically but inefficiently.
         
         
    def HuangRhyes(self):

        self.HuangRhyes = 0
        self.Delta_R = 0
        self.Delta_Q = 0
        self.IPR = []
        self.q = []
        self.q_amu = []
        self.S = []
        
        atomunit=np.sqrt(constant.physical_constants["atomic mass constant"][0])/1e10
        
        for i in range(self.numModes):
            
            q_i = 0
            IPR_i = 0
            participation = 0
            #频率为振子能量（eV）           
            max_Delta_r = 0
            
            if "large" in self.method:
                for a in range(self.numAtoms):
                    # Normalize r:
                    #participation = self.Modes[i][a][0] * self.Modes[i][a][0] + \
                        #self.Modes[i][a][1] * self.Modes[i][a][1] + self.Modes[i][a][2] * self.Modes[i][a][2]
                    participation=np.dot(self.Modes[i][a],self.Modes[i][a])
                    IPR_i += participation**2#所有位移的平方和
                    if i > 2:
                        q_i +=  np.dot(self.fe[a]-self.fg[a], self.Modes[i][a])/np.sqrt(self.m[a])\
                        /(self.frequencies[i]/constant.physical_constants["Planck constant over 2 pi in eV s"][0])**2
                    else:
                        q_i =0
                
            else:
                    
                for a in range(self.numAtoms):
                    # Normalize r:
                    #participation = self.Modes[i][a][0] * self.Modes[i][a][0] + \
                        #self.Modes[i][a][1] * self.Modes[i][a][1] + self.Modes[i][a][2] * self.Modes[i][a][2]
                    participation=np.dot(self.Modes[i][a],self.Modes[i][a])
                    IPR_i += participation**2#所有位移的平方和
                    q_i += np.sqrt(self.m[a]) * np.dot(self.D_R[a],self.Modes[i][a]) * 1e-10 #Angstrom
                    if i==0 or np.dot(self.D_R[a],self.D_R[a])>0.1:
                        print(a,np.dot(self.D_R[a],self.D_R[a]))



            IPR_i = 1.0 / IPR_i#inverse participation ratio
            S_i = self.frequencies[i] * q_i**2 / 2 * 1.0 / \
                (constant.physical_constants["Planck constant over 2 pi"][0] * constant.physical_constants["Planck constant over 2 pi in eV s"][0])

            #print(S_i)
            self.IPR += [IPR_i]
            self.q += [q_i]#位形坐标
            self.q_amu +=[q_i/atomunit]
            self.S += [S_i]#部分黄昆因子
            self.HuangRhyes += S_i#黄昆因子 
            
        if "large" not in self.method:
            f=open("D(e-g).data","w")
            for a in range(self.numAtoms):
                self.Delta_R += sum((self.D_R[a])**2)
                self.Delta_Q += sum((self.D_R[a])**2 * self.m[a])
                f.write(str(self.D_R[a])+"\t")
                f.write("\n")
            f.close()
        
        
        self.Delta_R = (self.Delta_R)**0.5

        self.Delta_Q = (self.Delta_Q / constant.physical_constants["atomic mass constant"][0]) ** 0.5

        self.max_energy = 5.0
        self.min_energy = -0.0
        #print("S computed successfully!")
        return self.HuangRhyes
    


    def print_table(self,k,visualmode):
        atomunit=np.sqrt(constant.physical_constants["atomic mass constant"][0])/1e10
        f=open("main_modes.data","w")
        f2=open("modes.data","w")
        f.write("Modes NO.\t DeltaQ(A^2)\t IPR\t energy(meV)\t partial HuangRhyes\n ")
        f2.write("Modes NO.\t DeltaQ(A^2)\t IPR\t energy(meV)\t partial HuangRhyes\n ")
        for i in range(self.numModes):
            
            f2.write(" "+str(i+1)+"\t"+str(self.q[i]/atomunit)+"\t"+str(self.IPR[i])+" \t"+str(self.frequencies[i]*1000)+" \t"+str(self.S[i])+"\n")
            if(self.S[i]>k and self.frequencies[i]>0.04 or i in visualmode):
                print("mode No.",i+1," q=",self.q[i]/atomunit," IPR=",round(self.IPR[i],2)," energy=",round(self.frequencies[i]*1000,3),"meV  Sk=",self.S[i])
                f.write(" "+str(i+1)+"\t"+str(self.q[i]/atomunit)+"\t"+str(self.IPR[i])+" \t"+str(self.frequencies[i])+" \t"+str(self.S[i])+"\n")
                #print("IPR\t", i, "\tSk\t", self.S[i], "\tenergy\t",\
                #      self.frequencies[i], "\t=\t", self.IPR[i], "\twith localization ratio beta =\t", 64 / self.IPR[i])
        f.close()   
        f2.close()
