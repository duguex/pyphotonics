# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 11:52:39 2023

@author: 87930
"""

import numpy as np
import matplotlib.pyplot as plt
import math

class phonon:
    def __init__(self,index,energy,S,TDM=np.zeros(3,dtype=float)):
        self.index=index-1
        self.energy=energy
        self.S=S
        self.TDM=TDM
        
        self.lenth=math.ceil(4*S)+5
        self.solve()
        return
    
    def solve(self):
        lenth=math.ceil(5*self.S)+10
        mat=np.zeros(lenth**2,dtype=float).reshape(lenth,lenth)
        g=np.sqrt(self.S)
        for i in range(lenth-1):
            mat[i][i]=i
            mat[i+1][i]=g*np.sqrt(i+1)
            mat[i][i+1]=g*np.sqrt(i+1)
        mat[lenth-1][lenth-1] = lenth-1 
        #print(mat)
            
        E,a=np.linalg.eig(mat.tolist())
        a=a.T
        sort=np.argsort(E)
        E=np.array([E[i] for i in sort])
        a=np.array([a[i] for i in sort])   
        for i in range(len(a)):
            if a[i][0]<0:
                a[i]=-a[i]
        a=a.T
        self.FC=a[0][0:math.ceil(2*self.S)+5]#<0|a>=a[0]
        self.HT=a[1][0:math.ceil(2*self.S)+5]#<0|x|a>=<1|a>
        #plt.stem(self.FC,markerfmt="ro",label="FC")
        #plt.stem(self.HT,markerfmt="bo",label="HT")
        #plt.legend()
        #plt.show()
        self.FCenergy=self.energy*np.linspace(0,math.ceil(2*self.S)+4,math.ceil(2*self.S)+5,dtype=float)
        return

    def spectral(self):
        a=0.01
        k=1
        index=1
        gy=[]
        while k>0.001:
            k=np.exp(-a*index**2)
            gy.append(k)
            index+=1
            

        X=self.FCenergy
        Y=self.HT
           

        res=1
        maxe=4000
        offset=500
        inten=np.array(np.zeros(int(res*maxe),dtype=float))
        #inteng=np.array(np.zeros(int(res*5000),dtype=float))
        x=np.array(np.zeros(int(res*maxe),dtype=float))
        for i in range(len(x)):
            x[i]=-i+3500
        
        for i in range(len(X)):
            loca=round(res*X[i])+offset
            inten[loca]+=Y[i]
            
            for j in range(len(gy)):
                inten[loca+j+1]+=gy[j]*Y[i]
                inten[loca-j-1]+=gy[j]*Y[i]
                
        self.HTspectral=inten[::-1]
        #plt.plot(x,inten[::-1])
        #plt.show()
        return


class band:
    def __init__(self,TDM=np.array([0.0,0.0,0.0])):
        self.TDM=TDM
        self.ph={}
        self.count=0
        self.spectral()
        return
    
        
    def addphonon(self,index,energy,S,TDM=np.zeros(3,dtype=float)):
        self.ph[index-1]=phonon(index,energy,S,TDM)
        self.count+=1
        return
    
    def energy_add(self,a,b):
        return (np.tile(a,(len(b),1)).T+np.tile(b,(len(a),1))).reshape(1,-1)[0]
    
    def FC_luminescece(self):
        lumine=np.array([1.0])
        ener=np.array([0.0])
        for key in self.ph:
            #print(key,"fc",self.ph[key].FC)
            #print(lumine)
            lumine=np.kron(lumine,self.ph[key].FC)
            ener=self.energy_add(ener,self.ph[key].FCenergy)
            
        return ener,lumine**2
    
    def spectral(self):
        a=0.01
        k=1
        index=1
        gy=[]
        while k>0.001:
            k=np.exp(-a*index**2)
            gy.append(k)
            index+=1
            

        X=[0.0]
        Y=[1.0]
           
        res=1
        maxe=4000
        offset=500
        inten=np.array(np.zeros(int(res*maxe),dtype=float))
        #inteng=np.array(np.zeros(int(res*5000),dtype=float))
        
        for i in range(len(X)):
            loca=round(res*X[i])+offset
            inten[loca]+=Y[i]
            
            for j in range(len(gy)):
                inten[loca+j+1]+=gy[j]*Y[i]
                inten[loca-j-1]+=gy[j]*Y[i]
                
        self.HTspectral=inten[::-1]
        #plt.plot(inten[::-1])
        #plt.show()
        return
    
    
    def FC_HT_luminescece(self):
        lumine=np.tile(np.array([1.0]),(self.count+1,1)).tolist()
        print(lumine)
        ener=np.array([0.0])
        i=0
        tdm=[self.TDM]
        for key in self.ph:
            i+=1
            #print(i,self.count)
            #print(key,"fc",self.ph[key].FC)
            
            tdm.append(self.ph[key].TDM)
            
            lumine[i]=np.kron(lumine[0],self.ph[key].HT)
            lumine[0]=np.kron(lumine[0],self.ph[key].FC)
            for j in range(i-1):
                lumine[j+1]=np.kron(lumine[j+1],self.ph[key].FC)
            ener=self.energy_add(ener,self.ph[key].FCenergy)
            print(lumine)
            #print("111",len(ener),len(lumine[0]))
        
            
        tdm=np.array(tdm)
        lumine=np.array(lumine).T
        #print(lumine)
        lumineall=np.zeros(len(lumine),dtype=float)
        for i in range(len(lumineall)):
            #print("test",lumine[i],tdm,np.matmul(lumine[i],tdm))
            lumineall[i]=np.linalg.norm(np.matmul(lumine[i],tdm))**2
            #print("LUMINE[I]",lumineall[i])
        return ener,lumineall
 
    
def htt(file):

        
    ba=band(TDM=np.array([-0.0111,-0.01923,0.04134]))
    #ba.addphonon(1, 3, 0)
    #ba.addphonon(120, 296, 1e-16,TDM=np.array([-0.05406,-5e-7,-2e-7]))
    #b.addphonon(118, 296., 1e-16,TDM=np.array([-1.1e-6,0.05411,1e-7]))
    f=open(file,"r")
    line=f.readline()
    while line:
        s=line.split()
        ba.addphonon(int(s[0]),float(s[1]),float(s[2]),TDM=np.array([float(ss) for ss in s[3:]]))
        line=f.readline()
        print(int(s[0]),float(s[1]),float(s[2]),np.array([float(ss) for ss in s[3:]]))
    return ba


    



def odd(file):   
    a=0.02
    k=1
    index=1
    gy=[]
    while k>0.01:
        k=np.exp(-a*index**2)
        gy.append(k)
        index+=1
        
    print("gy",gy,len(gy))
        
      
    
    ba=band(TDM=np.array([0.0,0.0,0.01]))
    #ba.addphonon(1, 3, 0)
    #ba.addphonon(120, 296, 1e-16,TDM=np.array([-0.05406,-5e-7,-2e-7]))
    #b.addphonon(118, 296., 1e-16,TDM=np.array([-1.1e-6,0.05411,1e-7]))
    f=open(file,"r")
    line=f.readline()
    while line:
        s=line.split()
        ba.addphonon(int(s[0]),float(s[1]),float(s[2]),TDM=np.array([float(ss) for ss in s[3:]]))
        line=f.readline()
        print(int(s[0]),float(s[1]),float(s[2]),np.array([float(ss) for ss in s[3:]]))
    
    X,Y=ba.FC_HT_luminescece()
    
    print(len(X))
    print(sum(Y))
    
    
    
    zpl=13796+480
    res=1
    maxe=10000
    offset=1000
    inten=np.array(np.zeros(int(res*maxe),dtype=float))
    #inteng=np.array(np.zeros(int(res*5000),dtype=float))
    x=[-i/res+zpl+offset/res for i in range(int(res*maxe))]
    
    for i in range(len(X)):
        loca=round(res*X[i])+offset
        inten[loca]+=Y[i]
        
        for j in range(len(gy)):
            inten[loca+j+1]+=gy[j]*Y[i]
            inten[loca-j-1]+=gy[j]*Y[i]
           
    
    #for i in range(len(gy)):
        #inteng=gy[i]*inten[-i-1:len(inten)-i-1]
        #inteng=gy[i]*inten[i+1-len(inten):i+1]
    #print(len(inteng),len(inten))   
    plt.figure(figsize=(7, 3.5),dpi=300)
    plt.plot(x,inten)
    plt.xlim(-1000+zpl,zpl+100)
    plt.title(file)
    plt.vlines(zpl, 0, max(inten),color="y")
    plt.vlines(zpl-480, 0, max(inten),color="r")
    plt.show()    
    return inten[::-1]
        
#odd()      
        
if __name__ == '__main__':
    ba=band(TDM=np.array([1.0,0.0,0.0]))
    ba.addphonon(1, 1, 0.15,TDM=np.array([0.3,0.0,0.0]) )
    #ba.addphonon(2, 2.5, 0.2,TDM=np.array([0.0,0.0,0.0]) )   
    x,y=ba.FC_HT_luminescece()  
            



    plt.stem(x,y)
    plt.show()
            
    
    
        