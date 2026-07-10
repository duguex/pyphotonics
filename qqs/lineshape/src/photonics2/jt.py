# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 10:21:22 2022

@author: 87930
"""
#demensionless jt solver
from unittest import result
from joblib import PrintTime
import numpy as np
import matplotlib.pyplot as plt
import math 
import multiprocessing as mp
import copy 
import concurrent.futures as concurrent
import multiprocessing as mp
from photonics.ht import band,phonon
from sympy import false, minimal_polynomial
import photonics.hermite as hermite
from matplotlib.pyplot import MultipleLocator, tight_layout

def nl(k):
    #print("nl_input_k",k)
    n=round(math.ceil(np.sqrt(2*k+2.25)-1.5))
    l=round(2*k-(n+2)*n)
    return n,l

def kk(n,l):
    return round(((n+2)*n+l)/2)

def jtSolverT(kt,na,theta):#chiral phonon base
    A=np.sqrt(0.5)*np.exp(1j*theta)
    mat0=np.zeros(((na+1)*(na+2))**2,dtype=np.complex64).reshape((na+1)*(na+2),(na+1)*(na+2))
    for k in range((na+1)*(na+2)):
        mat0[k][k]=nl(k)[0]+1
        n,l=nl(k)
        if n > 0:
            if l+1 <= n-1:
                mat0[k][ kk( n-1 , l+1 ) ]=A*kt*np.sqrt((n-l)/2.0)
            if l-1 >= 1-n:
                mat0[k][ kk( n-1 , l-1 ) ]=A.conjugate()*kt*np.sqrt((n+l)/2.0)
               
        if n <= na:
            mat0[k][ kk( n+1 , l+1 ) ]=A*kt*np.sqrt((n+l+1)/2.0)
            mat0[k][ kk( n+1 , l-1 ) ]=A.conjugate()*kt*np.sqrt((n-l+2)/2.0)

    E,a=np.linalg.eig(mat0)
    a=a.T
    #print(E.tolist())
    E=np.array([e.real for e in E])
    
    sort=np.argsort(E)
    E=np.array([E[i] for i in sort])
    a=np.array([a[i] for i in sort])
    """
    for i in range(len(E)):
        print("E(p=",i,")=",E[i])
        print(a[i].__abs__())
    """
    return E,a


def jtSolverTOringinalBase_New(kt,na,**para):
    print("================T*e solver module start!======================")

    ratio=para.get("ratio",[1.0,1.0])
    gex=(ratio[0]**2-1)/2
    gey=(ratio[1]**2-1)/2    
    
    hwx=para.get("hwx",1.0)
    hwy=para.get("hwy",1.0)

    dim=round((na+1)*(na+2)/2)
    mat0=np.zeros(dim**2,dtype=np.float64).reshape(dim,dim)

    
    

    print("hwx hwy",hwx,hwy)
    print("ratio",ratio)
    for key in kt:
        print(key,kt[key].k)
    
    for k in range(dim):
        n,x,y=nxy(k,na)

        mat0[k][k]=(x+0.5)*hwx+(y+0.5)*hwy
        write_mat(mat0, 2, 0, x, y, x, y, na, gex, n, n)
        write_mat(mat0, 0, 2, x, y, x, y, na, gey, n, n)
        for key in kt:
            write_mat(mat0,kt[key].nx,kt[key].ny,x,y,x,y,na,kt[key].k,n,n)



    E,a=np.linalg.eig(mat0)
    a=a.T
    E=np.array([e.real for e in E])
    sort=np.argsort(E)
    E=np.array([E[i] for i in sort])
    #print(E.tolist())
    a=np.array([a[i] for i in sort])
    print("================T*e solver module End!======================")
    return E,a





def nxy(k,na):
    n= int(np.floor(k/(na+2)/(na+1)*2))
    #print(n)
    ne,l=nl(k-round(n*(na+2)*(na+1)/2))
    x=round((ne+l)/2)
    y=round((ne-l)/2)
    return n,x,y

def kkk(n,x,y,na):
    return int(round((n*round((na+2)*(na+1)/2)+kk(x+y,x-y))))

#write matrix H_jt
def write_mat(mat,n1,n2,x0,y0,x,y,na,k,pm1,pm2): #k_x^n1_y^n2|x,y,pm1>
    if n1>0:
        write_mat(mat,n1-1,n2,x0,y0,x+1,y,na,k*np.sqrt((x+1)/2),pm1,pm2)
        if x>0:
            write_mat(mat,n1-1,n2,x0,y0,x-1,y,na,k*np.sqrt(x/2),pm1,pm2) 
    elif n2>0:
        write_mat(mat,n1,n2-1,x0,y0,x,y+1,na,k*np.sqrt((y+1)/2),pm1,pm2)
        if y>0:
            write_mat(mat,n1,n2-1,x0,y0,x,y-1,na,k*np.sqrt(y/2),pm1,pm2)
    else:
        if y+x>na:
            return
        else:
            mat[kkk(pm1,x0,y0,na)][kkk(pm2,x,y,na)]+=k
    return



def jtSolver_Split_OringinalBase_New(ke,na,**para):#in original base ke=[xx,yy,xy]
    dim=(na+2)*(na+1)#na: max number of total phonons 
    mat0=np.zeros(dim**2,dtype=np.float64).reshape(dim,dim)
    Sx=para.get("Sx",0)
    Sy=para.get("Sy",0)
    hwx=para.get("hwx",1)
    hwy=para.get("hwy",1)
    ge=para.get("ge",[0,0,0])
    v=para.get("v",[0,0])
    
    print("================E*e solver module start!======================")
    #print("w=",ke[0],"v=",ke[1],"diag",v)
    for key in ke[0]:
        print("w",key,ke[0][key].k)
    for key in ke[1]:
        print("v",key,ke[1][key].k)
    for key in v:
        print("diag",key,v[key].k)    
    print("Sx",Sx,"Sy",Sy)
    
    for k in range(dim):
        #print("k,na",k,na,(na+2)*(na+1))
        n,x,y=nxy(k,na)
        #print("knxy",k,n,x,y)
        mat0[k][k]=Sx*(1-2*n)+x*hwx+y*hwy
        mat0[k][kkk(1-n,x,y,na)]=-Sy
        #third order 
        for key in v:
            write_mat(mat0,v[key].nx,v[key].ny,x,y,x,y,na,v[key].k,n,n)
        for key in ke[0]:
            if n==0:
                write_mat(mat0,ke[0][key].nx,ke[0][key].ny,x,y,x,y,na,ke[0][key].k,n,n)
            if n==1:   
                 write_mat(mat0,ke[0][key].nx,ke[0][key].ny,x,y,x,y,na,-ke[0][key].k,n,n)
        for key in ke[1]:    
            write_mat(mat0,ke[1][key].nx,ke[1][key].ny,x,y,x,y,na,ke[1][key].k,n,1-n)
        
    E,a=np.linalg.eig(mat0.tolist())
    
    a=a.T
    E=np.array([e.real for e in E])
    sort=np.argsort(E)
    E=np.array([E[i] for i in sort])
    #print(E.tolist())
    a=np.array([a[i] for i in sort])
    print("=================E*e solver module End!==================")
    return E,a

def jtSolver(k,na,m):#chiral phonon base
    na=na+1
    mat0=np.zeros(na*na,dtype=np.float64).reshape(na,na)
    for i in range(na):
        mat0[i][i]=abs(m)+1+i
    for i in range(int(na/2)):
        mat0[2*i][2*i+1]=k*np.sqrt(i+1+(abs(m)+m)/2)
        mat0[2*i+1][2*i]=k*np.sqrt(i+1+(abs(m)+m)/2)
    for i in range(int(na/2)-1):
        mat0[2*i+2][2*i+1]=k*np.sqrt(i+1+(abs(m)-m)/2)
        mat0[2*i+1][2*i+2]=k*np.sqrt(i+1+(abs(m)-m)/2)
    E,a=np.linalg.eig(mat0)
    a=a.T
    E=np.array([e.real for e in E])
    sort=np.argsort(E)
    E=np.array([E[i] for i in sort])
    a=np.array([a[i] for i in sort])
    """
    for i in range(len(E)):
        print("E(m=",m,",p=",i,")=",E[i])
        print(a[i].__abs__())
    """
    return E,a

def convertToFullBase(a,m,na):
    full=np.zeros((na+1)*(na+2),dtype=float)
    for i in range(len(a)):
        if kk( m+i , int(round(m-0.5*(-1)**i+0.5))) <(na+1)*(na+2):
            if(kk( m+i , int(round(m-0.5*(-1)**i+0.5)))<=(na+1)*(na+2)):
                full[ kk( m+i , int(round(m-0.5*(-1)**i+0.5))) ]=a[i].real

    return full
        

def TtoE(Ee,e_state,Et,t_state,na,me):
    intensity=[]
    E=[]
    for i in range(len(Ee)):
        intensity+=[abs(np.dot(t_state,convertToFullBase(e_state[i],me,na)))**2]
        E+=[Et-Ee[i]]
    return  E,intensity



##ax·X+ay·Y|phi>
def XdotPhi(phi,ax,ay,na):
    xphi=np.array(np.zeros(len(phi),dtype=float))#创建0向量
    phi2=np.array(copy.deepcopy(phi))
    for k in range(len(phi2)):
        n,x,y=nxy(k,na)
        if x>0:
            xphi[kkk(n,x-1,y,na)]+=ax/np.sqrt(2*x)*phi2[k]#降算符
        if y>0:
            xphi[kkk(n,x,y-1,na)]+=ay/np.sqrt(2*y)*phi2[k]#降算符
        if x+y<na:
            xphi[kkk(n,x+1,y,na)]+=ax/np.sqrt(2*(x+1))*phi2[k]#升算符
            xphi[kkk(n,x,y+1,na)]+=ay/np.sqrt(2*(y+1))*phi2[k]#升算符
    return xphi
        
def XpolydotPhi(phi,nx,ny,na):
    if nx>0:
        return XdotPhi(XpolydotPhi(phi,nx-1,ny,na),1,0,na)
    elif ny>0:
        return XdotPhi(XpolydotPhi(phi,nx,ny-1,na),0,1,na)
    else:    
        return phi 

def single_TtosplitE(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,index):
    #print("here is process!")
    intensity_FC_Ex=[]
    intensity_HT_Ex=[]
    intensity_FC_Ey=[]
    intensity_HT_Ey=[]
    ax=mixa[0]
    ay=mixa[1]#Q=ax·X+ay·Y
    E=[]
    mid=round(len(e_state[0])/2)
    for i in range(len(Ee)):#(len(Ee)):
        ##### FC term #####
        x2=ax*e_state[i][0:mid]+ay*e_state[i][mid:]
        y2z2=-ay*e_state[i][0:mid]+ax*e_state[i][mid:]
        intenEx=np.dot(t_state,x2)
        intenEy=np.dot(t_state,y2z2)
        #### HT term ####
        inten2Ex=copy.deepcopy(miu1[0])*0.0
        inten2Ey=copy.deepcopy(miu1[0])*0.0
        for k in range(len(miu2)-1-4):
            #print(nl(k+1))
            n,l=nl(k+1)
            nx=round((n-l)/2)
            ny=round((n+l)/2)
            dQ_to_dx=0.319876**n
            x_state=XpolydotPhi(e_state[i],nx,ny,na)
            x_x2=ax*x_state[0:mid]+ay*x_state[mid:]
            x_y2z2=-ay*x_state[0:mid]+ax*x_state[mid:]
            inten2Ex+=dQ_to_dx*miu1[k+1]*np.dot(t_state,x_x2)
            inten2Ey+=dQ_to_dx*miu2[k+1]*np.dot(t_state,x_y2z2)

        intensity_FC_Ex+=[intenEx]
        intensity_HT_Ex+=[inten2Ex]
        intensity_FC_Ey+=[intenEy]
        intensity_HT_Ey+=[inten2Ey]

        E+=[Et-Ee[i]]
    return index,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey



def TtosplitE_mp(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,polarization=np.array([0.,0.,0.]),maxpp=100):
    intensity_FC_Ex=[]
    intensity_HT_Ex=[]
    intensity_FC_Ey=[]
    intensity_HT_Ey=[]
    E=[]
    #print("direction ",polarization)
    if np.linalg.norm(polarization)>1e-10:
        polarization/=np.linalg.norm(polarization)
        for i in range(len(miu1)):
            print(miu1[i],polarization,np.dot(miu1[i],polarization))
            miu1[i]=np.dot(miu1[i],polarization)*polarization
            miu2[i]=np.dot(miu2[i],polarization)*polarization
            
    num_thread=32
    lenth=int(round(maxpp/num_thread))
    print("lenth",lenth)
    print("is that right?")
    with mp.Pool(processes=num_thread) as pool:
        result=[]
        for i in range(num_thread):
            pool.apply_async(single_TtosplitE, [Ee[lenth*i:lenth*(i+1)],e_state[lenth*i:lenth*(i+1)],Et,t_state,na,miu1,miu2,mixa,i],callback=result.append)

        pool.close()
        pool.join()

    for i in range(num_thread):
        for data in result:
            if data[0]==i:
                #print("data",data)
                E+=data[1]
                intensity_FC_Ex+=data[2]
                intensity_HT_Ex+=data[3]
                intensity_FC_Ey+=data[4]
                intensity_HT_Ey+=data[5]




    return np.array(E),np.array(intensity_FC_Ex),np.array(intensity_HT_Ex),np.array(intensity_FC_Ey),np.array(intensity_HT_Ey)



def single_splitEtoT(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,index):
    #print("here is process!")
    intensity_FC_Ex=[]
    intensity_HT_Ex=[]
    intensity_FC_Ey=[]
    intensity_HT_Ey=[]
    ax=mixa[0]
    ay=mixa[1]#Q=ax·X+ay·Y
    E=[]
    mid=round(len(e_state)/2)
    for i in range(len(Et)):#(len(Ee)):
        ##### FC term #####
        x2=ax*e_state[0:mid]+ay*e_state[mid:]
        y2z2=-ay*e_state[0:mid]+ax*e_state[mid:]
        intenEx=np.dot(t_state[i],x2)
        intenEy=np.dot(t_state[i],y2z2)
        #### HT term ####
        inten2Ex=copy.deepcopy(miu1[0])*0.0
        inten2Ey=copy.deepcopy(miu1[0])*0.0
        for k in range(len(miu2)-1-4):
            #print(nl(k+1))
            n,l=nl(k+1)
            nx=round((n-l)/2)
            ny=round((n+l)/2)
            dQ_to_dx=0.319876**n
            x_state=XpolydotPhi(e_state,nx,ny,na)
            x_x2=ax*x_state[0:mid]+ay*x_state[mid:]
            x_y2z2=-ay*x_state[0:mid]+ax*x_state[mid:]
            inten2Ex+=dQ_to_dx*miu1[k+1]*np.dot(t_state[i],x_x2)
            inten2Ey+=dQ_to_dx*miu2[k+1]*np.dot(t_state[i],x_y2z2)

        intensity_FC_Ex+=[intenEx]
        intensity_HT_Ex+=[inten2Ex]
        intensity_FC_Ey+=[intenEy]
        intensity_HT_Ey+=[inten2Ey]
        E+=[Et[i]-Ee]
    #print("here is end of process!")
    return index,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey

def splitEtoT_mp(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,polarization=np.array([0.,0.,0.]),maxpp=100):
    intensity_FC_Ex=[]
    intensity_HT_Ex=[]
    intensity_FC_Ey=[]
    intensity_HT_Ey=[]
    E=[]

    #print("direction ",polarization)
    if np.linalg.norm(polarization)>1e-10:
        polarization/=np.linalg.norm(polarization)
        for i in range(len(miu1)):
            print(miu1[i],polarization,np.dot(miu1[i],polarization))
            miu1[i]=np.dot(miu1[i],polarization)*polarization
            miu2[i]=np.dot(miu2[i],polarization)*polarization

    num_thread=32
    lenth=int(round(maxpp/num_thread))
    print("lenth",lenth)
    with mp.Pool(processes=num_thread) as pool:
        result=[]
        for i in range(num_thread):
            pool.apply_async(single_splitEtoT, [Ee,e_state,Et[i*lenth:(i+1)*lenth],t_state[i*lenth:(i+1)*lenth],na,miu1,miu2,mixa,i],callback=result.append)
        pool.close()
        pool.join()
    
    for i in range(num_thread):
        for data in result:
            if data[0]==i:
                #print("data",data)
                E+=data[1]
                intensity_FC_Ex+=data[2]
                intensity_HT_Ex+=data[3]
                intensity_FC_Ey+=data[4]
                intensity_HT_Ey+=data[5]

    return np.array(E),np.array(intensity_FC_Ex),np.array(intensity_HT_Ex),np.array(intensity_FC_Ey),np.array(intensity_HT_Ey)

def energy_add(a,b):
    return (np.tile(a,(len(b),1)).T+np.tile(b,(len(a),1))).reshape(1,-1)[0]

def direct_mul(a,b):
    la=len(a)
    lb=len(b)
    return np.matmul(a.reshape(la,1),b.reshape(1,lb))

def onlyEmode(miu1,miu2,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,process,polarization=np.array([0.,0.,0.])):
    print("onlyEmode!!!!!!!!!!!!!!!!!!")
    #FC term <QeQa|miu0|QeQa>
    inten_FC=direct_mul(intensity_FC_Ex,miu1[0])+direct_mul(intensity_FC_Ey,miu2[0])#miu0<QeQa|QeQa>

    #FC term <QeQa|miu(Qe)|QeQa>
    inten_HT_e=intensity_HT_Ex+intensity_HT_Ey#<Qe|miu(Qe)|Qe>

    print("shape of intens",inten_FC.shape,inten_HT_e.shape)
    inten=inten_FC+inten_HT_e
    
    total_inten=(inten)**2
    print("calculated polarization:",sum(total_inten.T[0]),sum(total_inten.T[1]),sum(total_inten.T[2]))
    total_inten=total_inten.T[0]+total_inten.T[1]+total_inten.T[2]
    
    #if np.linalg.norm(polarization)>0.1:
        #polar_inten=np.matmul(inten,polarization)
        #polar_inten=(inten)**2
    #else:   
    polar_inten=copy.deepcopy(total_inten)

    return E,total_inten,polar_inten,inten_FC,inten_HT_e


def combineA1mode(miu1,miu2,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,process,polarization=np.array([0.,0.,0.]),**para):
    print("combineA1mode!!!!!!!!!!!!!!!!!!")
    dmiu1=para.get("dmiu1")
    a1x=phonon(0,0.11/0.04,0.15,TDM=np.array(dmiu1[0]))
    print(np.linalg.norm(miu1[0]))
    dmiu2=para.get("dmiu2")
    a1y=phonon(1,0.11/0.04,0.15,TDM=np.array(dmiu2[0]))
    print(np.linalg.norm(miu2[0]))

    if "emission" in process:
        E=energy_add(E,-(a1x.FCenergy))
    elif "absorption" in process:
        E=energy_add(E,a1x.FCenergy)
    #FC term <QeQa|miu0|QeQa>
    intenx=np.kron(intensity_FC_Ex,a1x.FC)
    inteny=np.kron(intensity_FC_Ey,a1y.FC)#<QeQa|QeQa>
    inten_FC=direct_mul(intenx,miu1[0])+direct_mul(inteny,miu2[0])#miu0<QeQa|QeQa>
    #HT term <QeQa|miu(Qa)|QeQa>
    inten_HT_ax=np.kron(intensity_FC_Ex,a1x.HT)
    inten_HT_ay=np.kron(intensity_FC_Ey,a1y.HT)#<QeQa|Qa|QeQa>
    inten_HT_a=direct_mul(inten_HT_ax,a1x.TDM)+direct_mul(inten_HT_ay,a1y.TDM)#miuQa<QeQa|Qa|QeQa>
    #HT term <QeQa|miu(Qe)|QeQa>
    intensity_HT_E=intensity_HT_Ex+intensity_HT_Ey#<Qe|miu(Qe)|Qe>
    inten_HT_e=np.kron(intensity_HT_E,a1x.FC.reshape(len(a1x.FC),1))#<QaQe|miu(Qe)|QeQa>

    #test term miu=miuQaQe
    if para.get("AE",False):
        intensity_HTa_Ex=para.get("intensity_HTa_Ex")
        intensity_HTa_Ey=para.get("intensity_HTa_Ey")
        print("inten_HTa",intensity_HTa_Ex,"\n",intensity_HTa_Ey)
        intensity_HT_AE=intensity_HTa_Ex+intensity_HTa_Ey
        inten_HT_ae=np.kron(intensity_HT_AE,a1x.FC.reshape(len(a1x.FC),1))#<QaQe|miu(Qe)|QeQa>
        print("AE_EA")
    else:
        inten_HT_ae=0*copy.deepcopy(inten_HT_e)

    print("shape of intens",inten_FC.shape,inten_HT_a.shape,inten_HT_e.shape,inten_HT_ae.shape)
    inten=inten_FC+inten_HT_a+inten_HT_e+ inten_HT_ae
    total_inten=(inten)**2
    print("calculated polarization:",sum(total_inten.T[0]),sum(total_inten.T[1]),sum(total_inten.T[2]))
    total_inten=total_inten.T[0]+total_inten.T[1]+total_inten.T[2]
    print("shape of total_inten",total_inten.shape)
    if np.linalg.norm(polarization)>0.1:
        polar_inten=np.matmul(inten,polarization)
        polar_inten=(inten)**2
    else:
        polar_inten=copy.deepcopy(total_inten)
    return E,total_inten,polar_inten,inten_FC,inten_HT_a,inten_HT_e


def gaussian(omega, omega_k, sigma):#归一化高斯函数
    return 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(omega - omega_k)**2 / sigma**2 / 2)

def step(x):
    if x>0:
        return 1.0
    else: 
        return 0.0

def Spectra_single(E,I,sigma,x,hwe,offset,index):
    spectra=[]
    per=0.0
    for i in range(len(x)):
        xi=x[i]
        if(float(i)/len(x)>per):
            per+=0.1
        spec=0.0
        
        for i in range(len(I)):
            spec+=I[i]*gaussian(xi, E[i]*hwe+offset, sigma)#*(-0.5*E[i]+1+ 1*step(-E[i]-0.5) ) )#*(-1*E[i]+1+ 1*step(-E[i]-0.5) ) )
        spectra+=[spec] 

    return spectra,index

def Spectra(E,I,sigma,x,hwe,offset):
    print("convolve with gaussian start!")
    spectra=[]
    num_thread=32
    lenth=math.ceil(len(x)/num_thread)
    print("lenth",len(x)/num_thread,lenth)
    pool=mp.Pool(processes=num_thread)
    result=[]
    for i in range(num_thread):
        pool.apply_async(Spectra_single, [E,I,sigma,x[i*lenth:(i+1)*lenth],hwe,offset,i],callback=result.append)
    pool.close()
    pool.join()

    for j in range(num_thread):
        for data in result:
            if j==data[1]:
                spectra+=data[0]
    
    print("convolve with gaussian OK!")
    return np.array(spectra)

def vibplot(eigen):
    x=np.linspace(-5,5,100)
    y=np.linspace(-5,5,100)
    X,Y=np.meshgrid(x,y)
    out=np.zeros(len(x)*len(y),dtype=float).reshape(len(x),len(y))
    for i in range(len(x)):
        for j in range(len(y)):
            out[i][j]=vibfunc(eigen,x[i],y[j])**2
    #print(out)
    #Z=vibfunc(eigen,X,Y)
    fig, ax = plt.subplots(figsize=(8,8),dpi=300)
    cs=ax.contourf(X, Y, out, 30, alpha=.75, cmap='jet')
    fig.colorbar(cs)
    fig.show()
   # cs=ax.contour(x, y, eigen, 30, colors='black',)
   
    
def vibfunc(eigen,x,y):
    value=0
    for i in range(48):
        n,l=nl(i)
        nx=(n+l)/2
        ny=(n-l)/2
        value+=eigen[i]*hermite.vibration_wave_function(nx,x)*hermite.vibration_wave_function(ny,y)
    return value


def jtTtosplitE(kt,ke,**para):
    na=para.get("na",40)

    maxp=para.get("maxp",1000)

        
    miux=para.get("miux")
    miuy=para.get("miuy")

    dmiux=para.get("dmiux",np.array([]))
    dmiuy=para.get("dmiuy",np.array([]))
    mixa=para.get("mixa")
    Et,at=jtSolverTOringinalBase_New(kt,na,hwx=para.get("hwx",0),hwy=para.get("hwy",0),ratio=para.get("ratio",[1,1]))
    
    
    Ee,ae=jtSolver_Split_OringinalBase_New(ke,na,Sx=para.get("Sx",0),Sy=para.get("Sy",0),\
                                         hwx=para.get("hwx",1),hwy=para.get("hwy",1),v=para.get("v"))
    
        
    print("mix(1,x)(1,y)(2,x)(2,y)",sum((ae[0]**2)[:int(len(ae[0])/2)]),sum((ae[0]**2)[int(len(ae[0])/2):]),sum((ae[1]**2)[:int(len(ae[1])/2)]),sum((ae[1]**2)[int(len(ae[1])/2):]))
    np.savetxt("T.DAT",ae)    

    if np.linalg.norm(para.get("direction")) <1e-9:
        normal=1.0
    else:
        normal=np.linalg.norm(para.get("direction"))
    direction=para.get("direction",[.0,.0,.0])/normal

    if "emission" == para.get("process","emission"):  
        print("emission")
        E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey=TtosplitE_mp(Ee,ae,Et[0],at[0],na,miux,miuy,mixa,polarization=direction,maxpp=maxp)
        print("************len(E)=",len(E))
        E,IA,IALL,inten_FC,inten_HT_e=onlyEmode(miux,miuy,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,"emission")
        print("min,max",min(E),max(E))
        E-=max(E)
        print("min,max",min(E),max(E))
    elif "absorption" == para.get("process","emission"):
        print("absorption")
        E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey=splitEtoT_mp(Ee[0],ae[0],Et,at,na,miux,miuy,mixa,polarization=direction,maxpp=maxp)
        print("************len(E)=",len(E))
        E,IA,IALL,inten_FC,inten_HT_e=onlyEmode(miux,miuy,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,"absorption")
        E-=min(E)
    elif "absorption_and_A" == para.get("process","emission"):
        print("absorption_and_A")
        if len(dmiux)==0 or len(dmiuy)==0:
            exit("dmiux is empty!")
        E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey=splitEtoT_mp(Ee[0],ae[0],Et,at,na,miux,miuy,mixa,polarization=direction,maxpp=maxp)
        E,FCa_Ex,HTa_Ex,FCa_Ey,HTa_Ey=splitEtoT_mp(Ee[0],ae[0],Et,at,na,dmiux,dmiuy,mixa,polarization=direction,maxpp=maxp)
        E,IALL,inten_FC,inten_HT_a,inten_HT_e=combineA1mode(miux,miuy,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,"absorption_and_A",AE=True,intensity_HTa_Ex=HTa_Ex,intensity_HTa_Ey=HTa_Ey,dmiu1=dmiux,dmiu2=dmiuy)
        E-=min(E)
       
    elif "emission_and_A" == para.get("process","emission"):
        print("absorption_and_A")
        if len(dmiux)==0 or len(dmiuy)==0:
            exit("dmiux is empty!")
        E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey=TtosplitE_mp(Ee,ae,Et[0],at[0],na,miux,miuy,mixa,polarization=direction,maxpp=maxp)
        E,FCa_Ex,HTa_Ex,FCa_Ey,HTa_Ey=TtosplitE_mp(Ee,ae,Et[0],at[0],na,dmiux,dmiuy,mixa,polarization=direction,maxpp=maxp)
        E,IALL,inten_FC,inten_HT_a,inten_HT_e=combineA1mode(miux,miuy,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,"emission_and_A",AE=True,intensity_HTa_Ex=HTa_Ex,intensity_HTa_Ey=HTa_Ey,dmiu1=dmiux,dmiu2=dmiuy)   
        E-=max(E)
    else:
        print("proc input error,proc=",para.get("process"))
    np.savetxt("INTEN.DAT",np.array([E,IALL]).T)
    #print(polar)
    #I=I_HT
    #polar=polar_HT
    """ 
    IALL=[]
    for i in range(len(E)):
        a=I[i]+I_FCHT[i]+I_HT[i]
        IALL+=[np.linalg.norm(a)]

        
    print("I(FC)=",sum(I),sum(I)/sum(IALL),"I(FCHT)=",sum(I_FCHT),sum(I_FCHT)/sum(IALL),"I(HT)=",sum(I_HT),sum(I_HT)/sum(IALL))
    #print("Ix=",sum(Ix),"  Iy=",sum(Iy),"  Iz=",sum(Iz))

    print("direction=",direction)
    print("intensity=",sum(IALL))
    
    #print(Et[0],at[0])
  
    nn=np.linspace(0,len(at[0])-1,len(at[0]))
    plt.figure(figsize=(16, 8),dpi=300)
    plt.stem(nn,at[0])
    plt.show()
    
    vibplot(ae[0])
    vibplot(ae[1])
    vibplot(ae[2])
    vibplot(ae[3])
    vibplot(ae[4])
    vibplot(ae[5])
    print(Ee[:6])
    """
    
    
    xmin=E[0]-para.get("xmin",20)*para.get("hwx",1)
    xmax=E[0]+para.get("xmax",4)*para.get("hwx",1)
    
    
    res=para.get("resolution",0)
    if res>0:
        sigma=para.get("sigma",0.001)
        hwe=para.get("hwx",1)
        x=np.linspace(round(xmin*hwe,3),round((xmax)*hwe,3),int((round(xmax*hwe,3)-round(xmin*hwe,3))*res))
        xcm=x/0.0001204
        offset=0
        
        spectra=Spectra(E[0:maxp], IALL[0:maxp], sigma, x, hwe, offset)
        
        plt.figure(figsize=(16, 8),dpi=300)
        plt.plot(x,spectra)
        plt.title("T->E spectral function $A_e$  ($k_t=$"+str(kt)+" $\phi$="+para.get("thetastr","0")+"$\pi$"+" $k_e=$"+str(ke)+") ge="+str(para.get("Ge",[0])[0]))#+"   $k_e=$"+str(ke)+")")
        #n=spectra.tolist().index(max(spectra))
        plt.xlabel("phonon energy ($eV)")
        plt.savefig("TtoE.eps", dpi=300)
        plt.show()
    
    return E, IALL

def MnJT(kt,ke,**para):
    #Ae=jtTtoE(-2.4,3.0,sigma=2e-3,na=50,maxm=25,hwe=0.04,resolution=100,theta=2/3*np.pi,thetastr="2/3",xmax=2,xmin=40)
    sx=para.get("Sx",0)
    sy=para.get("Sy",0)
    Maxp=para.get("maxp",100)
    Na=para.get("na",30)
    hwX=para.get("hwx",1)
    hwY=para.get("hwy",1)
    ge=para.get("ge",[0.0,0.0])
    
     ###各T态到近简并E##


    E,I=jtTtosplitE(kt,ke,Sx=sx,na=Na,hwx=1.0,hwy=hwY/hwX,resolution=0,\
                      theta=para.get("theta",0)/3*np.pi,thetastr=str(para.get("theta",0)),maxp=Maxp,Ge=ge,\
                      ratio=para.get("ratio",[1,1]),v=para.get("v",[0.0,0.0]),\
                          process=para.get("process","emission"),mixa=para.get("mixa",[1.0,0.0]),\
                          miux=para.get("miux"), miuy=para.get("miuy"),direction=para.get("direction"),dmiux=para.get("dmiux",[]),dmiuy=para.get("dmiuy",[]))
    

    ##lancozs stem##
    print("split(cm-1)",(E[0]-E[1])*hwX/0.000124)
    print("amplite I/J I and J",I[0]/I[1],I[0],I[1])
    
    for i in range(len(I)):
        if I[i]>0.01 and I[i+1]>0.01:
            print("sasafddgfvdsbgsw")
            #E[i+1]=(2*E[i]+E[i+1])/3
    #E[24]=E[23]
            

    
    hwe=0.5*(hwX+hwY)
    xmin=para.get("xmin",(min(E)-3)*hwe)
    xmax=para.get("xmax",(max(E)+3)*hwe)
    res=para.get("resolution",100/hwe)
    sigma=para.get("sigma",1e-3*hwe)
    x=np.linspace(round(xmin,3),round((xmax),3),int((round(xmax,3)-round(xmin,3))*res))
    spectra=Spectra(E, I, sigma, x, hwX, 0)
    
    print("JT ok!\n intensity=",sum(I))
    

    #"S(hw) and Sk plot"
    fig0,ax0 = plt.subplots(1,1,figsize=(6.5,3.5),dpi=500)
    fig,ax,note=para.get("plot", (fig0,ax0,""))

    #plt.rcParams.update({'font.size': 20})
    xcm=x/0.000124
    a1=ax.plot(xcm,spectra,linewidth=1.2,label="A$_e$ ")
    #ax.set_title("$Cs_2SO_4:Mn^{6+} $ A$_e$($\hbar\omega$)") 
    ax.set_ylabel('A$_e$ abd stickes',fontsize=15)
    #ax.set_xlabel('energy ($cm^{-1}$) ' )#+" width="+str(p.delta_width))
    if "emission" in para.get("process"):
        ax.set_xlim(xmin*0.5/0.000124 , xmax*0.5/0.000124)
    elif "absorption" in para.get("process"):
        ax.set_xlim(xmin*0.5/0.000124 , xmax*0.5/0.000124)
    ax.set_ylim(min(spectra)*1.2-(max(spectra)-min(spectra))*0.3,max(spectra)*1.2)
    #ax.set_ylim(-0.5,3.0)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    #ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.set_yticks([])
    ax2=ax.twinx()
    Eev=[e*hwX/0.000124 for e in E]
    markerline, stemlines, baseline=ax2.stem(Eev,np.array(I),label="stikers",basefmt="k-") 
    markerline.set_markerfacecolor("tan")
    #markerline.set_size(1)
    stemlines.set_color("tan")
    #ax2.set_ylabel("stikers")
    ax2.set_ylim(0,max(I)*2)
    ax2.set_yticks([])

    ax2.legend(loc="upper right",frameon=False,fontsize=15)
    ax.legend(loc="upper left",frameon=False,fontsize=15)
    #ax2.set_ylim(0,0.008)
    #fig.legend(loc=[0.15,0.6],frameon=false)
    plt.savefig("JTspectra.eps", dpi=500,bbox_inches='tight')
    plt.savefig("JTspectra.png", dpi=500,bbox_inches='tight')
    #plt.show()
    
    

    return spectra,sum(I),E,I



    
    