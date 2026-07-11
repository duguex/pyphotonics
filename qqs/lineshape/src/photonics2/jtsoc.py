# -*- coding: utf-8 -*-
"""DEPRECATED/DORMANT — Jahn-Teller + spin-orbit coupling (SOC) solver.

Status: dormant. Not maintained, not imported by any active code path
(verified 2026-07-11). Import path fixed (`photonics` → `photonics2`)
but module still depends on `sympy` and `joblib`, which are not
installed. Future JT+SOC development would start here.

Originally written 2022-09-01 by original repo author (87930).
"""
#demensionless jt solver
from cProfile import label
from unittest import result
from joblib import PrintTime
import numpy as np
import matplotlib.pyplot as plt
import math
import multiprocessing as mp
import copy
import concurrent.futures as concurrent
import multiprocessing as mp
#from photonics2.ht import band,phonon
from sympy import LambertW, false, minimal_polynomial
#import photonics2.hermite as hermite
from matplotlib.pyplot import MultipleLocator, tight_layout

from photonics2.xyz import XYZ

def nl(k):
    #print("nl_input_k",k)
    n=round(math.ceil(np.sqrt(2*k+2.25)-1.5))
    l=round(2*k-(n+2)*n)
    return n,l

def kk(n,l):
    return round(((n+2)*n+l)/2)


def jtSolverTOringinalBase_New(kt,na,**para):
    print("================T*e solver module start!======================")

    ratio=para.get("ratio",[1.0,1.0])
    gex=(ratio[0]**2-1)/2
    gey=(ratio[1]**2-1)/2    
    
    hwx=para.get("hwx",1.0)
    hwy=para.get("hwy",1.0)
    dim=round((na+1)*(na+2)/2)
    print("#######################",na,dim)
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


    #print("OLD MAT0\n",mat0)
    E,a=np.linalg.eig(mat0)
    a=a.T
    E=np.array([e.real for e in E])
    sort=np.argsort(E)
    E=np.array([E[i] for i in sort])
    #print(E.tolist())
    a=np.array([a[i] for i in sort])
    print("================T*e solver module End!======================")
    return E,a


def jtSolverTOringinalBase_SOC(kt,na,lam,**para):#kt=[kt_yz,kt_zx,kt_xy]
    print("================T*e soc solver module start!======================")

    ratio=para.get("ratio",[[1.0,1.0],[1.0,1.0],[1.0,1.0]])
    print("ratio",ratio)
    hwx=para.get("hwx",[1.0,1.0,1.0])
    hwy=para.get("hwy",[1.0,1.0,1.0])
    print("hwx,hwy",hwx,hwy)
    print("soc lambda",lam)
    dim=round((na+1)*(na+2)/2)

    mat0=np.zeros(dim**2*9,dtype=np.complex64).reshape(dim*3,dim*3)
    
    
    
    for k in range(dim*3):
        n,x,y=nxy(k,na)
        #print(n,x,y)
        #(soc coupling constant) lambda/hwx
        #print(n,x,y)
        mat0[k][k]+=(x+0.5)*hwx[n]+(y+0.5)*hwy[n]
        write_mat(mat0, 2, 0, x, y, x, y, na, (ratio[n][0]**2-1)/2, n, n)
        write_mat(mat0, 0, 2, x, y, x, y, na, (ratio[n][1]**2-1)/2, n, n)
        #write H_jt
        for key in kt[n]:
            write_mat(mat0,kt[n][key].nx,kt[n][key].ny,x,y,x,y,na,kt[n][key].k,n,n)
        #write H_soc
        """
        H_SOC=
        YZ          i*lamba/2   -lamba/2
        -i*lamba/2      ZX      i*lamba/2
        -lamba/2    -i*lamba/2    XY
        """
        if n==0: #yz
            #print("yz")
            mat0[k,kkk(1,x,y,na)]+=complex(0,lam/2)
            mat0[k,kkk(2,x,y,na)]+=-lam/2
    
        elif n==1: #zx
            #print("zx")
            mat0[k,kkk(0,x,y,na)]+=complex(0,-lam/2)
            mat0[k,kkk(2,x,y,na)]+=complex(0,lam/2)

        elif n==2: #xy
            #print("xy")
            mat0[k,kkk(0,x,y,na)]+=-lam/2
            mat0[k,kkk(1,x,y,na)]+=complex(0,-lam/2)
        
 
    #print("NEW MAT0\n",mat0)
    E,a=np.linalg.eig(mat0)
    a=a.T
    E=np.array([e.real for e in E])
    sort=np.argsort(E)
    E=np.array([E[i] for i in sort])
    #print(E.tolist())
    a=np.array([a[i] for i in sort])
    print("================T*e soc solver module End!======================")
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
    print("here is process!")
    print(e_state.shape,t_state.shape)
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
        x_state=XpolydotPhi(e_state[0],1,0,na)
        x_x2=ax*x_state[0:mid]+ay*x_state[mid:]
        inten2Ex=copy.deepcopy(miu1[0]*np.dot(t_state,x_x2))*0.0
        inten2Ey=copy.deepcopy(miu1[0]*np.dot(t_state,x_x2))*0.0
        for k in range(len(miu2)-1):
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
    print("here is end of process!")
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
    #print("is that right?")
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
        #print("i",i)
        ##### FC term #####
        x2=ax*e_state[0:mid]+ay*e_state[mid:]
        y2z2=-ay*e_state[0:mid]+ax*e_state[mid:]
        intenEx=np.dot(t_state[i],x2)
        intenEy=np.dot(t_state[i],y2z2)
        #### HT term ####
        x_state=XpolydotPhi(e_state,1,0,na)
        x_x2=ax*x_state[0:mid]+ay*x_state[mid:]
        inten2Ex=copy.deepcopy(miu1[0]*np.dot(t_state[i],x_x2))*0.0
        inten2Ey=copy.deepcopy(miu1[0]*np.dot(t_state[i],x_x2))*0.0
        for k in range(len(miu2)-1):
            
            
            #print(nl(k+1))
            n,l=nl(k+1)
            nx=round((n-l)/2)
            ny=round((n+l)/2)
            dQ_to_dx=0.319876**n
            x_state=XpolydotPhi(e_state,nx,ny,na)
            
            x_x2=ax*x_state[0:mid]+ay*x_state[mid:]
            x_y2z2=-ay*x_state[0:mid]+ax*x_state[mid:]
            #print(dQ_to_dx*miu1[k+1],np.dot(t_state[i],x_x2))
            #print(dQ_to_dx*miu1[k+1]*np.dot(t_state[i],x_x2),(dQ_to_dx*miu1[k+1]*np.dot(t_state[i],x_x2)).shape)
            
            inten2Ex+=dQ_to_dx*miu1[k+1]*np.dot(t_state[i],x_x2)
            #print("k",k)
            inten2Ey+=dQ_to_dx*miu2[k+1]*np.dot(t_state[i],x_y2z2)
            

        intensity_FC_Ex+=[intenEx]
        intensity_HT_Ex+=[inten2Ex]
        intensity_FC_Ey+=[intenEy]
        intensity_HT_Ey+=[inten2Ey]
        E+=[Et[i]-Ee]
    #print("here is end of process!")
    #print("shape intensity_FC_Ex",np.array(intensity_FC_Ex).shape)
    return index,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey

def single_splitEtoT_soc(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,index):
    dim=round((na+1)*(na+2)/2)
    #print("test_dim",dim)
    print("single_splitEtoT_soc")
    print(t_state.shape)
    index,E,FC_Ex0,HT_Ex0,FC_Ey0,HT_Ey0=single_splitEtoT(Ee,e_state,Et,t_state[:,0*dim:(0+1)*dim],na,miu1[0],miu2[0],mixa,index)
    index,E,FC_Ex1,HT_Ex1,FC_Ey1,HT_Ey1=single_splitEtoT(Ee,e_state,Et,t_state[:,1*dim:(1+1)*dim],na,miu1[1],miu2[1],mixa,index)
    index,E,FC_Ex2,HT_Ex2,FC_Ey2,HT_Ey2=single_splitEtoT(Ee,e_state,Et,t_state[:,2*dim:(2+1)*dim],na,miu1[2],miu2[2],mixa,index)
    #print("shape",np.array(FC_Ex0).shape,np.array(FC_Ex1).shape,np.array(FC_Ex2).shape)
    #print("xy miu",miu1[2],miu2[2])
    #print("FC_Ex2",np.array(FC_Ex2),"HT_Ex2",np.array(HT_Ex2),np.array(FC_Ey2),np.array(HT_Ey2))
    intensity_FC_Ex=direct_mul(np.array(FC_Ex0),miu1[0][0])+direct_mul(np.array(FC_Ex1),miu1[1][0])+direct_mul(np.array(FC_Ex2),miu1[2][0])
    intensity_HT_Ex=np.array(HT_Ex0)+np.array(HT_Ex1)+np.array(HT_Ex2)
    intensity_FC_Ey=direct_mul(np.array(FC_Ey0),miu2[0][0])+direct_mul(np.array(FC_Ey1),miu2[1][0])+direct_mul(np.array(FC_Ey2),miu2[2][0])
    intensity_HT_Ey=np.array(HT_Ey0)+np.array(HT_Ey1)+np.array(HT_Ey2)
    #print("reduce shape",intensity_FC_Ex.shape)
    
    return index,E,intensity_FC_Ex.tolist(),intensity_HT_Ex.tolist(),intensity_FC_Ey.tolist(),intensity_HT_Ey.tolist()

def single_TtosplitE_soc(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,index):
    dim=round((na+1)*(na+2)/2)
    #print("test_dim",dim)
    print("single_TtosplitE_soc")
    print(t_state.shape)
    index,E,FC_Ex0,HT_Ex0,FC_Ey0,HT_Ey0=single_TtosplitE(Ee,e_state,Et,t_state[0*dim:(0+1)*dim],na,miu1[0],miu2[0],mixa,index)
    index,E,FC_Ex1,HT_Ex1,FC_Ey1,HT_Ey1=single_TtosplitE(Ee,e_state,Et,t_state[1*dim:(1+1)*dim],na,miu1[1],miu2[1],mixa,index)
    index,E,FC_Ex2,HT_Ex2,FC_Ey2,HT_Ey2=single_TtosplitE(Ee,e_state,Et,t_state[2*dim:(2+1)*dim],na,miu1[2],miu2[2],mixa,index)
    #print("shape",np.array(FC_Ex0).shape,np.array(FC_Ex1).shape,np.array(FC_Ex2).shape)
    #print("xy miu",miu1[2],miu2[2])
    #print("FC_Ex2",np.array(FC_Ex2),"HT_Ex2",np.array(HT_Ex2),np.array(FC_Ey2),np.array(HT_Ey2))
    intensity_FC_Ex=direct_mul(np.array(FC_Ex0),miu1[0][0])+direct_mul(np.array(FC_Ex1),miu1[1][0])+direct_mul(np.array(FC_Ex2),miu1[2][0])
    intensity_HT_Ex=np.array(HT_Ex0)+np.array(HT_Ex1)+np.array(HT_Ex2)
    intensity_FC_Ey=direct_mul(np.array(FC_Ey0),miu2[0][0])+direct_mul(np.array(FC_Ey1),miu2[1][0])+direct_mul(np.array(FC_Ey2),miu2[2][0])
    intensity_HT_Ey=np.array(HT_Ey0)+np.array(HT_Ey1)+np.array(HT_Ey2)
    #print("reduce shape",intensity_FC_Ex.shape)
    
    return index,E,intensity_FC_Ex.tolist(),intensity_HT_Ex.tolist(),intensity_FC_Ey.tolist(),intensity_HT_Ey.tolist()

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



def splitEtoT_mp_soc(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,polarization=np.array([0.,0.,0.]),maxpp=100):
    intensity_FC_Ex=[]
    intensity_HT_Ex=[]
    intensity_FC_Ey=[]
    intensity_HT_Ey=[]
    E=[]
    print("plitEtoT_mp miu1",miu1)
    print("plitEtoT_mp miu2",miu2)

    #print("direction ",polarization)
    if np.linalg.norm(polarization)>1e-10:
        polarization/=np.linalg.norm(polarization)
        for j in range(3):
            for i in range(len(miu1)):
                print(miu1[j][i],polarization,np.dot(miu1[j][i],polarization))
                miu1[j][i]=np.dot(miu1[j][i],polarization)*polarization
                miu2[j][i]=np.dot(miu2[j][i],polarization)*polarization

    num_thread=32
    lenth=int(round(maxpp/num_thread))
    print("lenth",lenth)
    with mp.Pool(processes=num_thread) as pool:
        result=[]
        for i in range(num_thread):
            pool.apply_async(single_splitEtoT_soc, [Ee,e_state,Et[i*lenth:(i+1)*lenth],t_state[i*lenth:(i+1)*lenth],na,miu1,miu2,mixa,i],callback=result.append)
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
    
    print(np.array(E).shape,np.array(intensity_FC_Ex).shape,np.array(intensity_HT_Ex).shape,np.array(intensity_FC_Ey).shape,np.array(intensity_HT_Ey).shape)

    return np.array(E),np.array(intensity_FC_Ex),np.array(intensity_HT_Ex),np.array(intensity_FC_Ey),np.array(intensity_HT_Ey)


def TtosplitE_mp_soc(Ee,e_state,Et,t_state,na,miu1,miu2,mixa,polarization=np.array([0.,0.,0.]),maxpp=100):
    intensity_FC_Ex=[]
    intensity_HT_Ex=[]
    intensity_FC_Ey=[]
    intensity_HT_Ey=[]
    E=[]
    print("plitEtoT_mp miu1",miu1)
    print("plitEtoT_mp miu2",miu2)

    #print("direction ",polarization)
    if np.linalg.norm(polarization)>1e-10:
        polarization/=np.linalg.norm(polarization)
        for j in range(3):
            for i in range(len(miu1)):
                print(miu1[j][i],polarization,np.dot(miu1[j][i],polarization))
                miu1[j][i]=np.dot(miu1[j][i],polarization)*polarization
                miu2[j][i]=np.dot(miu2[j][i],polarization)*polarization

    num_thread=32
    lenth=int(round(maxpp/num_thread))
    print("lenth",lenth)
    with mp.Pool(processes=num_thread) as pool:
        result=[]
        for i in range(num_thread):
            pool.apply_async(single_TtosplitE_soc, [Ee[i*lenth:(i+1)*lenth],e_state[i*lenth:(i+1)*lenth],Et,t_state,na,miu1,miu2,mixa,i],callback=result.append)
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
    
    print(np.array(E).shape,np.array(intensity_FC_Ex).shape,np.array(intensity_HT_Ex).shape,np.array(intensity_FC_Ey).shape,np.array(intensity_HT_Ey).shape)

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

def onlyEmode_soc(miu1,miu2,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,process,polarization=np.array([0.,0.,0.])):
    print("onlyEmode!!!!!!!!!!!!!!!!!!")
    #FC term <QeQa|miu0|QeQa>
    inten_FC=intensity_FC_Ex+intensity_FC_Ey#miu0<QeQa|QeQa>

    #FC term <QeQa|miu(Qe)|QeQa>
    inten_HT_e=intensity_HT_Ex+intensity_HT_Ey#<Qe|miu(Qe)|Qe>

    print("shape of intens",inten_FC.shape,inten_HT_e.shape)
    inten=inten_FC+inten_HT_e
    
    total_inten=inten*inten.conjugate()
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
   #print("11111111111111",Et,)
    if "soc" in para.get("process","emission"):
        Et,at,MIU1,MIU2=temp_soc(para.get("lam",0.0),na,para.get("test",1))
    else:
        Et,at=jtSolverTOringinalBase_New(kt,na,hwx=para.get("hwx",0),hwy=para.get("hwy",0),ratio=para.get("ratio",[1,1]))
    #print("22222222222222",Et)
    
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
    elif "absorption_soc" == para.get("process","emission"):
        print("absorption_soc")
        E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey=splitEtoT_mp_soc(Ee[0],ae[0],Et,at,na,MIU1,MIU2,mixa,polarization=direction,maxpp=maxp)
        print("\n************\n",E,"\n************\n",Ee,"\n************\n",Et,"\n************\n")
        print("************len(E)=",len(E))
        E,IA,IALL,inten_FC,inten_HT_e=onlyEmode_soc(miux,miuy,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,"absorption_soc")
        E-=min(E)
    elif "emission_soc" == para.get("process","emission"):
        print("emission_soc")
        E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey=TtosplitE_mp_soc(Ee,ae,Et[0],at[0],na,MIU1,MIU2,mixa,polarization=direction,maxpp=maxp)
        print("************len(E)=",len(E))
        print("************\n",E[:20],"************\n",Ee[:20],"************\n",Et[:20],"************\n")
        E,IA,IALL,inten_FC,inten_HT_e=onlyEmode_soc(miux,miuy,E,intensity_FC_Ex,intensity_HT_Ex,intensity_FC_Ey,intensity_HT_Ey,"emission_soc")
        E-=max(E)
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
        print("proc input error,proc=",proc)
    
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
                      ratio=para.get("ratio",[1,1]),v=para.get("v",[0.0,0.0]),lam=para.get("lam",1.0),\
                          process=para.get("process","emission"),mixa=para.get("mixa",[1.0,0.0]),test=para.get("test",1),\
                          miux=para.get("miux"), miuy=para.get("miuy"),direction=para.get("direction"),dmiux=para.get("dmiux",[]),dmiuy=para.get("dmiuy",[]))
            
    print(E,I)
    ##lancozs stem##
    print("split(cm-1)",(E[0]-E[1])*hwX/0.000124)
    print("amplite I/J I and J",I[0]/I[1],I[0],I[1])
    

            

    
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
    lam=para.get("lam",1.0)
    ax.set_title(f"$\lambda={lam}E(Q_1)$")

    ax2.legend(loc="upper right",frameon=False,fontsize=15)
    ax.legend(loc="upper left",frameon=False,fontsize=15)
    #ax2.set_ylim(0,0.008)
    #fig.legend(loc=[0.15,0.6],frameon=false)
    plt.savefig("JTspectra.eps", dpi=500,bbox_inches='tight')
    plt.savefig("JTspectra.png", dpi=500,bbox_inches='tight')
    #plt.show()
    
    

    return spectra,sum(I),E,I


class aps:
    def __init__(self,k,nx,ny,**para):
        self.k=k
        self.nx=nx
        self.ny=ny
        if para.get("hwe",0)>0:
            amu=1.66053886E-27
            e=1.602176634E-19
            angstrom=1e-10
            hbar=1.05457266E-34
            Q=np.sqrt(amu)*angstrom
            #print("Q",Q)
            hwe_in_SI=para.get("hwe")*e
            k_in_SI=k*e/Q**(nx+ny)
            k_dimensionless=k_in_SI*hbar**(nx+ny)/hwe_in_SI**(0.5*nx+0.5*ny+1)
            #print("nx,ny,k,k_in_SI,k_dim",nx,ny,k,k_in_SI,k_dimensionless)
            self.k=k_dimensionless

def Cs(T):
    qx=0.04085
    qy=0.04093

    hw=qx
    
    w={}
    w["x"]=aps(-0.0791,1,0,hwe=hw)
    #w["x"]=aps(-0.624,1,0)
    w["x2"]=aps(-0.01501,2,0,hwe=hw)
    w["y2"]=aps(0.01201,0,2,hwe=hw)
    #w["x2"]=aps(-0.038,2,0)
    #w["y2"]=aps(0.038,0,2)
    #w["x3"]=aps(0.002,3,0,hwe=hw)
    #w["xy2"]=aps(0.0015,1,2,hwe=hw)
    
    v={}
    v["y"]=aps(-0.0708,0,1,hwe=hw)
    #v["y"]=aps(-0.553,0,1)
    v["xy"]=aps(0.028,1,1)
    #v["x2y"]=aps(0.001,2,1,hwe=hw)
    #v["y3"]=aps(0.0032,0,3,hwe=hw)
    
    
    diag={}
    diag["x3"]=aps(-0.012,3,0,hwe=hw)
    diag["xy2"]=aps(0.028,1,2,hwe=hw)
    #diag["x3"]=aps(-0.038,3,0)
    #diag["xy2"]=aps(0.038*3,1,2)
    #diag["x4"]=aps(5.29E-04,4,0)
    #diag["x2y2"]=aps(3.27E-03,2,2)
    #diag["y4"]=aps(-2.42E-03,2,2)
        

    if "yz" in T or "zx" in T:
        exqx=0.03997
        exqy=0.04094
        
        diag2={}
        diag2["x"]=aps(-0.21512,1,0,hwe=hw)
        diag2["y"]=aps(-0.37683,0,1,hwe=hw)
        diag2["xy"]=aps(0.01863,1,1,hwe=hw)
        diag2["x3"]=aps(-0.007661,3,0,hwe=hw)
        



        diag2["x2y"]=aps(0.0075,2,1,hwe=hw)
        diag2["xy2"]=aps(0.03875,1,2,hwe=hw)
        diag2["y3"]=aps(0.00652,0,3,hwe=hw)

        x2=[-0.08592,0.03444,-0.3178,-0.15365,-0.27984,-0.52495]
        z2=[0.10938,-0.09185,0.26181,0.08945,0.21194,0.39719]
        y2=[0.14196,-0.03642,0.34265,0.10917,0.2905,0.89176]
        x1=[-0.01336,-0.44519,-0.09491,0.10005,-0.26749,0.02826]
        z1=[0.05177,0.35201,0.06572,-0.08301,0.21263,0.01141]
        y1=[0.04362,0.65859,0.07906,-0.12145,0.39467,-0.14354]

        x2=[-0.09498,0.07426,-0.29843,-0.24411,-0.30453,-0.47361]
        z2=[0.11679,-0.12442,0.24596,0.16345,0.23213,0.35519]
        y2=[0.15442,-0.09112,0.31604,0.23344,0.32441,0.82122]
        x1=[-0.0069,-0.47358,-0.10872,0.16453,-0.2499,-0.00834]
        z1=[0.04838,0.36688,0.07295,-0.1168,0.20341,0.03059]
        y1=[0.03628,0.69085,0.09475,-0.19473,0.37468,-0.10195]


        #hs equal ex occ


        x2=[-1.47E-01,-0.16548,-0.43323,0.01039,-0.04912,-0.05784,-0.02221,-0.06437,-0.07574,-0.02108]
        z2=[0.07681,6.12E-02,0.34235,0.01936,0.03772,0.02987,0.00747,0.03068,0.06409,0.01905]
        y2=[0.15473,3.11E-01,0.4838,-0.27669,-0.04424,0.23621,0.19706,0.01295,0.20436,-0.07229]
        x1=[-0.08096,-2.61E-01,-0.19654,-0.01199,-0.03886,-0.00187,0.01471,0.01915,-0.01374,-0.05335]   
        z1=[-1.34E-02,1.51E-01,0.08373,0.00688,0.11166,0.03874,-0.00605,-0.0264,-0.01318,0.0364]    
        y1=[0.05771,0.31189,0.15099,-0.01186,0.07643,0.04304,-0.0032,-0.02397,0.0199,0.04186] 
        #exhs minus zero
        dmiu1=np.array([x1,y1,z1]).T    
        dmiu2=np.array([x2,y2,z2]).T*-1

        x2=[0.09266,-0.25093,-0.5402,0.09827,0.0662,-0.04789,-0.05133,-0.1556,-0.14893,-0.03555]
        z2=[-0.03403,0.11772,0.41814,-0.04996,-0.04292,0.03182,0.02901,0.09462,0.10991,0.02788]
        y2=[-0.10754,0.29882,0.71541,-0.11369,-0.07673,-0.00363,0.07104,0.16628,0.21035,0.0679]
        x1=[-0.00621,-0.23387,-0.1708,-0.09308,-0.17228,0.00813,3.39E-02,0.07895,0.0108,-0.09231]
        z1=[0.02303,0.23474,0.1415,0.01663,0.12436,-0.01202,-0.01333,-0.04404,0.00107,0.0744]
        y1=[-0.02099,-0.35003,-0.10178,-0.0829,-0.18948,-0.04194,0.04976,0.05569,0.01906,-0.08273]
        #exhs exocc

        x2=[0.09014 ,-0.17832 ,-0.49187 ,0.03553 ,-0.01843 ,-0.08191 ,-0.03288 ,-0.13521 ,-0.11507 ,-0.03032 ]
        z2=[-0.02800 ,0.04746 ,0.39352 ,-0.00091 ,0.02345 ,0.06685 ,0.01354 ,0.05507 ,0.08858 ,0.01768 ]
        y2=[-0.10885 ,0.24223 ,0.68523 ,-0.05469 ,-0.00140 ,0.02273 ,0.05352 ,0.12696 ,0.19047 ,0.06713 ]
        x1=[0.00531 ,-0.22846 ,-0.14218 ,-0.08210 ,-0.17049 ,-0.02412 ,0.02641 ,0.09827 ,0.01668 ,-0.07388 ]
        z1=[0.00729 ,0.26111 ,0.15961 ,-0.00123 ,0.10301 ,-0.03957 ,-0.00892 ,-0.03298 ,0.00851 ,0.08338 ]
        y1=[-0.01571 ,-0.38198 ,-0.09259 ,-0.07864 ,-0.21052 ,-0.07224 ,0.05220 ,0.05726 ,0.02033 ,-0.06572 ]
        #hs exocc 

        

        x2=[-0.15543,-0.15173,-0.43562,-0.00209,-0.05705,-0.06139,-0.0189,-0.06837,-0.07082,-0.02148]
        z2=[0.07897,0.04721,0.34976,0.01211,0.0418,0.03943,0.01199,0.02331,0.06265,0.01704]
        y2=[0.1694,0.21281,0.51726,-1.99E-04,0.06566,0.09366,0.0225,0.06349,0.11011,0.02369]
        x1=[-0.08761,-0.25981,-0.16944,-0.01007,-0.02433,-0.03678,0.00958,0.03178,-0.01751,-0.03813]
        z1=[-0.01061,0.16553,0.06779,0.0188,0.10153,0.0629,-0.01353,-0.02706,-0.00673,0.02723]
        y1=[0.0559,0.31342,0.16452,0.0124,0.12018,0.01174,-0.02291,-0.00763,-0.00686,0.06059]
        #hs minus zero
        
        if "test" in T:
            x1=[1.0,0.0]
            y1=[0.0,0.0]
            z1=[1.0,0.0]
            x2=[-1.77,0.0]
            y2=[0.0,0.0]
            z2=[-1.77,0.0]





        #Csyz结果
    elif "xy" in T:
        exqx=0.04123
        exqy=0.04094
        
        diag2={}
        
        diag2["split"]=aps(0.06,0,0,hwe=hw)  
        diag2["x"]=aps(0.42788,1,0,hwe=hw)
        diag2["y"]=aps(0,0,1,hwe=hw)
        diag2["x3"]=aps(-0.01984,3,0,hwe=hw)
        diag2["xy2"]=aps(0.02261,1,2,hwe=hw)
       
        x2=[0.00037 ,-0.00136 ,0.14712 ,-0.00204 ,-0.06306 ,0.00002 ,-0.00089 ,-0.08935 ,0.00008 ,-0.00505 ]
        z2=[0.00075 ,-0.00249 ,0.26126 ,-0.00321 ,-0.06549 ,0.00004 ,-0.00134 ,-0.10763 ,0.00021 ,-0.01045 ]
        y2=[-0.04726 ,0.09632 ,0.00055 ,0.08154 ,0.00147 ,-0.01935 ,0.02754 ,0.00079 ,-0.01391 ,-0.00006 ]
        x1=[-0.25411 ,0.54830 ,-0.00001 ,-0.05780 ,0.00084 ,0.02113 ,0.06613 ,0.00046 ,0.00766 ,-0.00014 ]
        z1=[-0.24196 ,0.71373 ,-0.00006 ,-0.10070 ,0.00109 ,0.04005 ,0.08735 ,0.00060 ,0.00793 ,-0.00027 ]
        y1=[-0.00015 ,0.00058 ,0.03498 ,0.00065 ,0.03036 ,-0.00010 ,0.00022 ,0.00763 ,-0.00004 ,0.00008 ]
        #Csxy结果
        if "test" in T:
            x1=[0,0.0]
            y1=[2,0.0]
            z1=[0,0.0]
            x2=[0,0.0]
            y2=[0.0,0.0]
            z2=[0,0.0]

    if "zx" in T:
        print("mirror yz for zx")
        for key in diag2:
            if diag2[key].ny%2==1:
                diag2[key].k*=-1
        y2=[-a for a in y2]
        x1=[-a for a in x1]
        z1=[-a for a in z1]
        for i in range(len(x2)):
            n,x,y=nxy(i,4)
            if x%2==1:
                x2[i]*=-1
                y2[i]*=-1
                z2[i]*=-1
                x1[i]*=-1
                y1[i]*=-1
                z1[i]*=-1
        if "test" in T:
            x1=[1.0,0.0]
            y1=[0.0,00]
            z1=[-1.0,0.0]
            x2=[1.77,0.0]
            y2=[0.0,0.0]
            z2=[-1.77,0.0]

            

        

    
    miu1=np.array([x1,y1,z1]).T
    miu2=np.array([x2,y2,z2]).T*-1

    for i in range(1,len(miu1)):
        miu1[i]*=1
        miu2[i]*=1
    
    
    #dmiu1 = (dmiu1-miu1)/0.54
    #dmiu2 = (dmiu2-miu2)/0.54#0.54 for S=0.15 d(miu)/dQa

    #print("miu1",miu1)
    #print("miu2",miu2)
    #print("dmiu1",dmiu1)
    #print("dmiu2",dmiu2)

    return w,v,diag,diag2,miu1,miu2,qx,qy,exqx,exqy#,dmiu1,dmiu2
    


    
    
def temp_soc(lam,na,test):
    np.set_printoptions(
        linewidth=500,
        precision=1,
        suppress=True,
        threshold=1000,
        edgeitems=100,
    )
    kt=[]
    rat=[]
    QX=[]
    QY=[]
    MIU1=[]
    MIU2=[]
    w,v,diag,diag2,miu1,miu2,qx,qy,exQx,exQy=Cs("yz")
    w["x"].k*=0.85
    v["y"].k*=0.85
    diag2["x"].k*=0.85
    diag2["y"].k*=0.85
    kt.append(diag2)
    QX.append(qx)
    QY.append(qy)
    rat.append([exQx/qx,exQy/qx])
    MIU1.append((miu1*test).tolist())
    MIU2.append((miu2*test).tolist())
    w,v,diag,diag2,miu1,miu2,qx,qy,exQx,exQy=Cs("zx")
    w["x"].k*=0.85
    v["y"].k*=0.85
    diag2["x"].k*=0.85
    diag2["y"].k*=0.85
    kt.append(diag2)
    QX.append(qx)
    QY.append(qy)
    rat.append([exQx/qx,exQy/qx])
    MIU1.append((miu1*test).tolist())
    MIU2.append((miu2*test).tolist())
    w,v,diag,diag2,miu1,miu2,qx,qy,exQx,exQy=Cs("xy")
    w["x"].k*=0.85
    v["y"].k*=0.85
    diag2["x"].k*=0.85
    diag2["y"].k*=0.85
    kt.append(diag2)
    QX.append(qx)
    QY.append(qy)
    rat.append([exQx/qx,exQy/qx])
    MIU1.append((miu1*(1-test)).tolist())
    MIU2.append((miu2*(1-test)).tolist())
    print("MIU1,MIU2",MIU1,MIU2)
    print("AAAAAAAAAA",rat)
    dim=round((na+1)*(na+2)/2)
    kt[2]["split"]=aps((858-557)*0.000124,0,0,hwe=qx)  
    Et,at=jtSolverTOringinalBase_SOC(kt,na,lam,ratio=rat)
    print(at)
    for i,k in enumerate(kt):
        print("i=",i)
        for key in k:
            print(key,k[key].k,k[key].nx,k[key].ny) 

    #print("Et",Et)
    #print("at",at)
    data=[]
    y=[]
    x=[]

    for i,a in enumerate(at):
        if i>41:
            break
        print("i=",i,end="\t")
        print("Et",Et[i],end="\t")
        print("yz=",np.dot(a[:dim],a[:dim].conjugate()).real,"zx=",np.dot(a[dim:2*dim],a[dim:2*dim].conjugate()).real,"xy=",np.dot(a[2*dim:],a[2*dim:].conjugate()).real)
        data.append((np.dot(a[:dim],a[:dim].conjugate()).real))
        y.append(1)
        x.append(i)
        data.append((np.dot(a[dim:2*dim],a[dim:2*dim].conjugate()).real))
        y.append(2)
        x.append(i)
        data.append((np.dot(a[2*dim:3*dim],a[2*dim:3*dim].conjugate()).real))
        y.append(3)
        x.append(i)
        
        #print("max imag=",max(a.imag))
        if abs(Et[i]-Et[i+1])<0.1:
            print("split",Et[i]-Et[i+1])
    fig=plt.figure(figsize=(8,4),dpi=500)
    ax=fig.add_subplot()
    sc=ax.scatter(x,y,c=data,cmap="binary")
    plt.colorbar(sc,label="|<T|t>|$^2$")
    ax.set_yticks([1,2,3],["yz","zx","xy"])
    fig.savefig("soc.png")
    return Et,at,np.array(MIU1),np.array(MIU2)
    
    
    
