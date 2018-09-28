#!/usr/bin/env python

import numpy as np
import healpy as hp
import random
import sys
import cosmo_func as cosmo
from scipy.special import spherical_jn
from scipy.integrate import trapz, quad
from scipy.interpolate import splrep, splev, splint
from scipy.interpolate import interp1d




#function to read input pk
def readPk_input(indir,z_vec,single_cl):
    pk_vec=[]
    z_list=[]
    if single_cl==False:
        for z1 in range(len(z_vec)):
            for z2 in range(z1,len(z_vec)):
                filename=indir +'pk_z'+str("%4.3f"% z_vec[z1])+'_z'+str("%4.3f"% z_vec[z2])+'.txt'
                z_list.append([z_vec[z1],z_vec[z2]])
                pk_vec.append(np.loadtxt(filename))
        pk_vec=np.array(pk_vec)
        z_list=np.array(z_list)
    else:
        filename=indir +'pk_z'+str("%4.3f"% z_vec[0])+'_z'+str("%4.3f"% z_vec[1])+'.txt'
        z_list.append([z_vec[0],z_vec[1]])
        pk_vec.append(np.loadtxt(filename))
        pk_vec=np.array(pk_vec)
        z_list=np.array(z_list)
    return z_list, pk_vec



######
#function for pk to cl intergation
#from eq. (6)-(7) of https://arxiv.org/pdf/astro-ph/0408515.pdf 

#Compute Bessel

def I_bes_r(l,nu0,k,Dnu,bres):
    r_f=cosmo.transverse_comoving_distance(cosmo.fq2z(nu0-Dnu/2))
    r_i=cosmo.transverse_comoving_distance(cosmo.fq2z(nu0+Dnu/2))
    R=np.linspace(r_i,r_f,bres)
    Dr=r_f-r_i
    #print(len(R),Dr)
    Ib=[]
    for kk in k:
        integrand=spherical_jn(l,kk*R)
        int_rep=splrep(R,integrand)
        Ib.append(splint(r_i,r_f,int_rep))
    
    return np.array(Ib)/Dr
    
def Pk2Cl_trapz_r(l,z1,z2,k_in,pk_in,Dnu,bres,kres):
    #intepolate Pk
    
    pk_interp=interp1d(k_in,pk_in)
    k=np.linspace(k_in[0],k_in[-1],kres)
    pk=pk_interp(k)
    integrand=pk*k**2*2/np.pi
    if (z1==z2):
        cl=trapz(integrand*I_bes_r(l,cosmo.z2fq(z1),k,Dnu,bres)**2,k) 
    else: 
        cl=trapz(integrand*I_bes_r(l,cosmo.z2fq(z1),k,Dnu,bres)*I_bes_r(l,cosmo.z2fq(z2),k,Dnu,bres),k) 
    return cl   

def Pk2Cl_trapz_r_log(l,z1,z2,k_in,pk_in,Dnu,bres,kres):
    #intepolate Pk
    
    pk_interp=interp1d(k_in,pk_in)
    k=np.logspace(np.log10(k_in[0]+1e-5),np.log10(k_in[-1]),kres)
    pk=pk_interp(k)
    integrand=pk*k**2*2/np.pi
    if (z1==z2):
        cl=trapz(integrand*I_bes_r(l,cosmo.z2fq(z1),k,Dnu,bres)**2,k,Dnu,bres),k) 
    else:
        cl=trapz(integrand*I_bes_r(l,cosmo.z2fq(z1),k,Dnu,bres)*I_bes_r(l,cosmo.z2fq(z2),k,Dnu,bres),k) 
    return cl   

#######




#function to generate the (lmax,nnu,nnu) matrix to use for correlations
def readCl_input(indir,z_vec,lmax):
    nz=len(z_vec)
    cl_matrix=np.zeros((nz,nz,lmax,2),dtype=float)
    for z1 in range(nz):
        for z2 in range(z1,nz):
            filename=indir +'Cl_z'+str("%4.3f"% z_vec[z1])+'_z'+str("%4.3f"% z_vec[z2])+'.txt'
            cl=np.loadtxt(filename)
            #check cl vs nside
            if len(cl[:,1])!=lmax: 
                 print "Cl length do not match nside"
                 sys.exit(1)
            cl_matrix[z1,z2,:,0]=cl[:,0]
            if z1!=z2: cl_matrix[z2,z1,:,0]=cl[:,0]
            cl_matrix[z1,z2,:,1]=cl[:,1]
            if z1!=z2: cl_matrix[z2,z1,:,1]=cl[:,1]
    
    return cl_matrix


#function to generate the sky cube
def genSKYcube_corr(nnu, nside,seed, cov):
   
    
    #cov is (lmax,nnu,nnu)
    
    ell=np.arange(3*nside,dtype=float)

    count=0
    lmax=3*nside
    alm_wo=np.zeros((lmax*(lmax+1)/2,nnu),dtype=complex)
    
    ##set seed 
    np.random.seed(seed)
    for l in range(lmax):
            
            alm_wo[count:count+l+1,:]=1./np.sqrt(2)*(np.random.multivariate_normal(np.zeros(nnu),cov[:,:,l],l+1)+1j*np.random.multivariate_normal(np.zeros(nnu),cov[:,:,l],l+1))
            count+=l+1
    
    #rearrange alm
    alm=np.zeros((lmax*(lmax+1)/2,nnu),dtype=complex)
    count=0
    for m in range(lmax):
        for l in range(m,lmax):
            alm[count,:]=alm_wo[l*(l+1)/2+m,:]
            count+=1
    cube=np.zeros((nnu,hp.nside2npix(nside)),dtype=float)
    for nu in range(nnu):
        alm_C=alm[:,nu]
        cube[nu,:]=hp.sphtfunc.alm2map(np.ascontiguousarray(alm_C), nside)
        

    return cube
