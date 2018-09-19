"""
Extrapolating Pk with matter Pk

"""

import numpy as np
import healpy as hp
import os,sys
import random
from astropy.io import fits
import scipy.optimize 
from scipy import stats
import cosmo_func as cosmo
import EoR_func as Ef
from scipy.interpolate import interp1d

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] indir')
    o.set_description(__doc__)
    o.add_option('--zmin',dest='zmin',default=8,help='zmin of the input Pk')
    o.add_option('--zmax',dest='zmax',default=8.500,help='zmax of the input Pk')
    o.add_option('--Dz',dest='Dz',default=0.125,help='Delta z of the input Pk')
    o.add_option('--outdir',dest='outdir',help='Path of Output extrapolated Pk file directory')
    o.add_option('--mPk_file',dest='mPk_file',default='matter_Pk.txt', help='filename of the input matter Pk')
    o.add_option('--single_pk',action='store_true',dest='single_pk',default=False,help='If True only on Pk is extrapolated: the one specified by zmin and zmax')

    opts, args = o.parse_args(sys.argv[1:]) 
    
    print opts, args 

    indir = args[0]  
    zmin,zmax,Dz = np.float(opts.zmin),np.float(opts.zmax),np.float(opts.Dz)
    
    outdir=str(opts.outdir)
    os.mkdir(outdir) #create out directory and check if exits #comment this line if you want to write on the same folder
    


    if (opts.single_pk==False): 
        print "generation z list"
        z_vec=np.arange(zmin,zmax+Dz,Dz)
        nz=len(z_vec)
    else: 
        print "generating single cl"
        z_vec=[zmin,zmax]
        nz=1   

    print "Reading Pks in" + indir
    #fileformat: pk_z#_z#.txt
    z_list,pk_vec=Ef.readPk_input(indir,z_vec,opts.single_pk)
    npk=len(pk_vec)
    print "Input pks read"

    print "reading matter Pk"
    mk=np.loadtxt(str(opts.mPk_file))[0,:] #units are Mpc-1
    mPk=np.loadtxt(str(opts.mPk_file))[1,:] #units are Mpc^3
   
    Pcamb = interp1d(mk,mPk) #create matter Pk func
    Dk=pk_vec[0,1,0]-pk_vec[0,0,0]
    k_low=np.arange(mk[0],pk_vec[0,1,0]+Dk/10,Dk/10) #refining k space: factor 10 can be changed
    ikmax=-199 #doing the kmax cut at 5
    k_high=pk_vec[0,1:ikmax,0]
    k_long=np.concatenate((k_low,k_high))
   

    print "extrapolating pk in "+ outdir
    pk_extrapol=np.zeros((len(z_list),len(k_long),2))
    index=0
    for z1 in range(len(z_vec)):
        for z2 in range(z1,len(z_vec)):
            pk_extrapol[index,0:len(k_low),0]=k_low
            pk_extrapol[index,len(k_low):len(k_long),0]=k_high
        
            pk_extrapol[index,0:len(k_low),1]=Pcamb(k_low)/Pcamb(pk_vec[index,1,0])*pk_vec[index,1,1]
            pk_extrapol[index,len(k_low):len(k_long),1]=pk_vec[index,1:ikmax,1]    
        
        
            filename=outdir +'pk_z'+str("%4.3f"% z_vec[z1])+'_z'+str("%4.3f"% z_vec[z2])+'.txt'
            print filename 
        
            np.savetxt(filename,pk_extrapol[index,:,:])
            index+=1

    print " done."
