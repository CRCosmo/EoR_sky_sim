
"""
Generate Cl from Pk

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



  
if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] indir')
    o.set_description(__doc__)
    o.add_option('--zmin',dest='zmin',default=8,help='zmin of the input Pk')
    o.add_option('--zmax',dest='zmax',default=8.500,help='zmax of the input Pk')
    o.add_option('--Dz',dest='Dz',default=0.125,help='Delta z of the input Pk')
    o.add_option('--nside',dest='nside',default=128,help='to define lmax')
    o.add_option('--lmin',dest='lmin',default=50,help='min l for integration. lower will be matter Cl')
    o.add_option('--mCl_file',dest='mCl_file',default='matter_Cl.txt', help='filename of the input matter Cl for extrapolarion')
    o.add_option('--outdir',dest='outdir',help='Path of Output Cl file directory')
    o.add_option('--Dnu',dest='Dnu',default=2.,help='Width in MHz for Bessel integration')
    o.add_option('--bres',dest='bres',default=200,help='Number of bins for Bessel integration. 200 is the lower limit for good behaviour')
    o.add_option('--kres',dest='kres',default=10000,help='resolution for pk integration. Lower values generate noisy results.')
    o.add_option('--log_flag',action='store_true',dest='log_flag',default=False,help='If True the integration is in log')
    o.add_option('--single_cl',action='store_true',dest='single_cl',default=False,help='If True only Cl is computed: the one specified by zmin and zmax')
    o.add_option('--verbose',action='store_true',dest='verbose',default=False)
    

    opts, args = o.parse_args(sys.argv[1:]) 

    print opts, args 

    indir = args[0]  
    nside = int(opts.nside)
    log_flag = opts.log_flag
    lmin=int(opts.lmin)
    lmax=3*nside
    outdir=str(opts.outdir)
    os.mkdir(outdir) #create out directory and check if exits #comment this line if you want to write on the same folder
    single_cl=bool(opts.single_cl)
    zmin,zmax,Dz = np.float(opts.zmin),np.float(opts.zmax),np.float(opts.Dz)

    if (single_cl==False): 
        print "generation z list"
        z_vec=np.arange(zmin,zmax+Dz,Dz)
        nz=len(z_vec)
    else: 
        print "generating single cl"
        z_vec=[zmin,zmax]
        nz=1


    #reading input for integration
    Dnu = float(opts.Dnu)*1e6 #reading Dnu in Hz
    bres = float(opts.bres)
    kres = float(opts.kres)
    verbose = bool(opts.verbose)
        
    print "Reading Pks in" + indir
    #fileformat: pk_z#_z#.txt
    z_list,pk_vec=Ef.readPk_input(indir,z_vec,single_cl)
    npk=len(pk_vec)
    print "Input pks read"
    
    
    print "Computing Cl.."
    ell=np.arange(0,3*nside,1)
    for ipk in range(npk):
        print "Computing cl ", ipk+1 ," of",npk
        k_in=pk_vec[ipk,:,0]
        pk_in=pk_vec[ipk,:,1]
        C_l=np.zeros((len(ell),2),dtype='float')
        C_l[:,0]=ell
        if (z_list[ipk][1]-z_list[ipk][0])<=Dz: 
  
            for l in range(lmin,len(ell)):
                if verbose==True: print('ell= ',l)
                if log_flag==True: 
                    C_l[l,1]=Ef.Pk2Cl_trapz_r_log(l,z_list[ipk][0],z_list[ipk][1],k_in,pk_in,Dnu,bres,kres)
                else: 
                    C_l[l,1]=Ef.Pk2Cl_trapz_r(l,z_list[ipk][0],z_list[ipk][1],k_in,pk_in,Dnu,bres,kres)
    
            print "extrapolating Cl to low l"
            mCl=np.loadtxt(opts.mCl_file) #shape is (3,lmax) 0: ell 1: auto 2:cross
            if z_list[ipk][0]==z_list[ipk][1]: C_l[0:lmin,1]=mCl[1,0:lmin]/mCl[1,lmin]*C_l[lmin,1] #auto 
            else: C_l[0:lmin,1]=mCl[2,0:lmin]/mCl[2,lmin]*C_l[lmin,1] #cross
        
        print "Storing Cl"
        fname=outdir +'Cl_z'+str("%4.3f"% z_list[ipk][0])+'_z'+str("%4.3f"% z_list[ipk][1])+'.txt'
        np.savetxt(fname, C_l)
  
    print "Cl generated."
      
