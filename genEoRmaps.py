
"""
Generate EoR fullsky maps from Cl

"""

import numpy as np
import healpy as hp
import os,sys
import random
from astropy.io import fits
import scipy.optimize 
from scipy import stats
from scipy.integrate import quad
from scipy.interpolate import interp1d
import scipy.special as ss
import scipy.integrate as integrate
import cosmo_func as cosmo
import EoR_func as Ef



  
if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] indir')
    o.set_description(__doc__)
    o.add_option('--seed',dest='seed',default=1001,help='seed for generating maps')
    o.add_option('--zmin',dest='zmin',default=8,help='zmin of the input Cl')
    o.add_option('--zmax',dest='zmax',default=8.500,help='zmax of the input Cl')
    o.add_option('--Dz',dest='Dz',default=0.125,help='Delta z of the input Cl')
    o.add_option('--nside',dest='nside',default=128,help='Resolution of the maps')
    o.add_option('--outpath',dest='outpath',help='Path of Output map fits files')

    opts, args = o.parse_args(sys.argv[1:]) 

    print opts, args 

    indir = args[0]  
    nside = int(opts.nside)
    npix=hp.nside2npix(nside)
    lmax=3*nside
    seed=int(opts.seed)
    outpath=str(opts.outpath)
    os.mkdir(outpath) #create out directory and check if exits
    zmin,zmax,Dz = np.float(opts.zmin),np.float(opts.zmax),np.float(opts.Dz)
    z_vec=np.arange(zmin,zmax+Dz,Dz)
    nz=len(z_vec)
    np.savez(open(outpath+'input_info.npz','wb'),z=z_vec,nside=nside,seed=seed)
 
        
    print "Reading Cls in" + indir
    #fileformat: Cl_z#_z#.txt
    cl_matrix=Ef.readCl_input(indir,z_vec,lmax)
    print "Input read"
    
    
    print "Generating full sky maps.."

    SKY_cube=Ef.genSKYcube_corr(int(nz), nside,seed, cl_matrix[:,:,:,1])
    
    print "Storing maps"

   
    for zi in range(nz):
        
        fname=outpath+'map_nside'+str(nside)+'_seed'+str(seed)+'_z'+str(zi)+'.fits'
        hp.fitsfunc.write_map(fname, SKY_cube[zi])
  
    print "Eor full sky maps generated."
      
