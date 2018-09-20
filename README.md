# EoR_sky_sim
Full-sky simulations of the 21cm EoR signal
CRC-UWC September 2018

Read the pdf for more details.

Python 2 code. Requires: numpy, astropy, scipy and healpy.

Input EoR simulation created from SimFast21 simulations (Santos et al 2005)
The box is 500 Mpc each side. To extrapolate to lower k_min use:

python extrapol_kmin.py indir --zmin=8 --zmax=8.5 --Dz=0.125 --outdir=outdir/ 


First run:
python pk2cl_integration.py indir/ --zmin=8 --zmax=8.5 --Dz=0.125 --nside=128 --lmin=50 --outdir=outdir/ --Dnu=2 --bres=200 --kres=15000 --log_flag --verbose
 


Then run the map generation using as indir the folder created from the previous part of the code:

python genEoRmaps.py indir --zmin=8 --zmax=8.5 --Dz=0.125 --nside=128 --outpath=outpath/
