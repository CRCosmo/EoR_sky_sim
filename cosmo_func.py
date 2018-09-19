#! /usr/bin/env python

import numpy as np
import astropy.cosmology as CS

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1420405751.77 # FREQUENCY OF 21 CM HYDROGEN LINE
COSMO = CS.FlatLambdaCDM(H0=70.0, Om0=0.3)  # Using H0 = 100 km/s/Mpc

def fq2z(fq):
   '''
   ****************************************************************************************************** 
   Redshift corresponding the specified frequency

   Input(s)
      fq :  [scalar] frequency in Hz
   ******************************************************************************************************
   '''
   z = F21/fq-1
   return z

def z2fq(z):
   '''
   ****************************************************************************************************** 
   Frequency (in Hz) corresponding the specified redshift

   Input(s)
      z :  [scalar] redshift 
   ******************************************************************************************************
   '''
   fq = F21/(z+1)
   return fq

   
def transverse_comoving_distance(z):
   '''
   ****************************************************************************************************** 
   Transverse comoving distance at redshift z corresponding to an angular separation of 1 radian in Mpc/h 

   Input(s)
      z :  [scalar] redshift
   ******************************************************************************************************
   '''
   Dz =  COSMO.comoving_distance(z).value # Mpc/h 
   return Dz



