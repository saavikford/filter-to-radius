#! /usr/bin/env python
#filter_to_AGNradius_v00.py

#######
#
#Program purpose:
# Let's be clever about AGN disk radii and filters
#  -add up continuum as a function of wavelength
#  -find total continuum flux in specified filter
#  -find contribution to total continuum flux in specified filter from annulus of radius R as a function of R
#
#To use:
# User serviceable parts limited to:
#   --variables listed (and documented) between 'BEGIN PHYSICAL INPUTS' and 'END PHYSICAL INPUTS'
#   --variables preceded by 'USER' in the DISPLAY INPUTS section
#   --output file names (preceded by USER) near end of program
#
# To find planned future improvements, search '!!!'
#
# No warranty express or implied. Use at your own risk. Does not prevent drowning.
#
#Descended in part from mcd_gap_v4_0_b1.py (by K. E. Saavik Ford)
#
####
#
# Author: K. E. Saavik Ford
#

from pylab import *

#import sys, os, time, string, math, commands, subprocess
import sys, os, time, string, math, subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import AxesGrid

#SI units
Msun=1.99e30 #kg per solar mass
Rsun=6.95e8 #meters per solar radius
G=6.67e-11
c=3e8
sigma_SB=5.7e-8 #stefan-boltzmann const
yr=3.15e7 #seconds per year
pc=3.086e16 #meters per parsec
AU=1.496e11 #meters per AU
h=6.626e-34 #planck const
kB=1.38e-23 #boltzmann const
m_p=1.67e-27 #mass of proton
sigma_T=6.65e-29 #Thomson xsec
PI=3.1415926

def find_B_lambda(Temp, lam):
    #Planck function
    #BB intensity=2hc^2/lam^5 * 1/(exp(hc/lamkT)-1)
    I=(2*h*c**2/pow(lam,5))/(exp(h*c/(lam*kB*Temp))-1)

    return I

def find_Temp(epsilon,M_SMBH,dotm_edd,radius,r_in,model):
    #find temp as a fn of radius
    #
    #Options:
    #Sirko & Goodman 2003 (SG): T_SG
    #Zero torque at ISCO (ZT): T_ZT  *****DEFAULT
    #Non-zero torque at ISCO (NZT): T_NZT
    #
    #Sirko & Goodman 2003 is similar to non-zero torque at the ISCO
    #---but differs by factor to approximate spectral hardening
    #assume sigma_SB T^4=(3/8pi) dotM Omega^2
    #use Omega=sqrt(GM/r^3); put r in units r_g; dotM in units dotM_edd
    prefactor=pow(((3.0/2.0)*(c**5)*m_p/(epsilon*G*M_SMBH*Msun*sigma_T*sigma_SB)), 0.25)
    T_SG=prefactor*pow(dotm_edd, 0.25)*pow(radius,-0.75)
    
    #Non-zero-torque at ISCO temp profile:
    #include factor f, O(1), which approximates spectral hardening from pure BB
    #set f=2
    #prefactor is otherwise identical to SG
    #See e.g. Gammie 1999; Krolik & Agol 2000; Narayan et al. 1997; Ashfordi & Paczynski 2003
    f=2.0
    T_NZT=f*T_SG
    
    #Zero torque at the ISCO temp profile:
    #see e.g. Zimmerman et al. 2005
    #To implement:
    if (radius >= r_in):
        T_ZT=T_NZT*pow((1.0-pow(radius/r_in, -0.5)),0.25)
    else:
        T_ZT=0.0

    if (model=='SG'):
        T=T_SG
    elif (model=='ZT'):
        T=T_ZT
    elif (model=='NZT'):
        T=T_NZT
    else:
        print 'WARNING: temperature profile model not recognized, using zero-torque model'
        T=T_ZT
        
    return T

def find_area(r1,r2,r_g):
    #find area of annulus
    #2pi R deltaR
    area=2*PI*r1*r_g*(r2*r_g-r1*r_g)

    return area

def get_filter(path, filename):
    #assumes file contains ONLY wavelength and transmission frac, in appropriate units, no headers, etc.
    
    #open the file for reading
    file1 = open(path+filename, 'r')
    #".\data.csv"
    #read in the data line by line, split into float lists of wavelength and transmission fraction
    lam1list=[]
    trans1list=[]
    for line in file1:
        line=line.strip()
        columns=line.split()
        lam1list.append(float(columns[0]))
        trans1list.append(float(columns[1]))

    #close file
    file1.close()

    #re-cast as arrays (from lists) for manipulation
    lam1 = np.array(lam1list)
    trans1 = np.array(trans1list)

    return lam1, trans1

if __name__ == "__main__":
    #F_lam=sum over all blackbodies, assuming T(r), annuli of area=2piRdeltaR
    #integrating over solid angle=pi
    #BB intensity=2hc^2/lam^5 * 1/(exp(hc/lamkT)-1)
    #log_lam, lam=wavelength range of interest
    
    #BEGIN PHYSICAL INPUTS:
    #M_SMBH=mass of supermassive black hole in units of solar masses
    M_SMBH=1.0e7
    #dotm_edd=accretion rate in units of Eddington accretion
    dotm_edd=0.01
    #redshift of object
    redshift=0.0
    #Choose temperature profile model
    #see also 'find_Temp' for more details
    #options are:
    #   SG=Sirko & Goodman 2003
    #   ZT=zero-torque at ISCO model
    #   NZT=non-zero-torque at ISCO model
    #   If none of the above, default is ZT with warning
    tempmod='ZT'

    #choose observing filter (specify filename and path)
    #files should contain ONLY wavelengths in Angstroms (column 0) and transmission fraction (column 1)
    filterpath=''
    filter_filenames=['Palomar_ZTF.g.dat', 'Palomar_ZTF.r.dat']
    
    #USER can, but probably should not, alter the following physical
    #variables without thinking really hard:
    #
    #inner and outer disk radii in units of r_g of SMBH, unperturbed disk
    #(inner radius depends on spin, connect to epsilon later)!!!
    radius_in=6.0
    radius_out=1.0e4
    #epsilon=radiative efficiency, depends on spin of SMBH, assume 0.1
    epsilon=0.1
    #END PHYSICAL INPUTS

    #BEGIN DISPLAY INPUTS:
    #divide display box for graphs
    #params for axes
    #format=left, bottom, width, height
    rect1=0.1,0.1,0.75,0.75
    
    #make figure
    fig1=plt.figure(1)
    #add axes & label them
    ax1=fig1.add_axes(rect1)
    ax1.set_ylabel(r"$F(R_{disk})/F_{tot}$")
    ax1.set_xlabel(r"$R_{disk} \ (R_{g})$")
    #Title include variable values 
    ax1.set_title(r"$M_{SMBH}=%.1e \ M_{\odot}, \ \dotM=%.2f \ \dotM_{edd}, \ z=%.1f$" %(M_SMBH,dotm_edd,redshift))

    #set up ranges
    #Get filter function
    filterlam=[]
    filterlamSI=[]
    filtertrans=[]
    for i in range(len(filter_filenames)):
        fn=filter_filenames[i]
        profilelam, profiletrans=get_filter(filterpath, fn)
        #convert wavelengths from Angstroms to nm
        profilelam=profilelam/10.0
        filterlam.append(profilelam)
        profilelam=profilelam*1.0e-9 #converting nm to m
        filterlamSI.append(profilelam)
        filtertrans.append(profiletrans)

    #USER: choose range for y-axis
    plt.ylim(0.0,0.01)

    #END DISPLAY INPUTS

    #housekeeping:
    #compute r_g for SMBH:
    r_g_SMBH=G*M_SMBH*Msun/c**2
    #set up x-axis variables
    #divvy up disk radii
    log_radius=np.arange(log10(radius_in),log10(radius_out),0.01)
    radius=pow(10,log_radius)
    #initialize final arrays/vars
    F_tot=np.zeros(len(filter_filenames))
    F_tot_ann_iter=[]
    #for each filter
    for j in range(len(filter_filenames)):
        F_lam_tot=np.zeros(len(filterlamSI[j]))
        F_tot_ann=np.zeros(len(radius))
        #compute temp, area, emitted spectrum at each radius; then sum
        for i in range((len(radius)-1)):
            #find disk surface temperature at radius 
            Temp=find_Temp(epsilon,M_SMBH,dotm_edd,radius[i],radius_in,tempmod)
            #find blackbody intensity as a fn of wavelength for all wavelengths in filter, at temp
            B_lambda=find_B_lambda(Temp, filterlamSI[j])
            #compute area of emission
            Area=find_area(radius[i],radius[i+1],r_g_SMBH)
            #flux is pi*A*B_lambda per annulus, pi from integrating over solid angle
            F_lam_ann=PI*Area*B_lambda
            #multiply flux from annulus by filter transmission curve
            F_filter_ann=F_lam_ann*filtertrans[j]
            #sum over all wavelengths for total flux (in filter) coming from radius[i]
            F_tot_ann[i]=sum(F_filter_ann)
            #add annulus to total contributions to filter from this disk
            F_tot[j]=F_tot[j]+F_tot_ann[i]
        #keep track of radial dependence by filter
        F_tot_ann_iter.append(F_tot_ann)

    plt.xlim(0.0, 2000.0)
    for i in range(len(filter_filenames)):
        #plot
        ax1.plot(radius, F_tot_ann_iter[i]/F_tot[i], ls='solid', linewidth=2)

    savefig('filter_to_AGNradius.pdf')
    show()
