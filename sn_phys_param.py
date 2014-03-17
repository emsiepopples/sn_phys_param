#!/usr/bin/env python
# sn_phys_param.py by Emma Walker
# assumes:
# ergs /s
# Vega magnitudes

import numpy as np
import argparse
import __future__
from scipy import interpolate
from scipy.integrate import simps

def locate(x, a):
    t = abs(np.subtract(x,a))
    return min(xrange(len(t)), key=t.__getitem__)


def int_flux(ww, ff, tw, tf):
    a1 = locate(ww, min(tw))
    a2 = locate(ww, max(tw))
    wcut = ww[a1:a2]
    fcut = ff[a1:a2]
    t = interpolate.splrep(tw, tf, s=0)
    tfint = interpolate.splev(wcut, t, der=0)
    fluxtrans = np.multiply(fcut, tfint)
    return simps(fluxtrans, wcut)

class Supernova(object):
    
    __doc__ = "Supernova class"
    
    def __init__(self, inarr):
        
        vals = inarr.split()
        self.name = vals[0]
        self.distance = float(vals[1])
        self.filt = vals[2]
        self.lc_width = float(vals[3])
        self.lc_width_err = float(vals[4])
        self.lc_max = float(vals[5)
        self.lc_max_err = float(vals[6])
        self.velocity = float(vals[7])
        self.velocity_err = float(vals[8])
        self.rise = float(vals[9])
        self.rise_err = float(vals[10])
        self.luminosity = float('nan')
        self.luminosity_err = float('nan')
        self.mej = float('nan')
        self.mej_err = float('nan')
        self.ke = float('nan')
        self.ke_err = float('nan')
        self.ni = float('nan')
        self.ni_err = float('nan')
  
class Bolometric(Supernova):
    
    __doc__="Bolometric sub-class of Supernova"
    
    def mag_to_bol(self):

        #turn bolometric mag to luminosity (ergs/s)
        self.luminosity = 3.839e33 * 10.0**((4.75 - self.lc_max)/2.5)
        self.luminosity_err = np.sqrt((self.luminosity * np.log(10.)/2.5)**2 * self.lc_max_err**2) 
        #DO CALC PROPERLY
                 
        return
            
class Monochromatic(Supernova):
    
    __doc__ = "Monochromatic LC sub-class"
    
    def mag_to_lum(self, filterfile, vegazpt):                                                                                                                                                                                                                                                                                             ):
        
        #first need Vega and filter

        vega = np.genfromtxt('vega.dat', unpack = True, dtype = None)
        
        trans = np.genfromtxt(filterfile, unpack = True, dtype = None)
        
        vegaflux = int_flux(vega[0], vega[1], trans[0], trans[1])
        flux = vegaflux * 10.0 **(-(self.lc_max - vegazpt)/-2.5)
        flux_err = np.sqrt( (flux * np.log(10.)/-2.5)**2 * self.lc_max_err**2)
        
        self.luminosity = 4.0 * np.pi * (self.distance * 1e6 * 3.086e18)**2 * flux
        self.luminosity_err = np.sqrt( (4.0 * np.pi * (self.distance * 1e6 * 3.086e18)**2)**2 * flux_err**2)
    
        return
        
def ejecta_mass(new, comp):
    
    mass = new.lc_width**2 * new.velocity / (comp.lc_width**2 * comp.velocity) * comp.mej
    
    e1 = (2 * new.lc_width / comp.lc_width**2 * new.velocity / comp.velocity * comp.mej)**2 * new.lc_width_err**2
	e2 = (-2.0 * new.lc_width**2 / comp.lc_width**3 * new.velocity / comp.velocity)**2 * comp.lc_width_err**2
	e3 = ( new.lc_width**2 / comp.lc_width**2 * new.velocity / comp.velocity)**2 * comp.mej_err**2
	#velocity terms
	e4 = ((new.lc_width/comp.lc_width)**2 * 1./comp.velocity * comp.mej)**2 * new.velocity_err**2
	e5 = ((new.lc_width/comp.lc_width)**2 * (new.velocity/(comp.velocity)**2)*comp.mej)**2 * comp.velocity_err**2
	error = np.sqrt(e1 + e2 + e3 + e4 + e5)
	return (mass, error)
    
def kinetic_energy(new, comp):
    
	ke = (new.lc_width / comp.lc_width)**2 * (new.velocity / comp.velocity)**3 * comp.ke
	e1 = (2. * new.lc_width / comp.lc_width**2 * (new.velocity / comp.velocity)**3 * comp.ke)**2 * new.lc_width_err**2
	e2 = (-2. * new.lc_width**2 / comp.lc_width**3 * (new.velocity / comp.velocity)**3 * comp.ke)**2 * comp.lc_width_err**2
	e3 = ( new.lc_width**2 / comp.lc_width**2 * (new.velocity / comp.velocity)**3) **2 * comp.ke_err**2
	
	#velocity terms
	
	e4 = ((new.lc_width/comp.lc_width)**2 * 3./comp.velocity * (new.velocity/comp.velocity)**2.*comp.ke)**2 * comp.velocity_err**2
	
	e5 = ((new.lc_width/comp.lc_width)**2 * -3.0/new.velocity * (new.velocity/(comp.velocity)**4)*comp.ke)**2 * comp.velocity_err**2
	
	error = np.sqrt(e1 + e2 + e3 + e4 + e5)
	return (ke, error)
    
def nickel_mass(new, comp):
    
	tau_ni = 6.08 / np.log(2)
	tau_co = 77.23 / np.log(2)

	eni = 1.7
	eco = 3.67

	kA = eni / tau_ni
	kB = eco / tau_co * tau_co / (tau_co - tau_ni)

	f1 = kA * np.exp(-new.rise / tau_ni) + kB * (np.exp(-new.rise / tau_co) - np.exp(-new.rise/ tau_ni))

	f2 = kA * np.exp(-comp.rise / tau_ni )+ kB * (np.exp(-comp.rise / tau_co) - np.exp(-comp.rise/ tau_ni))


	mass = comp.ni * new.luminosity / comp.luminosity  * f2 / f1

	e1 = (new.luminosity / comp.luminosity * f1 / f2)**2 * comp_ni_err**2

	e2 = (comp.ni/comp.luminosity * f1 / f2)**2 * new.luminosity_err**2

	e3 = (-comp.ni / comp.luminosity**2 * new.luminosity * f1 / f2)**2 * comp.luminosity_err**2

	e4 = (-comp.ni *new.luminosity / comp.luminosity * f2 / f1**2 *(kA/(-tau_ni) * np.exp(-new.rise/tau_ni) + kB/(-tau_co)*np.exp(-new.rise/tau_co) - kB/(-tau_ni) * np.exp(-new.rise/tau_ni)))**2 * new.rise_err**2

	e5 = (-comp.ni *new.luminosity / comp.luminosity * 1./f1 * (kA/(-tau_ni) * np.exp(-comp.rise/tau_ni) + kB/(-tau_co)*np.exp(-comp.rise/tau_co) - kB/(-tau_ni) * np.exp(-comp.rise/tau_ni)))**2 * comp.rise_err**2

	#e5 = 0 for 06aj and 98bw as we have the XRF/GRB to time the explosion

	err = np.sqrt(e1 + e2 + e3 + e4 + e5)

	return (mass, err)
        



def main():
    
    parser = argparse.ArgumentParser(description='SN Physical Parameter Analogues')
    parser.add_argument('filename', type = str, help = 'Input file')
    parser.add_argument("-f", "--filter", type = str, help = 'Path to filter if not using bol')
    parser.add_argument("-z", "--zeropt", type = float, default = 0.0, help = 'Magntiude of Vega in this band.  Default = 0.0')
    
    args = vars(parser.parse_args())
    
    #first read in the first line of the file
    
    f = open(args['filename'], 'r')
    lines = f.readlines()
    lines = [x for x in lines if not x.startswith('#')] #remove comment lines
    f.close()
    
    #the first line of file MUST be the new supernova
    
    new_sn = Supernova(lines[0])
    
    n_comp = np.size(lines) - 1
    
    final_vals = {'meject':[]}
    

if __name__ == "__main__":
    main()