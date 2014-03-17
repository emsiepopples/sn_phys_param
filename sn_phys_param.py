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
    
    
    

if __name__ == "__main__":
    main()