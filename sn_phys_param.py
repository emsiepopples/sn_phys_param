#!/usr/bin/env python
# sn_phys_param.py by Emma Walker

import numpy as np
import argparse
import __future__

class Supernova():
    
    __doc__ = "Supernova class"
    
    name=''
    distance = float('nan')
    lc_width = float('nan')
    lc_width_err = float('nan')
    rise = float('nan')
    rise_err = float('nan')
    velocity = float('nan')
    velocity_err = float('nan')
    filt = ''
    lc_max = float('nan')
    lc_max_err = float('nan')
    
class Bolometric(Supernova):
    
    __doc__="Bolometric sub-class of Supernova"
    
    luminosity = float('nan')
    luminosity_err = float('nan')
    
    def mag_to_bol(lc_max, lc_max_err):

            #turn bolometric mag to luminosity (ergs/s)
            luminosity = 3.839e33 * 10.0**((4.75 - lc_max)/2.5)
            luminosity_err = luminosity * 0.05
            #DO CALC PROPERLY
            
            return
            
class Monochromatic(Supernova):
    
    __doc__ = "Monochromatic LC sub-class"

    luminosity = float('nan')
    luminosity_err = float('nan')
    
    #ADD CALCULATION


def main():
    
    parser = argparse.ArgumentParser(description='SN Physical Parameter Analogues')
    parser.add_argument('filename', type = str, help = 'Input file')
    parser.add_argument("-f", "--filter", type = str, help = 'Path to filter if not using bol')
    
    args = vars(parser.parse_args())
    
    #first read in the first line of the file
    
    f = open(args['filename'], 'r')
    lines = f.readlines()
    lines = [x for x in lines if not x.startswith('#')] #remove comment lines
    f.close()
    
    #the first line of file MUST be the new supernova
    

if __name__ == "__main__":
    main()