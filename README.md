sn_phys_param
=============

#### Overview

One piece of code for running SN analogue calculations.

Reading in 1 file, the code will output the ejecta mass, kinetic energy and nickel mass for a supernova based on any number of analogues.

Based on Arnett's equations, plus the equation for the nickel mass in Walker et al (submitted).

#### Options

- -h help
- -f path to filter transmission curve if using monochromatic data
- -z the zero point value for Vega in the filter.  Default = 0
- -v verbose printing of the data for each analogue. Default is to have this turned off and only the weighted mean of each quantity is printed at the end

####Examples

- run sn_phys_param.py 'testfile.dat' -v  
This uses the bolometric option as no filter file is supplied.  Includes the verbose output option where the values of the new SN are shown based on each analogue, rather than just printing the weighted mean at the end.
- run sn_phys_param.py 'testfile2.dat' -f 'bess-r.pass' -z 0.03  
This uses the monochromatic option by supplying an R-band filter. It also specifies a Vega zeropoint = 0.03 for this filter

####Input File Format

Requires an input file with the following format:

**First line**

name, distance, filter, lc_width, lc_width_err, lc_max magnitude, lc_max_err, velocity, vel_err, rise time, rise_err

name: This should be obvious

distance: the distance to the SN in Mpc.  This avoids the need to incorporate any cosmology in the programme: decide before you use it

filter: this should be either 'bol' or the name of the filter you are using

lc_width(_err): a measurement that increases for brighter objects eg the width of the lightcurve at max - 0.5mags or 1/dm15®

lc_max(_err): magnitude at maximum light.  Either bolometric or apparent if using monochromatic data

velocity(_err): the photospheric velocity at maximum

rise(_err): the rise time to maximum

**Following Lines**

In addition to this you can have any number of analogues to detail on the next lines.  The code will loop over them so there is no maximum number.  The format for these lines is:

name distance filter lc_width lc_width_err lc_max lc_max_err, vel, vel_err, rise, rise_err, eject, eject_err, ke, ke_err, ni, ni_err

which is the same as above with the additional columns

eject(_err): the ejecta mass of the supernova

ke(_err): kinetic energy

ni(_err): synthesised 56Ni mass

**Comments**

Comments should be added only with the '#' symbol.

