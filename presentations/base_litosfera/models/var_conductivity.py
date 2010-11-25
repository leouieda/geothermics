"""
Make a model of sublithospheric crust temperature profile considering a variant
thermal conductivity.
"""

import pylab
import numpy

import geothermics.subcrust as subcrust

# Define the model parameters
conductivity = 3.0          # [W/(m K)]
rad_heat = 1.*10**(-7)       # radiogenic head generation [W m^-3]
cond_var = 0.001            # linear coefficient for conductivity [K^-1] 
crust_temp = 400. + 273.    # temperature at the base of the crust [K]
crust_flux = 17.*10**(-3)   # heat flux at the base of the crust [W m^-2]
crust_depth = 35*10**3      # depth of the base of the crust [m]

# Make temperature profiles
depths = numpy.arange(35000., 200000., 1000.)

pylab.figure(figsize=(11,6))

for cond_var in [-10**(-6)*10**i for i in xrange(4)]:
    
    temps = subcrust.synthetic_temp_profile(depths, rad_heat, cond_var, 
                                            crust_flux, crust_temp, 
                                            conductivity, crust_depth)
    pylab.plot(temps - 273., depths*0.001, '-.', 
               label=r"B=%2.1e $K^{-1}$" % (cond_var), linewidth=3)

temps = subcrust.synthetic_temp_profile(depths, rad_heat, 0, 
                                        crust_flux, crust_temp, 
                                        conductivity, crust_depth)

pylab.plot(temps - 273., depths*0.001, '-k', label=r"B=0 $K^{-1}$", linewidth=3)

for cond_var in [10**(-6)*10**i for i in xrange(5)]:
    
    temps = subcrust.synthetic_temp_profile(depths, rad_heat, cond_var, 
                                            crust_flux, crust_temp, 
                                            conductivity, crust_depth)
    pylab.plot(temps - 273., depths*0.001, '--', 
               label=r"B=%2.1e $K^{-1}$" % (cond_var), linewidth=3)
    
pylab.ylim(0.001*depths.max(), 0.001*depths.min())
pylab.xlabel(r"Temperature [$^\circ$C]")
pylab.ylabel("Depth [km]")
pylab.grid()
pylab.legend(loc='upper right', prop={'size':12}, shadow=True)
pylab.savefig("var_cond.png")
pylab.show()