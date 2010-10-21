"""
Make a model of sublithospheric crust temperature profile considering a constant
thermal conductivity.
"""

import pylab
import numpy

import geothermics.subcrust as subcrust

# Define the model parameters
conductivity = 3.0          # [W/(m K)]
rad_heat = 1*10**(-7)       # radiogenic head generation [W m^-3]
cond_var = 0.000           # linear coefficient for conductivity [K^-1] 
crust_temp = 400. + 273.    # temperature at the base of the crust [K]
crust_flux = 17.*10**(-3)   # heat flux at the base of the crust [W m^-2]
crust_depth = 35*10**3      # depth of the base of the crust [m]

# Make a temperature profile
depths = numpy.arange(35000., 200000., 1000.)
temps = subcrust.synthetic_temp_profile(depths, rad_heat, cond_var, crust_flux, 
                                        crust_temp, conductivity, crust_depth)

#pylab.savetxt("const_cond_data.txt", numpy.array([depths, temps]).T)

pylab.figure()
pylab.plot(temps - 273., depths*0.001, '-k', linewidth=3, 
           label=r"A=%.1f $\mu W.m^{-3}$" % (10**6*rad_heat))
pylab.ylim(0.001*depths.max(), 0.001*depths.min())
pylab.xlabel(r"Temperature [$^\circ$C]")
pylab.ylabel("Depth [km]")
pylab.grid()
pylab.legend(loc='upper right', prop={'size':14}, shadow=True)
pylab.savefig("const_cond.png")
pylab.show()