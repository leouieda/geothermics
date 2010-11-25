"""
Invert the temperature profile of Russell et al. (2001)
"""

import math

import numpy
import pylab

import geothermics.subcrust as subcrust

# Define the model parameters
conductivity = 3.0          # [W/(m K)]
crust_temp = 461 + 273.    # temperature at the base of the crust [K]
crust_depth = 35*10**3      # depth of the base of the crust [m]

depths, temps = pylab.loadtxt('dados.txt', unpack=True)

temps += 273.
depths *= 1000

print "Read %d data" % (len(temps))

error = 15     # Temperature standard deviation [K]

# Set the initial estimate for the inversion
initial_radheat = 10**(-9)
initial_flux = 10**(-3)
initial_condvar = 10**(-3)

# Run the inversion
results = subcrust.invert_temp_profile(depths, temps, error, 
                                       initial_radheat, 
                                       initial_condvar, initial_flux, 
                                       ref_temp=crust_temp, 
                                       ref_cond=conductivity, 
                                       ref_depth=crust_depth,
                                       max_it=10000)

inv_radheat, inv_condvar, inv_flux, cov, adjusted, goals = results

sigma_radheat, sigma_condvar, sigma_flux = numpy.sqrt(numpy.diag(cov))

adjusted = subcrust.synthetic_temp_profile(depths, inv_radheat, inv_condvar, 
                                           inv_flux, crust_temp, conductivity, 
                                           crust_depth)

residuals = temps - adjusted

print "Inversion results:"
print "   A = %.5g +- %.2g" % (inv_radheat, sigma_radheat)
print "   B = %.5g +- %.2g" % (inv_condvar, sigma_condvar)
print "  qc = %.5g +- %.2g" % (inv_flux, sigma_flux)

# Plot the results
pylab.figure()
pylab.plot(temps - 273., depths*0.001, '.k', label="Data (Russell et al, 2001)")
pylab.ylim(210, 0.001*crust_depth) 
pylab.xlabel(r"Temperature [$^\circ$C]")
pylab.ylabel("Depth [km]")
pylab.grid()
pylab.legend(loc='upper right', prop={'size':14}, shadow=True)


pylab.figure(figsize=(16,6))

pylab.subplot(1,2,1)
pylab.plot(temps - 273., depths*0.001, '.k', label=r"Data")
pylab.plot(adjusted - 273., depths*0.001, '.-r', linewidth=1, label="Adjusted") 
pylab.ylim(210, 0.001*crust_depth)
pylab.xlabel(r"Temperature [$^\circ$C]")
pylab.ylabel("Depth [km]")
pylab.grid()
pylab.legend(loc='upper right', prop={'size':14}, shadow=True)

pylab.subplot(1,2,2)
plot = pylab.hist(residuals, bins=len(temps)/6, facecolor='gray')
pylab.xlabel(r"Residuals [$^\circ$C]")
pylab.ylabel("Number of occurrences")

pylab.show()