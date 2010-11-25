"""
Generate a synthetic temperature profile and use it in an inversion.
"""

import math

import numpy
import pylab

import geothermics.subcrust as subcrust

# Define the model parameters
conductivity = 3.0          # [W/(m K)]
rad_heat = 1*10**(-7)       # radiogenic head generation [W m^-3]
cond_var = -0.0005          # linear coefficient for conductivity [K^-1] 
crust_temp = 400. + 273.    # temperature at the base of the crust [K]
crust_flux = 20.*10**(-3)   # heat flux at the base of the crust [W m^-2]
crust_depth = 35*10**3      # depth of the base of the crust [m]

# Make a temperature profile
depths = numpy.arange(70000., 200000., 2000.)
temps = subcrust.synthetic_temp_profile(depths, rad_heat, cond_var, 
                                        crust_flux, crust_temp, 
                                        conductivity, crust_depth)

print "Generated %d synthetic data" % (len(temps))

# Contaminate the data with Gaussian noise
error = 15     # Temperature standard deviation [K]
temps_error = numpy.array([numpy.random.normal(temp, error) for temp in temps])

# Set the initial estimate for the inversion
initial_radheat = 10**(-9)
initial_flux = 10**(-6)
initial_condvar = 10**(-4)

# Run the inversion
results = subcrust.invert_temp_profile(depths, temps_error, error, 
                                       initial_radheat, 
                                       initial_condvar, initial_flux, 
                                       ref_temp=crust_temp, 
                                       ref_cond=conductivity, 
                                       ref_depth=crust_depth,
                                       max_it=10000)

inv_radheat, inv_condvar, inv_flux, cov, adjusted, goals = results

sigma_radheat, sigma_condvar, sigma_flux = numpy.sqrt(numpy.diag(cov))

residuals = temps_error - adjusted

print "Inversion results:"
print "   A = %.5g +- %.2g (true=%g)" % (inv_radheat, sigma_radheat, rad_heat)
print "   B = %.5g +- %.2g (true=%g)" % (inv_condvar, sigma_condvar, cond_var)
print "  qc = %.5g +- %.2g (true=%g)" % (inv_flux, sigma_flux, crust_flux)

# Save the results
#output = open("res_per_error.txt", 'a')
#output.write("%g %g %g %g %g %g %g %g %g %g\n" % (error, inv_radheat, 
#    sigma_radheat, rad_heat, inv_condvar, cond_var, sigma_condvar, inv_flux, 
#    crust_flux, sigma_flux))
#output.close()

# Plot the results
pylab.figure()
pylab.plot(temps - 273., depths*0.001, '--b', linewidth=1, label="Noise free")
pylab.ylim(210, 0.001*crust_depth) 
pylab.plot(temps_error - 273., depths*0.001, '.k', 
           label=r"Noisy ($\sigma=%g$ $^\circ$C)" % (error))
pylab.xlabel(r"Temperature [$^\circ$C]")
pylab.ylabel("Depth [km]")
pylab.grid()
pylab.legend(loc='upper right', prop={'size':14}, shadow=True)


pylab.figure(figsize=(16,6))

pylab.subplot(1,2,1)
pylab.plot(temps - 273., depths*0.001, '--b', linewidth=2, label="Noise free") 
pylab.plot(temps_error - 273., depths*0.001, '.k', 
           label=r"Noisy ($\sigma=%g$ $^\circ$C)" % (error))
pylab.plot(adjusted - 273., depths*0.001, '-r', linewidth=2, label="Adjusted") 
pylab.ylim(210, 0.001*crust_depth)
pylab.xlabel(r"Temperature [$^\circ$C]")
pylab.ylabel("Depth [km]")
pylab.grid()
pylab.legend(loc='upper right', prop={'size':14}, shadow=True)

pylab.subplot(1,2,2)
plot = pylab.hist(residuals, bins=len(temps)/10, facecolor='gray')
pylab.xlabel(r"Residuals [$^\circ$C]")
pylab.ylabel("Number of occurrences")

pylab.show()