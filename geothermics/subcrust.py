# Copyright 2010 Leonardo Uieda
#
# This file is part of Geothermics.
#
# Fatiando a Terra is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Geothermics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Geothermics.  If not, see <http://www.gnu.org/licenses/>.
"""
Subcrustal Lithosphere geothermal modeling and inversion.
"""
__author__ = 'Leonardo Uieda <leouieda@gmail.com>'


import time

import numpy


def synthetic_temp_profile(depths, radheat, condvar, ref_flux, ref_temp,
                           ref_cond):
    """
    Generate a synthetic temperature profile using radiogenic heat generation
    and a linearly temperature dependent thermal conductivity.

    For subcrustal lithospheric modeling, use the base of the crust as the
    reference surface.

    Parameters:

    * depths
        List of depths at which the temperature will be calculated

    * radheat
        Rate of radiogenic heat generation

    * condvar
        Rate with which the conductivity varies with temperature (a linear
        coefficient)

    * ref_flux
        The heat flux at the reference surface

    * ref_temp
        Temperature at the reference surface

    * ref_cond
        Thermal conductivity at the reference surface

    Returns:

    * temp_profile
        Temperatures at the given *depths*
    """

    pass



def _param_jacobian(estimate, data, depths, ref_temp, ref_cond, ref_depth):
    """
    Calculate the Jacobian matrix of the mathematical model with respect to the
    parameters.
    """

    qc_deriv = (ref_depth - depths)/ref_cond
    
    A_deriv = (depths - ref_depth)/(2*ref_cond)
    
    B_deriv = 0.5*(data - ref_temp)**2 

    jacobian = numpy.array([A_deriv, B_deriv, qc_deriv]).T
    
    return jacobian


def _data_jacobian(estimate, data, ref_temp):
    """
    Calculate the Jacobian matrix of the mathematical model with respect to the
    data.
    """

    B = estimate[1]
    
    jacobian = numpy.zeros((len(data), len(data)))
    
    for i in xrange(len(data)):
        
        jacobian[i][i] = 1 + B*(data[i] - ref_temp)

    return jacobian


def _model(estimate, data, depths, ref_temp, ref_cond, ref_depth):
    """
    Calculate the mathematical model for a given estimate and data.
    Ideally should be zero.
    """

    A = estimate[0]
    B = estimate[1]
    qc = estimate[2]

    f = (data - ref_temp + 0.5*B*(data - ref_temp)**2 - 
         qc*(depths - ref_depth)/ref_cond + 
         (0.5*A*(depths - ref_depth)**2)/ref_cond)
   
    return f


def invert_temp_profile(depths, temps, initial_radheat, initial_condvar,
                        initial_flux, ref_temp, ref_cond, ref_depth, max_it=50):
    """
    Invert a temperature profile for the radiogenic heat generation, linear
    thermal conductivity variation coefficient, and heat flux at the reference
    surface.

    Parameters:

    * depths
        List of depths at which the temperatures where measured

    * temps
        List of temperatures at the corresponding *depths*

    * initial_radheat
        Initial estimate for the radiogenic heat generation

    * initial_condvar
        Initial estimate for the linear thermal conductivity variation
        coefficient

    * initial_flux
        Initial estimate for the heat flux

    * ref_temp
        Temperature at the reference surface

    * ref_cond
        Conductivity at the reference surface
        
    * max_it
        Maximum iterations

    Returns:

    * list with [radheat, condvar, flux, residuals, goals]
        *radheat*, *condvar* and *flux* are the inversion results.
        *residuals* is a list with the residuals
        *goals* is a list with the goal function value per iteration
    """

    print "Inverting temperature profile:"
    print "  initial radiogenic heat=%g" % (initial_radheat)
    print "  initial conductivity variation=%g" % (initial_condvar)
    print "  initial heat flux=%g" % (initial_flux)
    print "  reference temperature=%g" % (ref_temp)
    print "  reference conductivity=%g" % (ref_cond)
    print "  reference depth=%g" % (ref_depth)
    print "  max iterations=%d" % (max_it)

    goals = []

    next = [initial_radheat, initial_condvar, initial_flux]

    # Need this to calculate the initial residuals
    f0 = _model(next, temps, depths, ref_temp, ref_cond, ref_depth)
    data_jac = _data_jacobian(next, temps, ref_depth)
    param_jac = _param_jacobian(next, temps, depths, ref_temp, ref_cond, 
                                ref_depth)

    # Also need to calculate the Lagrange multiplier
    system = numpy.dot(data_jac, data_jac.T)
    lagrange_mult = numpy.linalg.solve(system, f0)

    residuals = -1*numpy.dot(data_jac.T, lagrange_mult)

    # Since there is no regularization, the goal function is just the rms
    rms = (residuals*residuals).sum()
    goals.append(rms)

    # So that I can warn the user if exited because of max_it and not because of
    # convergence
    max_it_exit = True

    for iteration in xrange(max_it):

        start = time.time()

        prev = next

        f0 = _model(prev, temps, depths, ref_temp, ref_cond, ref_depth)
        data_jac = _data_jacobian(prev, temps, ref_depth)
        param_jac = _param_jacobian(prev, temps, depths, ref_temp, ref_cond, 
                                    ref_depth)

        data_jac_inv = numpy.linalg.inv(numpy.dot(data_jac, data_jac.T))

        aux = numpy.dot(param_jac.T, data_jac_inv)

        # Solve for the correction
        normal_eq_sys = numpy.dot(aux, param_jac)
        y = -1*numpy.dot(aux, f0)
        correction = numpy.linalg.solve(normal_eq_sys, y)

        misfit = numpy.dot(param_jac, correction) + f0

        lagrange_mult = numpy.dot(data_jac_inv, misfit)

        residuals = -1*numpy.dot(data_jac.T, lagrange_mult)

        rms = (residuals*residuals).sum()

        goals.append(rms)
        
        next = prev + correction

        end = time.time()
        
        print "  it %d: RMS=%g (%g s)" % (iteration + 1, rms, end - start)

        if abs((goals[-1] - goals[-2])/goals[-2]) <= 10**(-4):

            max_it_exit = False

            break

    if max_it_exit:

        print "WARNING! Exited due to reaching maximum number of iterations."

    return [next[0], next[1], next[2], residuals, goals]