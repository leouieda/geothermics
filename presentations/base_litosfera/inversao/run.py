"""
Invert a temperature profile for radiogenic heat generation, thermal 
conductivity variation coefficient and heat flux.
"""


import pylab
import numpy

import geothermics.subcrust as subcrust


def U(z, qc, lambc, zc, A):
    
    return qc*(z-zc)/lambc - (0.5*A*(z-zc)**2)/lambc


def T1(B, Tc, U):
    
    return ((B*Tc - 1) + numpy.sqrt(1 + 2*B*U))/B

def T2(B, Tc, U):
    
    return ((B*Tc - 1) - numpy.sqrt(1 + 2*B*U))/B


if __name__ == '__main__':
    
    qc = 13.*10**(-3)
    
    A = 0.0000000001
    
    B = 0.0000
    
    Tc = 400. + 273.
    
    lambc = 3.4
    
    zc = 35000.
    
    z = numpy.arange(35000., 250000., 10000.)
    
    u = U(z, qc, lambc, zc, A)
    print u
    t1 = T1(B, Tc, u)
    print t1 - 273
    
    t2 = T2(B, Tc, u)
    print t2 - 273
    
    t3 = subcrust.synthetic_temp_profile(z, A, B, qc, Tc, lambc, zc)
    print t3 - 273
    z *= 0.001
    
    pylab.figure()
    
    pylab.plot(t1 - 273, z, '.-b')
#    pylab.plot(t2 - 273, z, '.-r')
    pylab.plot(t3 - 273, z, '.-k')
    
    pylab.ylim(z.max(), z.min())
    
    pylab.show()