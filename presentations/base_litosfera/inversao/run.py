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
    
    qc = 15.*10**(-3)
    
    A = 0.1*10**(-6)
    
#    B = 1./(-1./0.00155 + 273)
    B = -0.002
    
    Tc = 455. + 273.
    
    lambc = 340.
    
    zc = 30000.
    
    z = numpy.arange(80000., 200000., 10000.)
    
    u = U(z, qc, lambc, zc, A)
    
    t1 = T1(B, Tc, u)
    
    t2 = T2(B, Tc, u)
    
    t3 = subcrust.synthetic_temp_profile(z, A, B, qc, Tc, lambc, zc)
    
    pylab.figure()
    
#    pylab.plot(t1, z, '.-b')
    pylab.plot(t2, z, '.-r')
    pylab.plot(t3, z, '.-k')
    
    pylab.ylim(z.max(), z.min())
    
    pylab.show()