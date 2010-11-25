import pylab
import numpy

def lithothick(x, eta, L0):

    return L0*(1 - (1./(1 + eta*numpy.sqrt(x))))


L0 = 90.
eta = 0.5
xmax = 500
xminplot = -20
lw = 2

# Make a continuous plot of the thickness of the lithosphere
xcont = numpy.arange(0, xmax, 0.1)
Lcont = lithothick(xcont, eta, L0)

pylab.figure(figsize=(10,5))
pylab.plot(xcont, Lcont, '-k', linewidth=lw)
pylab.plot([xminplot, xcont[-1]], [L0, L0], '--k', linewidth=lw)
pylab.xlabel('X')
pylab.ylabel('Z')
pylab.xlim(xminplot, xcont[-1])
pylab.ylim(1.2*L0, 0)

# Make a discrete plot of the thickness of the lithosphere
dx = 20
xdisc = numpy.arange(0.5*dx, xmax, dx)
Ldisc = lithothick(xdisc, eta, L0)

Lplot = [0]
xplot = [0]
for x, L in zip(xdisc, Ldisc):

    Lplot.append(L)
    Lplot.append(L)
    Lplot.append(0)
    xplot.append(x - 0.5*dx)
    xplot.append(x + 0.5*dx)
    xplot.append(x + 0.5*dx)

pylab.figure(figsize=(10,5))
pylab.plot(xplot, Lplot, '-k', linewidth=lw)
pylab.plot(xcont, Lcont, '--b', linewidth=lw)
pylab.plot([xminplot, xcont[-1]], [L0, L0], '--k', linewidth=lw)
pylab.xlabel('X')
pylab.ylabel('Z')
pylab.xlim(xminplot, xcont[-1])
pylab.ylim(1.2*L0, 0)

# Plot a regular grid of points
x = range(0, 10)
X, Y = pylab.meshgrid(x, x)

pylab.figure()
pylab.plot(X, Y, '.k')
pylab.xlim(-1, 11)
pylab.ylim(-1, 11)

pylab.show()