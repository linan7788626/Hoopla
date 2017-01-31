import numpy as np
import pylab as pl

bsz = 9.0  # arcsec
nnn = 300
dsx = bsz / nnn
xx01 = np.linspace(-bsz / 2.0, bsz / 2.0 - dsx, nnn) + 0.5 * dsx
xx02 = np.linspace(-bsz / 2.0, bsz / 2.0 - dsx, nnn) + 0.5 * dsx
xi1, xi2 = np.meshgrid(xx01, xx02)
kappa = np.loadtxt("./logkap_map.dat")

pl.xlabel("arcsec")
pl.ylabel("arcsec")
pl.contourf(xi1, xi2, kappa)
pl.colorbar()
pl.show()
