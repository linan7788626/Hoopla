import sie_lensing_simulation as sls
import numpy as np
import pylab as pl

if __name__ == '__main__':
    num_imgs = 1
    sourcpos = 0.0

    ysc1 = [0.4]
    ysc2 = [-0.3]
    zl = 0.298  # zl is the redshift of the lens galaxy.
    zs = 1.0
    vd = [320]  # Velocity Dispersion.
    q = [0.5]
    pha = [-45.0]

    i = 0
    xf1, xf2, al1, al2, mua_in = sls.cc_for_test(i, ysc1[i], ysc2[i], q[i], vd[i], pha[i], zl, zs)

    yf1 = xf1 - al1
    yf2 = xf2 - al2

    #levels = [-1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0]

    #kappa_out = np.loadtxt("./data_from_Quinn/bestfit_q0.5.logkappa")

    #pl.figure(figsize=(8, 8))
    #pl.contour(xf1, xf2, np.log10(kappa_in), levels, colors=('k',))
    #pl.contour(xf1, xf2, kappa_out, levels, colors=('r',))

    mua_out = np.loadtxt("./data_from_Quinn/nan_subhalo.mag")

    pl.figure(figsize=(8, 8))
    pl.contour(xf1, xf2, mua_in, colors=('k',))
    pl.contour(xf1, xf2, mua_out, colors=('r',))

    pl.figure(figsize=(8, 8))
    pl.contour(yf1, yf2, mua_in, colors=('k',))
    pl.contour(yf1, yf2, mua_out, colors=('g',))

    #vpix = np.loadtxt("./data_from_Quinn/src_pixel.dat")

    pl.show()
