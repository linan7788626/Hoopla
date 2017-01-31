#!/usr/bin/env python
import numpy as np
import libv4_cv as lv4
import mycosmology as mm
import astropy.io.fits as pyfits
from astropy.cosmology import Planck13
import scipy.interpolate as sci
import pot_ext_shears_kappa as psk
import pylab as pl


def a_b_bh(b, bh):
    res = np.sqrt(b * bh)
    return res


gain = 4.7
expsdss = 53.9
aa_sdss = -24.149
aa_cos = 25.523
kk = 0.156347
airmass = 1.201824


def mag2sdssccd(image):
    im_ccd = gain * expsdss * (10.0**((image + aa_sdss) / (-2.5)))
    return im_ccd


def cosccd2mag(image):
    im_mag = -2.5 * np.log10(image) + aa_cos
    return im_mag


def rebin_psf(input_psf, new_shape):
    nxo, nyo = np.shape(input_psf)
    nxn, nyn = new_shape

    xo = np.linspace(0, nxo - 1.0, nxo) + 0.5
    yo = np.linspace(0, nyo - 1.0, nyo) + 0.5
    xo, yo = np.meshgrid(xo, yo)
    xo = xo.reshape((nxo * nyo))
    yo = yo.reshape((nxo * nyo))
    zo = input_psf.reshape((nxo * nyo))

    xn = np.linspace(0, nxo - 1.0, nxn) + 0.5
    yn = np.linspace(0, nyo - 1.0, nyn) + 0.5
    xn, yn = np.meshgrid(xn, yn)

    res = sci.griddata(np.array([xo, yo]).T, zo, (xn, yn), method='linear')
    return res


nMgyCount_r = 0.004760406   # nanomaggies per count for SDSS detector.
sky_r = 5.98          # SDSS typical r band sky
softbias = 1000.0        # SDSS softbias
Mgy2nanoMgy = 10e+9         # nanoMaggy to Maggy
aa_r = -24.149
kk = 0.156347
skycount = sky_r / (nMgyCount_r)
expsdss = 53.9
gain = 4.7
airmass = 1.201824
factor = 10.0**(0.4 * (aa_r + kk * airmass))


def psf_gaussian_norm(x1, x2, mu, sigma):
    r = np.sqrt(x1 * x1 + x2 * x2)
    res = 1.0 / (sigma * np.sqrt(2.0 * np.pi)) * \
        np.exp(-(r - mu)**2.0 / 2.0 * sigma**2.0)
    return res


def noise_map(nx1, nx2, nstd, NoiseType):
    if NoiseType == 'Poisson':
        noise = np.random.poisson(nstd, (nx1, nx2)) - nstd
    if NoiseType == 'Gaussian':
        noise = nstd * np.random.normal(0.0, 1.0, (nx1, nx2))
    return noise


def make_r_coor(nc, dsx):

    bsz = nc * dsx
    x1 = np.linspace(0, bsz - dsx, nc) - bsz / 2.0 + dsx / 2.0
    x2 = np.linspace(0, bsz - dsx, nc) - bsz / 2.0 + dsx / 2.0

    x2, x1 = np.meshgrid(x1, x2)
    return x1, x2


def make_c_coor(nc, dsx):

    bsz = nc * dsx
    x1, x2 = np.mgrid[0:(bsz - dsx):nc * 1j, 0:(bsz - dsx):nc * 1j] \
        - bsz / 2.0 + dsx / 2.0
    return x1, x2


def lens_equation_sie(x1, x2, lpar):
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = lpar[0]
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]
    q = lpar[2]  # Ellipticity of lens.
    rc = lpar[3]  # Core size of lens (in units of Einstein radius).
    re = lpar[4]  # Einstein radius of lens.
    pha = lpar[5]  # Orintation of lens.

    phirad = np.deg2rad(pha)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1 - xc1) * cosa + (x2 - xc2) * sina
    xt2 = (x2 - xc2) * cosa - (x1 - xc1) * sina

    phi = np.sqrt(xt2 * xt2 + xt1 * q * xt1 * q + rc * rc)
    sq = np.sqrt(1.0 - q * q)
    pd1 = phi + rc / q
    pd2 = phi + rc * q
    fx1 = sq * xt1 / pd1
    fx2 = sq * xt2 / pd2
    qs = np.sqrt(q)

    a1 = qs / sq * np.arctan(fx1)
    a2 = qs / sq * np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 = -sina

    fx11 = xt11 / pd1 - xt1 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd1 * pd1)
    fx22 = xt22 / pd2 - xt2 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd2 * pd2)
    fx12 = xt12 / pd1 - xt1 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd1 * pd1)
    fx21 = xt21 / pd2 - xt2 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd2 * pd2)

    a11 = qs / (1.0 + fx1 * fx1) * fx11
    a22 = qs / (1.0 - fx2 * fx2) * fx22
    a12 = qs / (1.0 + fx1 * fx1) * fx12
    a21 = qs / (1.0 - fx2 * fx2) * fx21

    rea11 = (a11 * cosa - a21 * sina) * re
    rea22 = (a22 * cosa + a12 * sina) * re
    rea12 = (a12 * cosa - a22 * sina) * re
    rea21 = (a21 * cosa + a11 * sina) * re

    y11 = 1.0 - rea11
    y22 = 1.0 - rea22
    y12 = 0.0 - rea12
    y21 = 0.0 - rea21

    jacobian = y11 * y22 - y12 * y21
    mu = 1.0 / jacobian

    res1 = (a1 * cosa - a2 * sina) * re
    res2 = (a2 * cosa + a1 * sina) * re
    return res1, res2, mu


def lensing_signals_sie(x1, x2, lpar):
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = lpar[0]
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]
    q = lpar[2]   # Axis ratio of lens.
    rc = lpar[3]   # Core size of lens (in units of Einstein radius).
    re = lpar[4]   # Einstein radius of lens.
    pha = lpar[5]   # Orintation of lens.

    phirad = np.deg2rad(pha)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1 - xc1) * cosa + (x2 - xc2) * sina
    xt2 = (x2 - xc2) * cosa - (x1 - xc1) * sina

    phi = np.sqrt(xt2 * xt2 + xt1 * q * xt1 * q + rc * rc)
    sq = np.sqrt(1.0 - q * q)
    pd1 = phi + rc / q
    pd2 = phi + rc * q
    fx1 = sq * xt1 / pd1
    fx2 = sq * xt2 / pd2
    qs = np.sqrt(q)

    a1 = qs / sq * np.arctan(fx1)
    a2 = qs / sq * np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 = -sina

    fx11 = xt11 / pd1 - xt1 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd1 * pd1)
    fx22 = xt22 / pd2 - xt2 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd2 * pd2)
    fx12 = xt12 / pd1 - xt1 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd1 * pd1)
    fx21 = xt21 / pd2 - xt2 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd2 * pd2)

    a11 = qs / (1.0 + fx1 * fx1) * fx11
    a22 = qs / (1.0 - fx2 * fx2) * fx22
    a12 = qs / (1.0 + fx1 * fx1) * fx12
    a21 = qs / (1.0 - fx2 * fx2) * fx21

    rea11 = (a11 * cosa - a21 * sina) * re
    rea22 = (a22 * cosa + a12 * sina) * re
    rea12 = (a12 * cosa - a22 * sina) * re
    rea21 = (a21 * cosa + a11 * sina) * re

    kappa = 0.5 * (rea11 + rea22)
    shear1 = 0.5 * (rea12 + rea21)
    shear2 = 0.5 * (rea11 - rea22)

    y11 = 1.0 - rea11
    y22 = 1.0 - rea22
    y12 = 0.0 - rea12
    y21 = 0.0 - rea21

    jacobian = y11 * y22 - y12 * y21
    mu = 1.0 / jacobian

    alpha1 = (a1 * cosa - a2 * sina) * re
    alpha2 = (a2 * cosa + a1 * sina) * re

    return alpha1, alpha2, kappa, shear1, shear2, mu


def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
    return (xnew, ynew)


def gauss_2d(x, y, par):
    (xnew, ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = np.sqrt(((xnew**2) * par[4] + (ynew**2) / par[4])) / np.abs(par[1])
    res = par[0] * np.exp(-res0**2.0)
    return res


def re_sv(sv, z1, z2):
    res = 4.0 * np.pi * (sv**2.0 / mm.vc**2.0) * \
        mm.Da2(z1, z2) / mm.Da(z2) * mm.apr
    return res


def Brightness(Re, Vd):
    a = 1.49
    b = 0.2
    c = -8.778
    mag_e = ((np.log10(Re) - a * np.log10(Vd) - c) / b) + \
        20.09  # Bernardi et al 2003
    nanoMgy = Mgy2nanoMgy * 10.0**(-(mag_e - 22.5) / 2.5)
    counts = nanoMgy / nMgyCount_r

    return counts


def de_vaucouleurs_2d(x, y, par):
    # [I0, Re, xc1,xc2,q,pha]
    # print "I0",par[0]
    # print "Re",par[1]
    (xnew, ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = np.sqrt((xnew**2) * par[4] + (ynew**2) / par[4]) / par[1]
    # res = par[0]*np.exp(-par[1]*res0**0.25)
    res = par[0] * np.exp(-7.669 * (res0**0.25 - 1.0))
    soften = par[0] * np.exp(-7.669 * ((0.2)**0.25 - 1.0))
    res[res > soften] = soften
    return res

def cc_for_test(ind, ysc1, ysc2, q, vd, pha, zl, zs, lens_tag=1):
    # dsx_sdss     = 0.396         # pixel size of SDSS detector.
    R = 3.0000     #
    nnn = 300  # Image dimension
    bsz = 9.0  # arcsecs
    dsx = bsz / nnn         # pixel size of SDSS detector.
    nstd = 59  # ^2

    xi1, xi2 = make_r_coor(nnn, dsx)
# ----------------------------------------------------------------------
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = 0.0
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = 0.0
    rc = 0.0  # Core size of lens (in units of Einstein radius).
    re = re_sv(vd, zl, zs)  # Einstein radius of lens.
    re_sub = 0.05 * re
    a_sub = a_b_bh(re_sub, re)
    ext_shears = 0.1
    ext_angle = 0.0
    ext_kappa = 0.2

# ----------------------------------------------------------------------

    ai1, ai2 = psk.deflection_nie(xc1, xc2, pha, q, re, rc, ext_shears, ext_angle,
                                  ext_kappa, xi1, xi2)

    as1, as2 = psk.deflection_sub_pJaffe(0.0, -2.169, re_sub, 0.0, a_sub, xi1, xi2)

    al1 =  ai1 + as1
    al2 =  ai2 + as2

    al11,al12 = np.gradient(al1,dsx)
    al21,al22 = np.gradient(al2,dsx)

    mua = 1.0/(1.0-(al11+al22)+al11*al22-al12*al21)

    return xi1,xi2, al1, al2, mua


def single_run_test(ind, ysc1, ysc2, q, vd, pha, zl, zs, lens_tag=1):
    # dsx_sdss     = 0.396         # pixel size of SDSS detector.
    R = 3.0000     #
    nnn = 300  # Image dimension
    bsz = 9.0  # arcsecs
    dsx = bsz / nnn         # pixel size of SDSS detector.
    nstd = 59  # ^2

    xi1, xi2 = make_r_coor(nnn, dsx)
# ----------------------------------------------------------------------
    dsi = 0.03
    g_source = pyfits.getdata(
        "./gals_sources/439.0_149.482739_1.889989_processed.fits")
    g_source = np.array(g_source, dtype="<d") * 10.0
    g_source[g_source <= 0.0001] = 1e-6
# ----------------------------------------------------------------------
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = 0.0
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = 0.0
    rc = 0.0  # Core size of lens (in units of Einstein radius).
    re = re_sv(vd, zl, zs)  # Einstein radius of lens.
    re_sub = 0.05 * re
    a_sub = a_b_bh(re_sub, re)
    ext_shears = 0.1
    ext_angle = 0.0
    ext_kappa = 0.2

# ----------------------------------------------------------------------
    #lpar = np.asarray([xc1, xc2, q, rc, re, pha])
    #ai1, ai2, kappa_out, shr1, shr2, mua = lensing_signals_sie(xi1, xi2, lpar)
    #ar = np.sqrt(ai1 * ai1 + ai2 * ai2)

    # psi_nie = psk.potential_nie(xc1, xc2, pha, q, re, rc, ext_shears, ext_angle,
    #                            ext_kappa, xi1, xi2)
    #ai1, ai2 = np.gradient(psi_nie, dsx)

    ai1, ai2 = psk.deflection_nie(xc1, xc2, pha, q, re, rc, ext_shears, ext_angle,
                                  ext_kappa, xi1, xi2)

    as1, as2 = psk.deflection_sub_pJaffe(0.0, -2.169, re_sub, 0.0, a_sub, xi1, xi2)

    yi1 = xi1 - ai1 - as1
    yi2 = xi2 - ai2 - as2

    g_limage = lv4.call_ray_tracing(g_source, yi1, yi2, ysc1, ysc2, dsi)
    g_limage[g_limage <= 0.25] = 1e-6

    pl.figure()
    pl.contourf(g_limage)
    pl.colorbar()

    g_limage = cosccd2mag(g_limage)
    g_limage = mag2sdssccd(g_limage)

    pl.figure()
    pl.contourf(xi1,xi2,g_limage)
    pl.colorbar()
# -------------------------------------------------------------
    dA = Planck13.comoving_distance(zl).value * 1000. / (1.0 + zl)
    Re = dA * np.sin(R * np.pi / 180. / 3600.)
    counts = Brightness(Re, vd)
    vpar = np.asarray([counts, R, xc1, xc2, q, pha])
    g_lens = de_vaucouleurs_2d(xi1, xi2, vpar)

    g_clean_ccd = g_lens * lens_tag + g_limage
    output_filename = "./fits_outputs/clean_lensed_imgs.fits"
    #pyfits.writeto(output_filename, g_clean_ccd, clobber=True)
# -------------------------------------------------------------
    from scipy.ndimage.filters import gaussian_filter
    g_images_psf = gaussian_filter(g_clean_ccd, 2.0)
# -------------------------------------------------------------
    g_noise = noise_map(nnn, nnn, np.sqrt(nstd), "Gaussian")
    output_filename = "./fits_outputs/noise_map.fits"
    #pyfits.writeto(output_filename, g_noise, clobber=True)
    g_final = g_images_psf + g_noise
# -------------------------------------------------------------
    output_filename = "./fits_outputs/lensed_imgs_only.fits"
    #pyfits.writeto(output_filename, g_final, clobber=True)

    pl.figure()
    pl.contourf(xi1,xi2,g_final)
    pl.colorbar()

    return 0

if __name__ == '__main__':
    # from mpi4py import MPI
    # import sys
    # sourcpos = 10.0 # arcsecs
    # num_imgs = int(sys.argv[1])
    num_imgs = 1
    sourcpos = 0.0

    # comm = MPI.COMM_WORLD
    # size = comm.Get_size()
    # rank = comm.Get_rank()

    # ysc1 = np.random.random(num_imgs)*sourcpos-sourcpos/2.0
    # ysc2 = np.random.random(num_imgs)*sourcpos-sourcpos/2.0
    # q = np.random.random(num_imgs)*0.5+0.5
    # vd = np.random.random(num_imgs)*100.0+200.0
    # pha = np.random.random(num_imgs)*360.0
    # zl = 0.2
    # zs = 1.0

    ysc1 = [0.4]
    ysc2 = [-0.3]
    zl = 0.298    # zl is the redshift of the lens galaxy.
    zs = 1.0
    vd = [320]    # Velocity Dispersion.
    q = [0.5]  # 0.5
    pha = [-45.0]         # -45.0

    # for i in xrange(rank,num_imgs,size):
    for i in xrange(num_imgs):
        single_run_test(i, ysc1[i], ysc2[i], q[i], vd[i], pha[i], zl, zs)
    pl.show()
