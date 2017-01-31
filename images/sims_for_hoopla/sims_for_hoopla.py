#!/usr/bin/env python
import numpy as np
import libv4_cv as lv4
import mycosmology as mm
import astropy.io.fits as pyfits
# from astropy.cosmology import Planck13
# import scipy.interpolate as sci
import pot_ext_shears_kappa as psk
import pylab as pl


def a_b_bh(b, bh):
    res = np.sqrt(b * bh)
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


def de_vaucouleurs_2d(x1, x2, xc1, xc2, afactor, Reff, ell, ori):
    (xnew, ynew) = xy_rotate(x1, x2, xc1, xc2, ori)
    res0 = np.sqrt((xnew**2) * ell + (ynew**2) / ell) / Reff
    res = afactor * np.exp(-7.669 * (res0**0.25 - 1.0))
    soften = afactor * np.exp(-7.669 * ((0.1)**0.25 - 1.0))
    res[res > soften] = soften
    return res


def single_run_test(ind, ysc1, ysc2, q, vd, pha, zl, zs, lens_tag=1):
    nnn = 400  # Image dimension
    bsz = 9.0  # arcsecs
    dsx = bsz / nnn         # pixel size of SDSS detector.
    nstd = 0.001  # ^2

    xi1, xi2 = make_r_coor(nnn, dsx)
# ----------------------------------------------------------------------
    dsi = 0.03
    g_source = pyfits.getdata(
        "./gals_sources/439.0_149.482739_1.889989_processed.fits")
    g_source = np.array(g_source, dtype="<d")
    g_std = np.std(g_source)
    print g_std
    g_source[g_source <= 6.0*g_std] = 1e-6
# ----------------------------------------------------------------------
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = 0.0
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = 0.0
    rc = 0.0  # Core size of lens (in units of Einstein radius).
    re = re_sv(vd, zl, zs)  # Einstein radius of lens.
    re_sub = 0.0 * re
    a_sub = a_b_bh(re_sub, re)
    ext_shears = 0.06
    ext_angle = -0.39
    ext_kappa = 0.08

    # ext_shears = 0.0
    # ext_angle = 0.0
    # ext_kappa = 0.0
# ----------------------------------------------------------------------
    ai1, ai2 = psk.deflection_nie(xc1, xc2, pha, q, re, rc, ext_shears, ext_angle,
                                  ext_kappa, xi1, xi2)

    as1, as2 = psk.deflection_sub_pJaffe(0.0, -2.169, re_sub, 0.0, a_sub, xi1, xi2)

    yi1 = xi1 - ai1 - as1
    yi2 = xi2 - ai2 - as2

    g_limage = lv4.call_ray_tracing(g_source, yi1, yi2, ysc1, ysc2, dsi)
    # g_limage[g_limage <= 0.25] = 1e-6
# -------------------------------------------------------------
    afactor = 0.01
    Reff = 3.0
    ell = q
    ori = pha
    g_lens = de_vaucouleurs_2d(xi1, xi2, xc1, xc2, afactor, Reff, ell, ori)

    g_clean_ccd = g_lens + g_limage
    output_filename = "./fits_outputs/clean_lensed_imgs.fits"
    pyfits.writeto(output_filename, g_clean_ccd, overwrite=True)
# -------------------------------------------------------------
    from scipy.ndimage.filters import gaussian_filter
    g_images_psf = gaussian_filter(g_clean_ccd, 2.0)
# -------------------------------------------------------------
    g_noise = noise_map(nnn, nnn, np.sqrt(nstd), "Gaussian")
    output_filename = "./fits_outputs/noise_map.fits"
    pyfits.writeto(output_filename, g_noise, overwrite=True)
    g_final = g_images_psf + g_noise
# -------------------------------------------------------------
    g_clean_ccd = g_limage
    g_images_psf = gaussian_filter(g_clean_ccd, 2.0)
    g_final = g_images_psf + g_noise
    output_filename = "./fits_outputs/lensed_imgs_only.fits"
    pyfits.writeto(output_filename, g_final, overwrite=True)
# -------------------------------------------------------------
    output_filename = "./fits_outputs/full_lensed_imgs.fits"
    pyfits.writeto(output_filename, g_final + g_lens, overwrite=True)

    pl.figure()
    pl.contourf(g_final + g_lens)
    pl.colorbar()
# -------------------------------------------------------------
    # al11,al12 = np.gradient(al1,dsx)
    # al21,al22 = np.gradient(al2,dsx)
    # mua = 1.0/(1.0-(al11+al22)+al11*al22-al12*al21)

    return 0

if __name__ == '__main__':
    # from mpi4py import MPI
    # import sys
    # sourcpos = 10.0 # arcsecs
    # num_imgs = int(sys.argv[1])

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

    # import Image
    # I8 = (((I - I.min()) / (I.max() - I.min())) * 255.9).astype(np.uint8)
    # img = Image.fromarray(I8)
    # img.save("file.png")
    # for i in xrange(rank,num_imgs,size):
#-----------------------------------------------------------------------------------
    num_imgs = 1
    sourcpos = 0.0

    ysc1 = [0.4]    # X position of the source, arcsec
    ysc2 = [-0.3]   # Y position of the source, arcsec
    zl = 0.298      # the redshift of the lens galaxy
    zs = 1.0        # the redshift of the source galaxy
    vd = [320]      # Velocity Dispersion, km/s
    q = [0.64]      # Ellipticity
    pha = [73.0]    # Orintation, degree

    for i in xrange(num_imgs):
        single_run_test(i, ysc1[i], ysc2[i], q[i], vd[i], pha[i], zl, zs)

    pl.show()
