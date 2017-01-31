import numpy as np


def deflection_sub_pJaffe(x0, y0, re, rc, a, x, y):

    r = np.sqrt((x - x0)**2.0 + (y - y0)**2.0)
    res = re / r * (np.sqrt(rc * rc + r * r) - rc - np.sqrt(a * a + r * r) + a)
    return res * (x - x0) / r, res * (y - y0) / r


def deflection_nie(x0, y0, theta, ql, re, rc, ext_shears, ext_angle, ext_kappa, x, y):  # SIE lens model
    tr = np.pi * (theta / 180.0)  # + np.pi / 2.0

    sx = x - x0
    sy = y - y0

    cs = np.cos(tr)
    sn = np.sin(tr)

    sx_r = sx * cs + sy * sn
    sy_r = -sx * sn + sy * cs

    psi = np.sqrt(ql**2.0 * (rc**2.0 + sx_r**2.0) + sy_r**2.0)
    dx_tmp = (re * np.sqrt(ql) / np.sqrt(1.0 - ql**2.0)) * \
        np.arctan(np.sqrt(1.0 - ql**2.0) * sx_r / (psi + rc))
    dy_tmp = (re * np.sqrt(ql) / np.sqrt(1.0 - ql**2.0)) * \
        np.arctanh(np.sqrt(1.0 - ql**2.0) * sy_r / (psi + rc * ql**2.0))
    dx = dx_tmp * cs - dy_tmp * sn
    dy = dx_tmp * sn + dy_tmp * cs

    # external shear
    tr2 = np.pi * (ext_angle / 180.0)
    cs2 = np.cos(2.0 * tr2)
    sn2 = np.sin(2.0 * tr2)
    dx2 = ext_shears * (cs2 * sx + sn2 * sy)
    dy2 = ext_shears * (sn2 * sx - cs2 * sy)

    # external kappa
    dx3 = ext_kappa * sx
    dy3 = ext_kappa * sy
    return dx + dx2 + dx3, dy + dy2 + dy3


def potential_nie(x0, y0, theta, ql, re, rc, ext_shears,
                  ext_angle, ext_kappa, x, y):
    tr = np.pi * (theta / 180.0)  # + np.pi / 2.0
    sx = x - x0
    sy = y - y0
    cs = np.cos(tr)
    sn = np.sin(tr)
    sx_r = sx * cs + sy * sn
    sy_r = -sx * sn + sy * cs
    psi = np.sqrt(ql**2.0 * (rc**2.0 + sx_r**2.0) + sy_r**2.0)
    dx_tmp = (re * np.sqrt(ql) / np.sqrt(1.0 - ql**2.0)) * \
        np.arctan(np.sqrt(1.0 - ql**2.0) * sx_r / (psi + rc))
    dy_tmp = (re * np.sqrt(ql) / np.sqrt(1.0 - ql**2.0)) * \
        np.arctanh(np.sqrt(1.0 - ql**2.0) * sy_r / (psi + rc * ql**2.0))
    pot_SIE = sx_r * dx_tmp + sy_r * dy_tmp - 0.5 * re * \
        np.sqrt(ql) * rc * np.log((psi + rc)**2.0 +
                                  (1.0 - (ql**2.0)) * (sx_r**2.0))

    # external shear
    tr2 = np.pi * (ext_angle / 180.0)
    cs2 = np.cos(2.0 * tr2)
    sn2 = np.sin(2.0 * tr2)
    pot_exts = ext_shears * (sn2 * sx * sy + 0.5 * cs2 * (sx**2.0 - sy**2.0))

    # external kappa
    pot_kaps = ext_kappa * (sx**2.0 + sy**2.0) * 0.5
    return pot_SIE + pot_exts + pot_kaps


def pot_sub_pJaffe(x0, y0, re, a, x, y):
    # res =
    # re*(np.sqrt(s*s+r*r)-s-np.sqrt(a*a+r*r)+a)-re*s*np.log((s+np.sqrt(s*s+r*r)/(2.0*s)))+re*a*np.log((a+np.sqrt(a*a+r*r)/(2.9*a)))

    r = np.sqrt((x - x0)**2.0 + (y - y0)**2.0)

    res = re * (r - np.sqrt(a * a + r * r) + a) + re * a * \
        np.log((a + np.sqrt(a * a + r * r)) / (2.0 * a))
    return res

#if __name__ == '__main__':
    #import pylab as pl

    #r = np.linspace(0.0, 3.0, 100)
    #y  = 1.0 / r * (np.sqrt(r * r) - np.sqrt(1.0 * 1.0 + r * r) + 1.0)
    #pl.plot(r,y,'k-')
    #pl.show()
