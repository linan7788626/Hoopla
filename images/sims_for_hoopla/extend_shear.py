import numpy as np


def deflect_SIE(lens, x, y):  # SIE lens model
    x0 = 0.0
    y0 = 0.0
    pi = 3.14

    tr = pi * (th / 180.0) + pi / 2.0

    sx = x − x0
    sy = y − y0

    cs = cos(tr)
    sn = sin(tr)
    sx_r = sx * cs + sy * sn
    sy_r = −sx * sn + sy * cs
    psi = sqrt(lens.fl**2.0 * (lens.rc**2.0 + sx_r**2.0) + sy_r**2.0)
    dx_tmp = (lens.bl * sqrt(lens.fl) / sqrt(1.0−lens.fl**2.0)) * arctan(sqrt(1.0−lens.fl**2.0) * sx_r / (psi + lens.rc))
    dy_tmp = (lens.bl * sqrt(lens.fl) / sqrt(1.0−lens.fl**2.0)) * arctanh(sqrt(1.0−lens.fl**2.0) * sy_r / (psi + lens.rc * lens.fl**2.0))
    dx = dx_tmp * cs − dy_tmp * sn
    dy = dx_tmp * sn + dy_tmp * cs
    # external shear
    tr2 = pi * (lens.sa / 180.0)
    cs2 = cos(2.0 * tr2)
    sn2 = sin(2.0 * tr2)
    dx2 = lens.ss * (cs2 * sx + sn2 * sy)
    dy2 = lens.ss * (sn2 * sx − cs2 * sy)
    return array([dx + dx2, dy + dy2])

#############################################################################
# Potential for SIE + external shear #
#############################################################################


def potential(lens, x, y):  # SIE lens model
    tr = np.pi * (lens.th / 180.0) + np.pi / 2.0
    sx = x − lens.x0
    sy = y − lens.y0
    cs = cos(tr)
    sn = sin(tr)
    sx_r = sx * cs + sy * sn
    sy_r = −sx * sn + sy * cs
    psi = sqrt(lens.fl**2.0 * (lens.rc**2.0 + sx_r**2.0) + sy_r**2.0)
    dx_tmp = (lens.bl * sqrt(lens.fl) / sqrt(1.0−lens.fl**2.0)) * arctan(sqrt(1.0−lens.fl**2.0) * sx_r / (psi + lens.rc))
    dy_tmp = (lens.bl * sqrt(lens.fl) / sqrt(1.0−lens.fl**2.0)) * arctanh(sqrt(1.0−lens.fl**2.0) * sy_r / (psi + lens.rc * lens.fl**2.0))
    pot_SIE = sx_r * dx_tmp + sy_r * dy_tmp − 0.5 * lens.bl * sqrt(lens.fl) * lens.rc * log((psi + lens.rc)**2.0 + (1.0−(lens.fl**2.0)) * (sx_r**2.0))

    # external shear
    tr2 = pi * (lens.sa / 180.0)
    cs2 = cos(2.0 * tr2)
    sn2 = sin(2.0 * tr2)
    pot_exts = lens.ss * (sn2 * sx * sy + 0.5 * cs2 * (sx**2.0−sy**2.0))
    pot_kaps = lens.kk * (sx**2.0 + sy**2.0) * 0.5
    return pot_SIE + pot_exts


def phi_shears(x1, x2):
    res = shear1 / 2.0 * (x1 * x1 - x2 * x2) + shear2 * x1 * x2
    return res


def phi_kappa(x1, x2):
    res = kappa / 2.0 * (x1 * x1 + x2 * x2)
    return res


def potential_shear_kappa(x1, x2):
    phi_all = phi_sie(x1, x2) + phi_shear(x1, x2) + phi_kappa(x1, x2)
    return phi_all
