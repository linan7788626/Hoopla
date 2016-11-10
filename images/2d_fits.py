import scipy.optimize as opt
import numpy as np
import pylab as plt
from PIL import Image
import json
import sys

def load_guessed_pars(filename):
#with open('./test.JSON') as data_file:
    with open(filename) as data_file:
        data = json.load(data_file)

    #---------------------------------------------------
    # Guessed lens
    #

    lx_guess = data['components'][1]['x']
    ly_guess = data['components'][1]['y']
    lr_guess = data['components'][1]['theta_e']

    if data['components'][1]['ell'] > 1:
        lq_guess = 1.0/data['components'][1]['ell']
    else:
        lq_guess = data['components'][1]['ell']

    lpa_guess = data['components'][1]['ang']

    #---------------------------------------------------
    # Guessed source
    #

    sx_guess = data['components'][0]['x']
    sy_guess = data['components'][0]['y']
    sr_guess = data['components'][0]['size']

    if data['components'][0]['ell'] > 1:
        sq_guess = 1.0/data['components'][0]['ell']
    else:
        sq_guess = data['components'][0]['ell']

    spa_guess = data['components'][0]['ang']

    return (lx_guess, ly_guess, lr_guess, lq_guess, lpa_guess, sx_guess,
            sy_guess, sr_guess, sq_guess, spa_guess)

def save_improved_pars(filename1, par_array, filename2):

    with open(filename1) as data_file:
        data = json.load(data_file)

    data['name'] = "Improved Model of "+data['name']
    #---------------------------------------------------
    # Improved lens
    #

    data['components'][1]['x'] = par_array[0]
    data['components'][1]['y'] = par_array[1]
    data['components'][1]['theta_e'] = par_array[2]

    data['components'][1]['ell'] = par_array[3]

    data['components'][1]['ang'] = par_array[4]+90

    #---------------------------------------------------
    # Improved source
    #

    data['components'][0]['x'] = par_array[5]
    data['components'][0]['y'] = par_array[6]
    data['components'][0]['size'] = par_array[7]

    data['components'][0]['ell'] = par_array[8]
    data['components'][0]['ang'] = par_array[9]+90

    with open(filename2, "w") as outfile:
        json.dump(data, outfile, indent=4)


def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

def lensing_signals_sie(x1, x2, xc1, xc2, q, rc, re, pha):

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

    #xt11 = cosa
    #xt22 = cosa
    #xt12 = sina
    #xt21 = -sina

    #fx11 = xt11 / pd1 - xt1 * \
        #(xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd1 * pd1)
    #fx22 = xt22 / pd2 - xt2 * \
        #(xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd2 * pd2)
    #fx12 = xt12 / pd1 - xt1 * \
        #(xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd1 * pd1)
    #fx21 = xt21 / pd2 - xt2 * \
        #(xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd2 * pd2)

    #a11 = qs / (1.0 + fx1 * fx1) * fx11
    #a22 = qs / (1.0 - fx2 * fx2) * fx22
    #a12 = qs / (1.0 + fx1 * fx1) * fx12
    #a21 = qs / (1.0 - fx2 * fx2) * fx21

    #rea11 = (a11 * cosa - a21 * sina) * re
    #rea22 = (a22 * cosa + a12 * sina) * re
    #rea12 = (a12 * cosa - a22 * sina) * re
    #rea21 = (a21 * cosa + a11 * sina) * re

    #kappa = 0.5 * (rea11 + rea22)
    #shear1 = 0.5 * (rea12 + rea21)
    #shear2 = 0.5 * (rea11 - rea22)

    #y11 = 1.0 - rea11
    #y22 = 1.0 - rea22
    #y12 = 0.0 - rea12
    #y21 = 0.0 - rea21

    #jacobian = y11 * y22 - y12 * y21
    #mu = 1.0 / jacobian

    alpha1 = (a1 * cosa - a2 * sina) * re
    alpha2 = (a2 * cosa + a1 * sina) * re

    return alpha1, alpha2 #, kappa, shear1, shear2, mu

def create_images((x1,x2),xc1,xc2,re,ql,pha,yc1,yc2,sigs,qs,phs):
    sig1 = sigs/qs*0.693
    sig2 = sigs*qs*0.693
    al1, al2 = lensing_signals_sie(x1, x2, xc1, xc2, ql, 0.0, re, pha)

    y1 = x1 - al1
    y2 = x2 - al2
    res = twoD_Gaussian((y1, y2), 1.0, yc1, yc2, sig1, sig2, phs, 0.0)
    return res

def loadin_images(filename):
    a = Image.open(filename)
    a = np.array(a)
    #res = a[:,:,0]/np.max(a[:,:,0])
    #res[res<0] = 0.0
    res = a[:,:,0]/256.0
    return res.ravel()


if __name__ == '__main__':
    # Create x and y indices
    nnn = 400
    dsx = 0.0225
    bsz = dsx*nnn
    x1 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
    x2 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
    x1, x2 = np.meshgrid(x1, x2)

    fname = "./ering.jpg"
    data = loadin_images(fname)

    fname_json = sys.argv[1]
    initial_guess = load_guessed_pars(fname_json)

    popt, pcov = opt.curve_fit(create_images, (x1, x2), data, p0=initial_guess)

    data_guessed = create_images((x1, x2), *initial_guess)

    levels = [0.5]
    fig, ax = plt.subplots(1, 1)
    ax.hold(True)
    ax.imshow(data.reshape(nnn, nnn), cmap=plt.cm.winter, origin='bottom',
        extent=(x1.min(), x1.max(), x2.min(), x2.max()))
    ax.contour(x1, x2, data_guessed.reshape(nnn, nnn), levels, colors='r')

    data_fitted = create_images((x1, x2), *popt)

    save_improved_pars(fname_json, popt, "Improved_"+fname_json)

    fig, ax = plt.subplots(1, 1)
    ax.hold(True)
    ax.imshow(data.reshape(nnn, nnn), cmap=plt.cm.winter, origin='bottom',
        extent=(x1.min(), x1.max(), x2.min(), x2.max()))
    ax.contour(x1, x2, data_fitted.reshape(nnn, nnn), levels, colors='r')

    plt.figure()
    plt.contourf(x1,x2,(data_fitted-data).reshape(nnn, nnn), cmap=plt.cm.jet)
    plt.colorbar()
    plt.show()
