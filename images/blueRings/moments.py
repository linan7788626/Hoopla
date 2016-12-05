import numpy
def moments(image,coords):
    x, y = coords
    x=x.ravel()
    y=y.ravel()
    image=(image*1).ravel()

    mass= sum(image)
    x_bar = sum(x*image)/sum(image)
    y_bar = sum(y*image)/sum(image)
    m00 = sum(image)
    m01 = sum(x*image)
    m10 = sum(y*image)
    m11 = sum(y*x*image)
    m02 = sum(x**2*image)
    m20 = sum(y**2*image)
    m12 = sum(x*y**2*image)
    m21 = sum(x**2*y*image)
    m03 = sum(x**3*image)
    m30 = sum(y**3*image)

    u11 = (m11 - x_bar * m01) / mass 
    u20 = (m20 - x_bar * m10) / mass
    u02 = (m02 - y_bar * m01) / mass

    cov = np.array([[u20, u11], [u11, u02]])

    eigvals, eigvecs = numpy.linalg.eigh(cov)
    com=(x_bar,y_bar)
    return com,eigvals,eigvecs

import numpy as np
import matplotlib.pyplot as plt

def raw_moment(data, iord, jord):
        nrows, ncols = data.shape
        y, x = np.mgrid[:nrows, :ncols]
        data = data * x**iord * y**jord
        return data.sum()

def com(data):
    return main(data,comonly=True)
def comflux(data):
    return main(data,comflux=True)

def main(data,comonly=False,comflux=False):
    data_sum = data.sum()
    m10 = raw_moment(data, 1, 0)
    m01 = raw_moment(data, 0, 1)
    x_bar = m10 / data_sum
    y_bar = m01 / data_sum

    flux=data_sum
    mean=data_sum/float(len(data[data!=0]))

    if comonly:return (y_bar,x_bar)
    if comflux:return (y_bar,x_bar),flux,mean

    u11 = (raw_moment(data, 1, 1) - x_bar * m01) / data_sum
    u20 = (raw_moment(data, 2, 0) - x_bar * m10) / data_sum
    u02 = (raw_moment(data, 0, 2) - y_bar * m01) / data_sum
    cov = np.array([[u20, u11], [u11, u02]])

    eigvals, eigvecs = np.linalg.eigh(cov)

    #mean = np.array([x_bar, y_bar])

    a=eigvals[1]**.5
    b=eigvals[0]**.5
    q=b/a



    com=(y_bar,x_bar)
    return com,(a,b,q), eigvals,eigvecs, flux,mean


if __name__ == '__main__':
    com,(a,b,q), eigvecs=main()
    va=eigvecs[:,1]
    xl=50
    yl=50

    ca=(com[0]-xl,com[1]-yl)
    ca/=(ca[0]**2+ca[1]**2)

    theta=numpy.arccos(ca[0]*va[0]+ca[1]*va[1])
    print theta


