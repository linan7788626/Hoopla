import glob,numpy,pylab,pyfits,cPickle

class ColorImage:
    def __init__(self):
        self.rnorm = None
        self.gnorm = None
        self.I = None
        self.M = None
        self.nonlin = 5.
        self.m = 0.5
        self.bMinusr = 0.8
        self.bMinusg = 0.4


    def clip(self,arr,nsig=3.5):
        a = arr.flatten()
        a.sort()
        a = a[a.size*0.05:a.size*0.8]
        m,s,l = a.mean(),a.std(),a.size
        while 1:
            a = a[abs(a-m)<s*nsig]
            if a.size==l:
                return m,s
            m,s,l = a.mean(),a.std(),a.size


    def createModel(self,b,g,r):
        bMinusr = self.bMinusr
        bMinusg = self.bMinusg
        b0 = b.copy()
        g0 = g.copy()
        r0 = r.copy()
        
        w = r.shape[0]/2-5
        rb = r0/b0
        gb = g0/b0
        rnorm = numpy.median(rb[w:-w,w:-w])
        gnorm = numpy.median(gb[w:-w,w:-w])
        r0 /= rnorm
        g0 /= gnorm
        r0 *= 10**(0.4*bMinusr)
        g0 *= 10**(0.4*bMinusg)

        r0 /= 620.
        g0 /= 540.
        b0 /= 460.

        I = (r0+g0+b0)/3.
        self.I = I
        self.rnorm = rnorm
        self.gnorm = gnorm
        return self.colorize(b,g,r)

    def colorize(self,b,g,r,newI=False):
        bMinusr = self.bMinusr
        bMinusg = self.bMinusg
        rnorm = self.rnorm
        gnorm = self.gnorm
        m = self.m
        nonlin = self.nonlin
        I = self.I.copy()

        b = b.copy()
        g = g.copy()
        r = r.copy()

        w = r.shape[0]/2-5
        r /= rnorm
        g /= gnorm
        r *= 10**(0.4*bMinusr)
        g *= 10**(0.4*bMinusg)

        r /= 620.
        g /= 540.
        b /= 460.

        sdev = self.clip(I)[1]
        m = m*sdev
        if self.M is None:
            M = I[w:-w,w:-w].max()
        else:
            M = self.M
        nonlin = nonlin*sdev

        if newI==True:
            I = (b+g+r)/3.
        f = numpy.arcsinh((I-m)/nonlin)/numpy.arcsinh((M-m)/nonlin)
        f[I<m] = 0.
        f[I>M] = 1.
        R = r*f/I
        G = g*f/I
        B = b*f/I

        R[I<=0] = 0.
        G[I<=0] = 0.
        B[I<=0] = 0.

        R[R<=0] = 0.
        G[G<=0] = 0.
        B[B<=0] = 0.

        R[R>1] = 1.
        G[G>1] = 1.
        B[B>1] = 1.

        white = True
        if white:
            cond = (f==1)
            R[cond] = 1.
            G[cond] = 1.
            B[cond] = 1.

        arr = numpy.empty((R.shape[0],R.shape[1],3))
        arr[:,:,0] = R
        arr[:,:,1] = G
        arr[:,:,2] = B

        return arr
