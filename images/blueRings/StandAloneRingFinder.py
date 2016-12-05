import scipy.ndimage as  ndimage
import pylab as plt
#import moments
#import indexTricks as iT
import numpy

class RingFinder():
    def __init__(self,B,R,sB,sR,pB,pR,pixelsize=0.265,zerofluxB=1e12,zerofluxR=1e12,visualize=False,psfmode="worst"):
        self.pixelsize=pixelsize
        self.zerofluxB=zerofluxB
        self.zerofluxR=zerofluxR
        self.zeroflux=self.zerofluxB
        self.visualize=visualize

        x,y=numpy.mgrid[0:B.shape[0],0:B.shape[1]]
        xc=(B.shape[0]+1.)/2
        yc=(B.shape[1]+1.)/2
        rcent=((x-xc)**2+(y-yc)**2)**0.5
        dcent=rcent*pixelsize


        comR=getPeak(R,xc,yc)

        """
        #comB=self.getPeak(B,38,38)
        #if the peaks aren't aligned then we need to shift
        from scipy.ndimage import shift
        R=shift(R,(comR[0]-comB[0],comR[1]-comB[1]))
        """

        comB=comR# often no galaxy visible in the g-band. let's hope the astronemtry is good enough that we can do the centroiding in i
        print comR
        self.xl=comB[0]
        self.yl=comB[1]

        if numpy.abs(self.xl-38)>1 or numpy.abs(self.yl-38)>1:
            #sometimes the centroiding fails catastrophically - uncomment to imshow failures.
            #print "com fail"
            #plt.imshow(B,interpolation="none")
            #plt.show()
            #self.visualize=True
            pass


        #then match the psfs
        B,R,sB,sR,p=self.psfmatch(B,R,sB,sR,pB,pR,mode=psfmode)

        #import colorImage
        #color = colorImage.ColorImage()#sigma-clipped single band residual
        #cim=color.createModel(imgg,imgr,imgi)
        self.Rfull=R*1

        self.fitmask=dcent>5.
        R[dcent>5.]=0

        self.B=B
        self.R=R

        x,y=numpy.mgrid[0:B.shape[0],0:B.shape[1]]
        x-=int(self.xl)
        y-=int(self.yl)

        r2=(x**2+y**2)*(pixelsize**2)
        mask=((r2<2.7**2) & (r2>0.5**2))

        alpha=B[mask].sum()*1./R[mask].sum()

        self.D=B-alpha*self.Rfull
        self.S=(sB**2+(alpha*sR)**2)**.5

        self.SN=self.D/self.S
        self.Dshow=self.D*1.

        if self.visualize:
            plt.imshow(self.D,interpolation="none")
            plt.savefig("gminusi.png")
            plt.show()

    def psfmatch(self,B,R,sB,sR,pB,pR,mode="crossconvolve"):
        #not sure if this is totally valid, need to check
        if mode=="dont":
            return B,R,sB,sR,pR
        if mode=="crossconvolve":
            B=convolve(B,pR,True)[0]
            sB=convolve(sB,pR,True)[0]
            R=convolve(R,pB,True)[0]
            sR=convolve(sR,pB,True)[0]
            p=convolve(pB,pR,True)[0]
            return B,R,sB,sR,p
        if mode=="worst":
            pBcov=numpy.cov(pB)
            pRcov=numpy.cov(pR)
            plt.imshow(pBcov,interpolation=None)
            plt.show()

            Bf=convolve(B,pB,True)[1]
            Rf=convolve(R,pR,True)[1]

            #Bf[Bf==0]=1e-100+1e-100j
            #Rf[Rf==0]=1e-100+1e-100j

            plt.imshow(pB/pR,interpolation="None")
            plt.show()

            Bc=Bf.conj()
            BoR=Bf/Rf
            BoR[Bf==0]=-1e-99

            RoB=Rf/Bf
            #RoB2=numpy.real(RoB.conj()*RoB)
            #BoR2=numpy.real(BoR.conj()*BoR)
            #BB=numpy.real(Bf*Bc)
            #plt.imshow(RoB2,interpolation="None")
            #plt.show()
            #plt.imshow(BoR2,interpolation="None")
            #plt.show()

            R=convolve(R,BoR,False)[0]
            sR=convolve(sR,BoR,False)[0]
            p=pB
            return B,R,sB,sR,p

    def ringfind(self,cmax=2.7,vb=False):
        self.RegionFind(vb=vb)
        self.MomentChecks(cmax=cmax,vb=vb)
        if len(self.goodlabels)==0:
            if vb:print "found nothing"
            return False
        elif len(self.goodlabels)==1:
            if self.q<0.7:
                return True
            else:
                if vb:print "1 feature, not elongated"
                return False
        elif len(self.goodlabels)>1:
            return True
        else:
            if vb:print "shouldn't be called"
            return False

    def RegionFind(self,significancefloor=1.2,vb=False):

        self.significancefloor=significancefloor

        if self.visualize:
            plt.imshow(self.SN)
            plt.show()
        print self.SN.max()

        self.SN[self.SN<significancefloor]=0

        Amax=(7./self.pixelsize**2)

        regions, nlbl = ndimage.measurements.label(self.SN)

        if self.visualize:
            plt.imshow(self.SN)
            plt.show()


        lbls = numpy.arange(1, nlbl+1)

        for lbl in lbls:
            lblength= len(regions[regions==lbl])
            if lblength<10:
                self.SN[regions==lbl]=0
            elif lblength>Amax*1000:
                self.SN[regions==lbl]=0
            else:
                pass

        self.D[self.SN==0]=0

        self.regions, self.multiplicity = ndimage.measurements.label(self.D)
        self.labels=numpy.arange(1, self.multiplicity+1)


    def MomentChecks(self,cmax=2.7,cmin=0.5,vb=False):
        import moments
        import indexTricks as iT
        x,y=iT.coords(self.D.shape)
        x-=self.xl
        y-=self.yl
        x*=self.pixelsize
        y*=self.pixelsize

        goodlabels=[]
        unalignedlabels=[]

        for label in self.labels:
            D2=self.D*1
            D2[self.regions!=label]=0


            com,flux,mean=moments.comflux(D2) # we check this bit first since evaluating the shape is the slowest part
            ca=(com[1]-self.xl,com[0]-self.yl)
            com_dist=((ca[0]**2+ca[1]**2)**0.5)*self.pixelsize
            if com_dist < cmin:#arcseconds
                if vb:print "too close"
                continue
            if com_dist > cmax:#arcseconds
                if vb:print "too far"
                continue

            if self.FluxCheck(flux,mean)==False:
                if vb:print "too faint"
                continue

            com,(a,b,q), eigvals,eigvecs, flux,mean= moments.main(D2)

            ca_hat=ca/((ca[0]**2+ca[1]**2)**0.5)
            va=eigvecs[:,1]#normalized
            dp=ca_hat[0]*va[0]+ca_hat[1]*va[1]
            theta=(numpy.arccos(dp)*180./3.14159) #degrees

#this bit will let you visualize what's going on in each detected residual:
            if self.visualize:
                plt.imshow(D2,interpolation="None")
                plt.scatter(self.xl,self.yl,c="r",s=25)
                plt.plot([self.xl,self.xl+ca[0]],[self.yl,self.yl+ca[1]],c="r",lw=2)
                plt.scatter(com[1],com[0],c="k",s=25)
                plt.plot([com[1]-5*va[0],com[1]+5*va[0]],[com[0]-5*va[1],com[0]+5*va[1]],c="k",lw=3)
                plt.show(block=True)

            if b*self.pixelsize<0.2:#arcseconds
                if vb:print "too thin"
                continue

            if numpy.abs(theta) >30 and (theta<60 or theta >120):
                unalignedlabels.append(label)
            else:
                goodlabels.append(label)

            self.q=q

        self.goodlabels=goodlabels
        self.unalignedlabels=unalignedlabels

    def FluxCheck(self, flux,mean):
        zeroflux=self.zeroflux

        if flux<0:return False

        mag=10**(2.5*(flux*1./self.zeroflux))
        mag=-2.5*numpy.log10(flux*1./self.zeroflux)
        mean/=self.pixelsize**2 #in arcsec^2

        meanSB=10**(2.5*(mean*1./self.zeroflux))
        meanSB=-2.5*numpy.log10(mean*1./self.zeroflux)

        if mag > 25.5:return False
        if meanSB>26.3: return False
        return True

#============================================================================
# Functions from other parts of my code library that RingFinder needs.

def OLDgetPeak(img,x,y):
    from scipy import ndimage
    for i in range(3):
        d = ndimage.gaussian_filter(img[y-20:y+20,x-20:x+20],1.)[13:-13,13:-13]
        Y,X = numpy.mgrid[0:d.shape[0],0:d.shape[1]]
        Y -= Y.mean()
        X -= X.mean()
        d /= d.sum()
        x = x+int(numpy.round((d*X).sum()))
        y = y+int(numpy.round((d*Y).sum()))
    return x,y


def getPeak(img,x,y):
    from scipy import ndimage
    for i in range(3):
        d = ndimage.gaussian_filter(img,1.)
        Y,X = numpy.mgrid[0:d.shape[0],0:d.shape[1]]
        Y -= int(Y.mean())
        X -= int(X.mean())
        d /= d.sum()
        x = x+int(numpy.round((d*X).sum()))
        y = y+int(numpy.round((d*Y).sum()))
    return x,y


def convolve(image,psf,doPSF=True,edgeCheck=True):
    """
    A reasonably fast convolution routine that supports re-entry with a
    pre-FFT'd PSF. Returns the convolved image and the FFT'd PSF. Code written by  M. Auger.
    """
    datadim1 = image.shape[0]
    datadim2 = image.shape[1]
    if datadim1!=datadim2:
        ddim = max(datadim1,datadim2)
        s = numpy.binary_repr(ddim-1)
        s = s[:-1]+'0' # Guarantee that padding is used
    else:
        ddim = datadim1
        s = numpy.binary_repr(ddim-1)
    if s.find('0')>0:
        size = 2**len(s)
        if edgeCheck==True and size-ddim<8:
            size*=2
        boxd = numpy.zeros((size,size))
        r = size-datadim1
        r1 = r2 = r/2
        if r%2==1:
            r1 = r/2+1
        c = size-datadim2
        c1 = c2 = c/2
        if c%2==1:
            c1 = c/2+1
        boxdslice = (slice(r1,datadim1+r1),slice(c1,datadim2+c1))
        boxd[boxdslice] = image
    else:
        boxd = image

    if doPSF:
        # Pad the PSF to the image size
        boxp = boxd*0.
        if boxd.shape[0]==psf.shape[0]:
            boxp = psf.copy()
        else:
            r = boxp.shape[0]-psf.shape[0]
            r1 = r/2+1
            c = boxp.shape[1]-psf.shape[1]
            c1 = c/2+1
            boxpslice = (slice(r1,psf.shape[0]+r1),slice(c1,psf.shape[1]+c1))
            boxp[boxpslice] = psf.copy()
        # Store the transform of the image after the first iteration
        a = (numpy.fft.rfft2(boxp))
    else:
        a = psf
        # PSF transform and multiplication
    b = a*numpy.fft.rfft2(boxd)
    # Inverse transform, including phase-shift to put image back in center;
    #   this removes the requirement to do 2x zero-padding so makes things
    #   go a bit quicker.
    b = numpy.fft.fftshift(numpy.fft.irfft2(b)).real
    # If the image was padded, remove the padding
    if s.find('0')>0:
        b = b[boxdslice]

    return b,a


#============================================================================

if __name__ == '__main__':
 #all
 prefs=["DES2125+0001_3008740963","DES2139+0001_3019294847","DES2139-0041_3020241609","DES2139+0126_3018475310","DES2142-0124_3009067114","DES2145+0043_3020321773","DES2154+0043_3009679651","DES2154-0124_3019495911","DES2202+0001_3019592951","DES2211+0001_3009732561","DES2228-0124_3012190468","DES2237-0041_3013430036","DES2237-0124_3010953613","DES2254+0043_3012382083","DES2302+0043_3013612507","DES2311+0043_3014345997","DES2322-0041_3014901645","DES2328-0041_3016163832","DES2331+0001_3014942043","DES2331+0126_3015415030","DES2336-0041_3015515009","DES2356+0043_3016629464"]

 #nice
 prefs=["DES2125+0001_3008740963","DES2139+0001_3019294847","DES2139-0041_3020241609","DES2139+0126_3018475310","DES2142-0124_3009067114","DES2145+0043_3020321773","DES2154+0043_3009679651","DES2154-0124_3019495911","DES2202+0001_3019592951","DES2211+0001_3009732561","DES2228-0124_3012190468","DES2237-0124_3010953613","DES2254+0043_3012382083","DES2302+0043_3013612507","DES2311+0043_3014345997","DES2322-0041_3014901645","DES2328-0041_3016163832","DES2331+0001_3014942043","DES2331+0126_3015415030","DES2336-0041_3015515009","DES2356+0043_3016629464"]

 #trickys
 #prefs=["DES2142+0001_3009867289","DES2237-0041_3013430036","DES2311+0043_3014345997","DES2331+0001_3014942043","DES2336-0041_3015515009"]

 prefs=["DES2139-0041_3020241609","DES2154-0124_3019495911","DES2202+0001_3019592951","DES2211+0001_3009732561","DES2228-0124_3012190468","DES2237-0124_3010953613","DES2254+0043_3012382083","DES2302+0043_3013612507","DES2311+0043_3014345997","DES2322-0041_3014901645","DES2328-0041_3016163832","DES2331+0001_3014942043","DES2331+0126_3015415030","DES2336-0041_3015515009","DES2356+0043_3016629464"]



 count=0
 for pref in prefs:
  count+=1
  import pyfits
  xcent=[]
  ycent=[]
  import time
  t0=time.clock()

  a=75

  imgg=pyfits.open("../dessims/gal_fits/%s_gal_g.fits"%pref)[0].data.copy()[a:-a,a:-a]
  imgr=pyfits.open("../dessims/gal_fits/%s_gal_r.fits"%pref)[0].data.copy()[a:-a,a:-a]
  imgi=pyfits.open("../dessims/gal_fits/%s_gal_i.fits"%pref)[0].data.copy()[a:-a,a:-a]
  sigg=pyfits.open("../dessims/weight_fits/%s_weight_g.fits"%pref)[0].data.copy()[a:-a,a:-a]**-0.5
  sigr=pyfits.open("../dessims/weight_fits/%s_weight_r.fits"%pref)[0].data.copy()[a:-a,a:-a]**-0.5
  sigi=pyfits.open("../dessims/weight_fits/%s_weight_i.fits"%pref)[0].data.copy()[a:-a,a:-a]**-0.5


  psfg=pyfits.open("../dessims/psf_fits/%s_psf_g.fits"%pref)[0].data.copy()
  psfr=pyfits.open("../dessims/psf_fits/%s_psf_r.fits"%pref)[0].data.copy()
  psfi=pyfits.open("../dessims/psf_fits/%s_psf_i.fits"%pref)[0].data.copy()


  #plt.imshow(imgg)
  #plt.show()

  #make colour image:
  import colorImage
  color = colorImage.ColorImage()
  colorimage = color.createModel(imgg,imgr,imgi)
  #plt.imshow(colorimage,interpolation="none")
  #plt.savefig("colour.png"%pref)
  #plt.show()


  RF=RingFinder(imgg,imgi,sigg,sigi,psfg,psfi,0.265,1e12,1e12,visualize=False,psfmode="crossconvolve")
  RFres=RF.ringfind(vb=True)

  #ax=plt.subplot(1,3,count)

  ax=plt.subplot(1,1,1)
  size=colorimage.shape[0]

  b=a*1
  a=55
  rf=float(228-(75*2))/float(228-(25*2))
  print rf


  imgg1=pyfits.open("../dessims/gal_fits/%s_gal_g.fits"%pref)[0].data.copy()[a:-a,a:-a]
  imgr1=pyfits.open("../dessims/gal_fits/%s_gal_r.fits"%pref)[0].data.copy()[a:-a,a:-a]
  imgi1=pyfits.open("../dessims/gal_fits/%s_gal_i.fits"%pref)[0].data.copy()[a:-a,a:-a]
  color2 = colorImage.ColorImage()
  colorimage2 = color2.createModel(imgg1,imgr1,imgi1)
  plt.imshow(colorimage2,interpolation="none",extent=[0,size,0,size])

  RF.Dshow[RF.fitmask]=numpy.nan
  colorresid=color.colorize(RF.Dshow,RF.Dshow,RF.Dshow)
  d=0
  plt.imshow(RF.Dshow,interpolation="none",extent=[0-d,rf*size-d,0-d,rf*size-d])
  plt.xlim([0,size-20])
  plt.ylim([0,size-20])
  ax.xaxis.set_visible(False)
  ax.yaxis.set_visible(False)

  from numpy import cos,sin,pi
  t=numpy.linspace(0,2*pi)
  r=19.5*rf
  x=(size+1)/2.+r*cos(t)-1#-d
  y=(size-1)/2.+r*sin(t)#-d
  plt.plot(x,y,lw=2,c="w")

  plt.show()
  plt.cla()
  #if count>2:break
  continue
  xcent.append(RF.xl)
  ycent.append(RF.yl)
  #try:
  #    ax=plt.subplot(4,6,2*count-1)
  #except IndexError: break
  plt.imshow(colorimage,interpolation="none")
  ax.xaxis.set_visible(False)
  ax.yaxis.set_visible(False)
  ax=plt.subplot(4,6,2*count)
  ax.xaxis.set_visible(False)
  ax.yaxis.set_visible(False)

  print RF.fitmask
  print RF.Dshow.shape
  RF.Dshow[RF.fitmask]=numpy.nan

  colorresid=color.colorize(RF.Dshow,RF.Dshow,RF.Dshow)
  #plt.imshow(colorresid,interpolation="none")
  #plt.text(5,15,RFres)
  plt.imshow(RF.Dshow,interpolation="none")

  plt.show()



"""
  # check the peak finding here:
  #print numpy.median(xcent), numpy.median(ycent)
 #plt.scatter(xcent,ycent,s=1)
  #plt.show()
  #exit()


  import time
  t0=time.clock()
  for i in range(1000):
    i=457
    imgg=pyfits.open("../LensPop/MockImages/Candidate_%i_g_SDSS_img.fits"%i)[0].data.copy()
    sigg=pyfits.open("../LensPop/MockImages/Candidate_%i_g_SDSS_sig.fits"%i)[0].data.copy()
    psfg=pyfits.open("../LensPop/MockImages/Candidate_%i_g_SDSS_psf.fits"%i)[0].data.copy()
    imgi=pyfits.open("../LensPop/MockImages/Candidate_%i_i_SDSS_img.fits"%i)[0].data.copy()
    sigi=pyfits.open("../LensPop/MockImages/Candidate_%i_i_SDSS_sig.fits"%i)[0].data.copy()
    psfi=pyfits.open("../LensPop/MockImages/Candidate_%i_i_SDSS_psf.fits"%i)[0].data.copy()

    RF=RingFinder(imgg,imgi,sigg,sigi,psfg,psfi,0.265,1e12,1e12,visualize=False)
    RFres=RF.ringfind(vb=False)
  print (time.clock()-t0)*3*1000./(60**2)

  """
