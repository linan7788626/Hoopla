import StandAloneRingFinder as sarf
from astropy.io import fits as pyfits
import pylab as pl
import colorImage

pref=1
ddir="example" #Change to your data directory
imgg=pyfits.getdata("%s/clip_%s_g.fits"%(ddir,pref))
imgr=pyfits.getdata("%s/clip_%s_r.fits"%(ddir,pref))
imgi=pyfits.getdata("%s/clip_%s_i.fits"%(ddir,pref))
imgz=pyfits.getdata("%s/clip_%s_z.fits"%(ddir,pref))

sigg=pyfits.getdata("%s/clip_%s_weight_g.fits"%(ddir,pref))
sigr=pyfits.getdata("%s/clip_%s_weight_r.fits"%(ddir,pref))
sigi=pyfits.getdata("%s/clip_%s_weight_i.fits"%(ddir,pref))
sigz=pyfits.getdata("%s/clip_%s_weight_z.fits"%(ddir,pref))

psfg=pyfits.getdata("%s/clip_%s_psf_g.fits"%(ddir,pref))
psfr=pyfits.getdata("%s/clip_%s_psf_r.fits"%(ddir,pref))
psfi=pyfits.getdata("%s/clip_%s_psf_i.fits"%(ddir,pref))
psfz=pyfits.getdata("%s/clip_%s_psf_z.fits"%(ddir,pref))


color = colorImage.ColorImage()
colorimage = color.createModel(imgg,imgr,imgi)

print colorimage.shape
print colorimage.max()
print colorimage.min()

psfmode="dont"
RF=sarf.RingFinder(imgg,imgi,sigg,sigi,psfg,psfi,0.265,1e12,1e12,visualize=False,psfmode=psfmode)
RFgz=sarf.RingFinder(imgg,imgz,sigg,sigi,psfg,psfz,0.265,1e12,1e12,visualize=False,psfmode=psfmode)
RFrz=sarf.RingFinder(imgr,imgz,sigr,sigz,psfr,psfz,0.265,1e12,1e12,visualize=False,psfmode=psfmode)
RFiz=sarf.RingFinder(imgi,imgz,sigi,sigz,psfi,psfz,0.265,1e12,1e12,visualize=False,psfmode=psfmode)

pl.figure(figsize=(15,9))
ax=pl.subplot(1,1,1)
size=colorimage.shape[0]

pl.subplot(131)
color.bMinusr = 0.8
color.bMinusg = 0.4
color.nonlin = 1.
colored=color.createModel(imgg,imgr,imgi)
pl.imshow(colored,interpolation="none")
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
pl.subplot(132)
colorresid=color.colorize(RFgz.Dshow,RFrz.Dshow,RFiz.Dshow)
d=0

pl.imshow(colorresid,interpolation="none")
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

pl.subplot(133)
pl.imshow(RF.Dshow,interpolation="none")
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

pl.show()
