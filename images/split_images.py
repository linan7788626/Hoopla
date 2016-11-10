#import cv2
#import pylab as pl

#file_name = "./ering.jpg"
#re_img1 = cv2.imread(file_name)
#b, g, r = cv2.split(re_img1)

from PIL import Image
import numpy as np

a = Image.open("./ering.jpg")
a = np.array(a)
a[:,:,0] *=1
a[:,:,1] *=0
a[:,:,2] *=0
a = Image.fromarray(a)
a.show()
a.save('lensed_galaxy.jpg')

