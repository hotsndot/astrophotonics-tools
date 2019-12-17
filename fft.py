'''
Simulate EM field with fourier transforms


'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft2,fftshift


def createpupil(resolution=128, radius=50, obscuration=0):
    y,x=np.ogrid[-res/2:res/2+1,  -res/2:res/2+1]
    cy, cx=np.meshgrid(x,y)
    aper=np.sqrt(cx**2+cy**2) <= radius
    return aperture



aper=createpupil()
