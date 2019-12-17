'''
Simulate EM field propagation with fourier transforms


'''

import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from scipy.fftpack import fft2,fftshift


def create_pupil(res=512, radius=40, obscuration=0., norm=1., geom="tele"):
    '''
    Create Telescope like pupil.

    res = resolution elements in both x and y direction
    radius = radius of pupil in resolution elements
    obscuration = radius of central obscuration in resolution elements
    norm = normalize integrated intensity
    '''
    y,x=np.ogrid[-res/2:res/2,  -res/2:res/2]
    cy, cx=np.meshgrid(x,y)

    aper_dist=np.sqrt(cx**2+cy**2)
    ampl=np.zeros_like(cx)

    if (geom=="tele"):
        ampl[aper_dist <= radius] = 1.
        ampl[aper_dist < obscuration*radius] = 0.
    elif (geom=="gaussian"):
        ampl = np.exp(-1* aper_dist**2 / (2 * radius**2) )

    ampl*=norm/np.sum(ampl)
    phase=np.zeros_like(cx)

    pupil=ampl*np.exp(1j*phase*3)
    return pupil

def plotfield(field, fig=None):
    '''
    Plot EM field.

    field = field as complex array
    '''
    fig, ax= plt.subplots(2,2, num=fig)

    imampl = ax[0,0].imshow(np.abs(field))
    ax[0,0].set_label("Amplitude")
    plt.colorbar(imampl, ax=ax[0,0])

    imphas = ax[0,1].imshow(np.arctan2(np.imag(field),np.real(field)))
    ax[0,1].set_label("Phase")
    plt.colorbar(imphas, ax=ax[0,1])

    imreal = ax[1,0].imshow(np.real(field))
    ax[1,0].set_label("Real part")
    plt.colorbar(imreal, ax=ax[1,0])

    imimag = ax[1,1].imshow(np.imag(field))
    ax[1,1].set_label("Imag part")
    plt.colorbar(imimag, ax=ax[1,1])

def fft(field):
    farfield=fft2(fftshift(field))
    return(fftshift(farfield))

#### Create telescope pupil and focus it
#### Create gaussian and focus it

telepupil=create_pupil()
#plotfield(telepupil, fig="Pupil")
gausspupil=create_pupil(geom="gaussian")
#plotfield(gausspupil, fig="Gaussian")
#fft!
telefocus=fft(telepupil)
#plotfield(telefocus, fig="Focus Pupil")
gaussfocus=fft(gausspupil)
#plotfield(gaussfocus, fig="Focus Gaussian")


#### Do fiber coupling efficiency calculations for different obscurations and gaussian sizes. Also plot it

obsclist=np.linspace(0.0,1.0,50)
gausslist=np.linspace(20,40,25)
coupling=np.zeros_like(gausslist)
couplingobsc=np.zeros_like(obsclist)
for iio, obsc in enumerate(obsclist):
    telepupil=create_pupil(obscuration=obsc)
    #couplingf=np.zeros_like(gausslist)
    for iis, sigma in enumerate(gausslist):
        gausspupil=create_pupil(geom="gaussian", radius=sigma)
        # gaussfocus=fft(gausspupil)
        coupling[iis]=np.real(np.sum(telepupil*np.conjugate(gausspupil))**2 / np.sum(telepupil**2) / np.sum(gausspupil**2))
        # couplingf[iis]=np.real(np.sum(telefocus*np.conjugate(gaussfocus))**2 / np.sum(telefocus**2) / np.sum(gaussfocus**2))
        # print("Coupling Pupil/Focus",
        #     np.sum(telepupil*np.conjugate(gausspupil))**2 / np.sum(telepupil**2) / np.sum(gausspupil**2),
        #     np.sum(telefocus*np.conjugate(gaussfocus))**2 / np.sum(telefocus**2) / np.sum(gaussfocus**2)  )

    plt.figure("Coupling efficiency")
    plt.plot(gausslist, coupling, label="obsc "+str(obsc))
    couplingobsc[iio]=np.max(coupling)
    #plt.plot(gausslist, couplingf)
plt.legend()

#Plot it like Ruilier 1998
plt.figure("Coupling vs. Obscuration")
plt.plot(obsclist, couplingobsc)

plt.show()
