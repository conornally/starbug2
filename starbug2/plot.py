"""
A collection of plotting functions
"""
import os
import numpy as np
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from astropy.wcs import WCS
import astropy.units as u

import starbug2
from starbug2.utils import printf, perror, wget, get_version


def load_style(fname):
    """
    Load a pyplot style sheet
    
    Parameters
    ----------
    fname : str
        Filename of the style sheet
    """
    if os.path.exists(fname): 
        plt.style.use(fname)
    else:
        perror("Unable to load style sheet \"%s\"\n"%fname)

def plot_test(ax, **kwargs):
    """
    Just plot the starbug image
    
    Parameters
    ----------
    ax : plt.axes
        Ax to plot into

    Returns
    -------
    ax : plt.axes
        The working axes
    """

    if not wget("https://raw.githubusercontent.com/conornally/starbug2/refs/heads/main/docs/source/_static/images/starbug.png",fname="/tmp/sb.png"):
        ax.imshow(mpimg.imread("/tmp/sb.png"))
        ax.set_title("starbug2 v%s"%get_version())
        ax.set_axis_off()
    return ax

def plot_cmd(*args,**kwargs):
    pass

def plot_inspectsource(src, images):
    """
    Show a source in an array of images

    Parameters
    ----------
    src : astropy.Table Row
        Input source to look at

    images : list 
        List of fits images to inspect

    Returns
    -------
    fig : plt.figure
        The figure
    """

    n=len(images)
    fig,axs=plt.subplots(1,n, figsize=(1.7*n,2))
    if n==1: axs=[axs]
    images=sorted(images,key=lambda a: list(starbug2.filters.keys()).index(a.header["FILTER"]))


    size=0.1#arcsec?
    for n,(im,ax) in enumerate(zip(images,axs)):
        wcs=WCS(im)
        x,y=wcs.all_world2pix(np.ones(2)*src["RA"], np.array([1,1+(size/3600)])*src["DEC"], 0)
        dp=np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)#abs(y[1]-y[0])#20#

        xmin=max(0,int(np.floor(x[0]-dp)))
        xmax=min(im.data.shape[1]-1,int(np.ceil(x[0]+dp)))
        ymin=max(0,int(np.floor(y[0]-dp)))
        ymax=min(im.data.shape[0]-1,int(np.ceil(y[0]+dp)))


        dat=im.data[min(ymin,ymax):max(ymin,ymax),min(xmin,xmax):max(xmin,xmax)]
        #print(xmin,xmax,ymin,ymax,dp,dat.shape)
        #if not all(dat.shape): 
            #print("invert")
            #dat=im.data[xmin:xmax,ymin:ymax] #HACK because MIRI has inverted coords?

        if all(dat.shape):
            ax.imshow(ZScaleInterval()(dat), cmap="Greys_r", origin="lower")
            ax.text(0,0,im.header.get("FILTER"),c="white")



        #ax.tick_params(labelleft=False,labelbottom=False)
        ax.set_axis_off()
        fig.suptitle(src["Catalogue_Number"][0])
    fig.tight_layout()
        
    return fig






