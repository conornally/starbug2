"""
A collection of plotting functions
"""
import os
import numpy as np
from astropy.visualization import ZScaleInterval
from scipy.interpolate import interp2d,RegularGridInterpolator
from multiprocessing import Pool

try: import matplotlib.pyplot as plt
except:
    from matplotlib import use; use("TkAgg")
    import maplotlib.pyplot
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap

from astropy.wcs import WCS
import astropy.units as u

import starbug2
from starbug2 import utils


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
        utils.perror("Unable to load style sheet \"%s\"\n"%fname)


def get_point_density(x,y,bins=30):
    hist,_x,_y=np.histogram2d(x,y,bins=bins)
    xx=np.linspace(min(x),max(x),bins)
    yy=np.linspace(min(y),max(y),bins)
    f=RegularGridInterpolator((xx,yy),hist)

    with Pool(processes=8) as pool:
        dens=pool.map(f,zip(x,y))
    return dens

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

    if not utils.wget("https://raw.githubusercontent.com/conornally/starbug2/refs/heads/main/docs/source/_static/images/starbug.png",fname="/tmp/sb.png"):
        ax.imshow(mpimg.imread("/tmp/sb.png"))
        ax.set_title("starbug2 v%s"%utils.get_version())
        ax.set_axis_off()
    return ax

def plot_cmd(tab, colour, mag, ax=None,col=None,hess=True,
        xlim=None, ylim=None, **kwargs):
    tt=utils.colour_index(tab,(colour,mag))
    mask=~(tt[colour].mask|tt[mag].mask)
    cc=tt[colour][mask]
    mm=tt[mag][mask]

    if not ax: fig,ax=plt.subplots(1)

    xm=np.nanmean(cc)
    dx=np.nanstd(cc)
    ym=np.nanmean(mm)
    dy=np.nanstd(mm)

    if xlim is None: xlim=(np.nanmin(cc),np.nanmax(cc))#( xm-(3*dx),xm+(5*dx))
    if ylim is None: ylim=(np.nanmin(mm),np.nanmax(mm))#( ym-(5*dy),ym+(3*dy))

    mask=((cc>=xlim[0])&(cc<=xlim[1]) & (mm>=ylim[0])&(mm<=ylim[1]))
    cc=cc[mask]
    mm=mm[mask]

    if col is None: col=plt.rcParams["axes.prop_cycle"].by_key()["color"][0]
    cmap=LinearSegmentedColormap.from_list("",[plt.rcParams["axes.prop_cycle"].by_key()["color"][0],col])

    if hess:
        bins=100
        hist,_x,_y=np.histogram2d(cc,mm,bins=bins)
        xx=np.linspace(min(cc),max(cc),bins)
        yy=np.linspace(min(mm),max(mm),bins)
        f=RegularGridInterpolator((xx,yy),hist)
        col=[f([X,Y]) for X,Y in zip(cc,mm)]
    pyplot_kw={"lw":0,"s":3}
    pyplot_kw.update(kwargs)
    ax.scatter(cc,mm,c=col,cmap=cmap,**pyplot_kw)




    ax.set_xlabel(colour)
    ax.set_ylabel(mag)
    ax.set_xlim(xlim)
    ax.set_ylim(*ylim[::-1])
    return ax




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
        if all(dat.shape):
            ax.imshow(ZScaleInterval()(dat), cmap="Greys_r", origin="lower")
            ax.text(0,0,im.header.get("FILTER"),c="white")

        ax.set_axis_off()
        fig.suptitle(src["Catalogue_Number"][0])
    fig.tight_layout()
        
    return fig

if __name__=="__main__":
    from astropy.table import Table
    fig,ax=plt.subplots(1)
    t=Table().read("/home/conor/sci/proj/ngc6822/overview/dat/ngc6822.fits")
    plot_cmd(t,"F115W-F200W","F200W",ax=ax)
    plt.show()





