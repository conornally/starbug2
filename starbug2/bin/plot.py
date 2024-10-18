"""StarbugII Plotting Scripts
usage: starbug2-plot [-vhX] [-I CN000] [-o outfile] images.fits ..
    -h  --help           : show help screen
    -o  --output   fname : output filename
    -v  --verbose        : verbose mode

    -I  --inspect  CN000 : inspect a source in an array of images
    -X  --test           : plot a test image

        --style    fname : load a custom pyplot style sheet
        --dark           : plot in dark mode
"""
import os,sys,getopt
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

import starbug2.bin as scr
import starbug2
from starbug2.plot import load_style, plot_test, plot_inspectsource
from starbug2.utils import printf, perror, warn

VERBOSE =0x01
SHOWHELP=0x02
STOPPROC=0x04
KILLPROC=0x08

DARKMODE=0x10

PTEST=   0x1000
PINSPECT=0x2000


def plot_parseargv(argv):
    options=0
    setopt={}
    cmd,argv = scr.parsecmd(argv)
    opts,args = getopt.gnu_getopt(argv, "hvXI:d:o:",("help","verbose","test",
                                                "inspect=",
                                                "output=", "style=", "dark"))

    for opt,optarg in opts:
        match(opt):
            case "-h"|"--help":     options|=(SHOWHELP|STOPPROC)
            case "-v"|"--verbose":  options|=VERBOSE
            case "-o"|"--output":   setopt["OUTPUT"]=optarg
            case "-d"|"--apfile":   setopt["APFILE"]=optarg

            case "-I"|"--inspect":
                options|=PINSPECT
                setopt["INSPECT"]=optarg
            case "-X"|"--test": options|=PTEST

            case "--style": setopt["STYLESHEET"]=optarg
            case "--dark":  options|=DARKMODE

    return options, setopt, args


def plot_onetimeruns(options, setopt, args):
    if options&SHOWHELP:
        scr.usage(__doc__,verbose=options&VERBOSE)

        if options & PINSPECT: perror(fn_pinspect.__doc__)

        return scr.EXIT_EARLY

    if (_fname:=setopt.get("STYLESHEET")):
        load_style(_fname)

    if options&DARKMODE: 
        load_style("%s/extras/dark.style"%starbug2.__path__[0])
    
    if options & STOPPROC: return scr.EXIT_EARLY
    if options & KILLPROC:
        perror("..killing process\n")
        return scr.EXIT_FAIL

    return scr.EXIT_SUCCESS

def fn_pinspect(options, setopt, images=None, tables=None):
    """
    Inspect Source
    --------------

    Plot at a source position cutouts in a range of images. 
    This requires a source list to be loaded, a list of image 
    file and the source catalogue number to be given. This will 
    take the form::

        $~ starbug2-plot -I CN123 sourcelist.fits image*.fits
    """
    """
    Parameters
    ----------
    options : int
        The starbug2.bin.plot options integar

    setopt : dict
        The starbug2.bin.plot setopt dictionary

    images : list
        The list of fits image HDUs to cut out from 

    tables : list
        The source list to pull the source from. Must have a column
        with the name "Catalogue_Number"

    Returns
    -------
    fig : plt.figure
        The output figure
    """
    fig=None
    if (cn:=setopt.get("INSPECT")) and images and tables:
        if "Catalogue_Number" in tables[0].colnames and cn in tables[0]["Catalogue_Number"]:
            i=np.where(tables[0]["Catalogue_Number"]==cn)[0]
            fig=plot_inspectsource(tables[0][i], images)

    else: perror("Must include the source Catalogue_Number, a list of images and a sourcelist \n")
    return fig

def plot_main(argv):
    warn("Still in development\n\n")
    options,setopt,args=plot_parseargv(argv)
    load_style("%s/extras/starbug.style"%starbug2.__path__[0])
    exit_code=0

    if options or setopt:
        if (exit_code:=plot_onetimeruns(options, setopt, args)):
            return exit_code

    images=[]
    tables=[]
    for arg in args:
        if (_fname:=os.path.exists(arg)):
            fp=fits.open(arg)
            _filter=fp[0].header.get("FILTER") # THIS IS A HACK
            for hdu in fp:
                if hdu.header.get("XTENSION")=="IMAGE":
                    images.append(hdu)
                    break
                if hdu.header.get("XTENSION")=="BINTABLE":
                    tables.append(Table(hdu.data))
                    break
            hdu.header["FILTER"]=_filter

    fig=None
    if options& PTEST:
        fig,ax=plt.subplots(1,figsize=(3,2.5))
        ax=plot_test(ax)

    if options& PINSPECT: fig=fn_pinspect(options, setopt, images=images, tables=tables)

    if fig is not None:
        fig.tight_layout()
        if (output:=setopt.get("OUTPUT")):
            fig.savefig(output, dpi=300)
        else:
            plt.show()

def plot_mainentry():
    """Command Line entry point"""
    return plot_main(sys.argv)
