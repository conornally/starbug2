"""
Miscillaneous functions...
"""
import os,stat,sys,numpy as np
import starbug2
from starbug2.utils import *
from starbug2.matching import sort_exposures
from astropy.io import fits



##########################
# One time run functions #
##########################

def generate_psfs(dname):
    """
    Generate the psf files inside a given directory
    INPUT: 
        dname : directory to generate into
                if None then it will generate them into the folder given
                in loaded parameter file "PSFDIR"
    """
    import webbpsf

    dname=os.path.expandvars(dname)
    if not os.path.exists(dname):
        os.mkdir(dname)

    printf("Generating PSFs --> %s\n"%dname)

    load=loading(len(starbug2.filters))
    for fltr,line in starbug2.filters.items():
        load.msg=fltr
        load.show()
        
        nc = webbpsf.NIRCam() if line[5]==starbug2.NIRCAM else webbpsf.MIRI()
        nc.filter=fltr
        psf=nc.calc_psf()
        load()
        fits.PrimaryHDU(data=psf[1].data, header=psf[1].header).writeto("%s/%s.fits"%(dname, fltr), overwrite=True)
    load.show()

def generate_runscript(fnames, args="starbug2 "):
    """
    """
    RUNFILE="./run.sh"
    fitsfiles=[]

    fp=open(RUNFILE,"w")
    fp.write("#!/bin/bash\n")
    for fname in fnames:
        if os.path.exists(fname):
            dname,name,ext=split_fname(fname)
            if ext==".fits":
                fitsfile=fits.open(fname)
                fitsfile[0].header["FILENAME"]=fname
                fitsfiles.append(fitsfile)

            else: perror("file %s must be type '.fits' not '%s'\n"%(name,ext))
        else: perror("file \x1b[1;31m%s\x1b[0m not found\n"%fname)

    sorted=sort_exposures(fitsfiles)

    for band, obs in sorted.items():
        for ob,visits in obs.items():
            for visit,dets in visits.items():
                for det,exps in dets.items():
                    s=args +"-n%d "%len(exps)
                    for exp in exps: s+="%s "%exp[0].header["FILENAME"]
                    fp.write("%s\n"%s)
    fp.close()
    os.chmod(RUNFILE,stat.S_IXUSR|stat.S_IWUSR|stat.S_IRUSR|stat.S_IRGRP|stat.S_IROTH)
    printf("->%s\n"%RUNFILE)

