"""
Miscillaneous functions...
"""
import os,stat,sys,numpy as np
import starbug2
from starbug2.utils import *
from starbug2.matching import sort_exposures, GenericMatch
from astropy.io import fits



##########################
# One time run functions #
##########################
def init_starbug():
    """
    Initialise Starbug..
        - generate PSFs
        - download crds files
    INPUT:
        dname : data directory
    """
    printf("Initialising StarbugII\n")

    dname=starbug2.DATDIR
    printf("-> using %s=%s\n"%( "STARBUG_DATDIR" if os.getenv("STARBUG_DATDIR") else "DEFAULT_DIR", dname))
    generate_psfs()

    _miri_apcorr=  "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_miri_apcorr_0010.fits"
    _nircam_apcorr="https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_nircam_apcorr_0004.fits"

    printf("Downloading APPCORR CRDS files. NB: \x1b[1mTHESE MAY NOT BE THE LATEST!\x1b[0m\n")
    printf("-> %s\n"%_miri_apcorr)
    printf("-> %s\n"%_nircam_apcorr)
    wget(_miri_apcorr,   "%s/apcorr_miri.fits"%dname)
    wget(_nircam_apcorr, "%s/apcorr_nircam.fits"%dname)

    printf("Downloading ABVEGA offsets.\n")
    wget("https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_miri_abvegaoffset_0001.asdf","%s/abvegaoffset_miri.asdf"%dname)
    wget("https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_nircam_abvegaoffset_0002.asdf","%s/abvegaoffset_nircam.asdf"%dname)
    puts("Downloading The Junior Colour Encyclopedia of Space\n")


def generate_psfs():
    """
    Generate the psf files inside a given directory
    INPUT: 
        dname : directory to generate into
                if None then it will generate them into the folder given
                in loaded parameter file starbug2.DATDIR
    """
    dname=starbug2.DATDIR
    if os.getenv("WEBBPSF_PATH"): 
        dname=os.path.expandvars(dname)
        if not os.path.exists(dname):
            os.makedirs(dname)

        printf("Generating PSFs --> %s\n"%dname)

        load=loading(145, msg="initialising")
        load.show()
        for fltr,_f in starbug2.filters.items():
            if _f.instr==starbug2.NIRCAM:
                if _f.length==starbug2.SHORT: detectors=["NRCA1","NRCA2","NRCA3","NRCA4","NRCB1","NRCB2","NRCB3","NRCB4"]
                else: detectors=["NRCA5","NRCB5"]
            else:
                detectors=[None]

            for det in detectors:
                load.msg="%6s %5s"%(fltr,det)
                load.show()
                psf=generate_psf( fltr, det, None)
                if psf: 
                    psf.writeto("%s/%s%s.fits"%(dname, fltr, "" if det is None else det), overwrite=True)
                load()
                load.show()

    else:perror("WARNING: Cannot generate PSFs, no environment variable 'WEBBPSF_PATH', please see https://webbpsf.readthedocs.io/en/latest/installation.html\n")


def generate_psf(fltr, detector=None, fov_pixels=None):
    """
    Generate a single PSF for JWST
    INPUT:  fltr=JWST filter e.g. F444W
            detector=Instrument detector module e.g. NRCA1
            fov_pixels=size of PSF
    RETURNS:fits.HDUlist containing PSF
    """
    import webbpsf
    psf=None
    model=None
    if fov_pixels is not None and fov_pixels<=0: fov_pixels=None

    if fltr in list(starbug2.filters.keys()):
        _f=starbug2.filters.get(fltr)
        if detector is None:
            if _f.instr==starbug2.NIRCAM and _f.length==starbug2.SHORT: detector="NRCA1"
            elif _f.instr==starbug2.NIRCAM and _f.length==starbug2.LONG: detector="NRCA5"
            elif _f.instr==starbug2.MIRI: detector="MIRIM"

        if _f.instr==starbug2.NIRCAM: model=webbpsf.NIRCam()
        elif _f.instr==starbug2.MIRI: model=webbpsf.MIRI()

        if model:
            model.filter=fltr
            if detector: model.detector=detector
            try: 
                psf=model.calc_psf(fov_pixels=fov_pixels)["DET_SAMP"]
                psf=fits.PrimaryHDU(data=psf.data,header=psf.header)
            except: perror("\x1b[2KSomething went from with: %s %s\n"%(fltr,detector))
        else: perror("Unable to determing instrument from fltr '%s'\n"%fltr)
    else: perror("Unable to locate '%s' in JWST filter list\n"%fltr)
    return psf



def generate_runscript(fnames, args="starbug2 "):
    """
    """
    RUNFILE="./run.sh"
    fitsfiles=[]

    fp=open(RUNFILE,"w")
    fp.write("#!/bin/bash\n")
    fp.write("CMDS=\"-vf\"\n")
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
                    s=args +"${CMDS} -n%d "%len(exps)
                    for exp in exps: s+="%s "%exp[0].header["FILENAME"]
                    fp.write("%s\n"%s)
    fp.close()
    os.chmod(RUNFILE,stat.S_IXUSR|stat.S_IWUSR|stat.S_IRUSR|stat.S_IRGRP|stat.S_IROTH)
    printf("->%s\n"%RUNFILE)



def calc_instrumental_zeropint(psftable, aptable, fltr=None ):
    """

    """
    if fltr is None and not (fltr:=psftable.meta.get("FILTER")):
        perror("Unable to determine filter, set with '--set FILTER=F000W'.\n")
        return None
    printf("Calculating instrumental zeropoint %s.\n"%fltr)

    m=GenericMatch(threshold=0.1, colnames=["RA","DEC",fltr])
    matched=m([psftable,aptable], join_type="and")
    dist=np.array((matched["%s_2"%fltr]-matched["%s_1"%fltr]).value)
    zp=np.nanmedian(dist)
    std=np.nanstd(dist)
    printf("-> zp=%.3f +/- %.2g\n"%(zp,std))
    return (zp,std)
