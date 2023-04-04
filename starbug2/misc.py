"""
Miscillaneous functions...
"""
import os,stat,sys,numpy as np
import pkg_resources
import starbug2
from starbug2.utils import *
from starbug2.matching import sort_exposures
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

    printf("Downloading APPCORR CRDS files. NB: \x1b[1mTHESE MAY NOT BE THE LATEST!\x1b[0m\n")
    wget("https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_miri_apcorr_0005.fits", "%s/apcorr_miri.fits"%dname)
    wget("https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_nircam_apcorr_0004.fits", "%s/apcorr_nircam.fits"%dname)
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
    import webbpsf
    if os.getenv("WEBBPSF_PATH"): 
        dname=os.path.expandvars(dname)
        if not os.path.exists(dname):
            os.makedirs(dname)

        printf("Generating PSFs --> %s\n"%dname)

        load=loading(145, msg="initialising")
        load.show()
        for fltr,_f in starbug2.filters.items():
            if _f.instr==starbug2.NIRCAM:
                instr=webbpsf.NIRCam
                if _f.length==starbug2.SHORT: detectors=["NRCA1","NRCA2","NRCA3","NRCA4","NRCB1","NRCB2","NRCB3","NRCB4"]
                else: detectors=["NRCA5","NRCB5"]
            else:
                instr=webbpsf.MIRI
                detectors=[None]

            for det in detectors:
                model=instr()
                model.filter=fltr
                if det: model.detector=det
                load.msg="%6s %5s"%(fltr,det)
                #if det=="NRCA5": det="NRCALONG"
                #if det=="NRCB5": det="NRCBLONG"
                if det is None: det=""
                model.calc_psf()[1].writeto("%s/%s%s.fits"%(dname,fltr,det), overwrite=True)
                if '5' in det: 
                    _det=det.replace('5',"LONG")
                    model.calc_psf()[1].writeto("%s/%s%s.fits"%(dname,fltr,_det), overwrite=True)
                load();load.show()
                
            #nc = webbpsf.NIRCam() if line[5]==starbug2.NIRCAM else webbpsf.MIRI()
            #nc.filter=fltr
            #psf=nc.calc_psf()
            #load()
            #fits.PrimaryHDU(data=psf[1].data, header=psf[1].header).writeto("%s/%s.fits"%(dname, fltr), overwrite=True)
    else:perror("WARNING: Cannot generate PSFs, no environment variable 'WEBBPSF_PATH', please see https://webbpsf.readthedocs.io/en/latest/installation.html\n")

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


def update_paramfile(fname):
    """
    When the local parameter file is from an older version, add or remove the
    new or obselete keys
    INPUT: fname=local file to update
    """
    default_fname="%s/default.param"%pkg_resources.resource_filename("starbug2","param/")
    default_param=load_params(default_fname)
    current_param=load_params(fname)

    if default_fname==fname:
        #perror("cannot change default parameter file")
        return 

    if os.path.exists(fname):
        printf("Updating \"%s\"\n"%fname)
        fpi=open(fname, 'r')
        fpd=open(default_fname, 'r')
        fpo=open("/tmp/starbug.param",'w')

        add_keys=set(default_param.keys())-set(current_param.keys())
        del_keys=set(current_param.keys())-set(default_param.keys())
        if add_keys: printf("-> adding: %s  \n"%(', '.join(add_keys)))
        if del_keys: printf("-> removing: %s\n"%(', '.join(del_keys)))
        
        if not len(add_keys|del_keys): 
            printf("-> No updates needed\n")
            return 

        for inline in fpd.readlines():
            if inline[0] in "# \t\n":
                fpo.write(inline)
                continue

            key,value,comment=parse("{}={}//{}\n",inline)
            key=key.strip().rstrip()

            if key not in add_keys:
                value=current_param[key]

            outline="%-24s"%("%-12s"%key+"= "+str(value))+" //"+comment+"\n"
            fpo.write(outline)

        os.system("mv /tmp/starbug.param %s"%fname)
    else: perror("local parameter file '%s' does not exist\n"%fname)


