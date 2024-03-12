import os
from parse import parse
from starbug2.utils import printf,perror,get_version

default="""## STARBUG CONFIG FILE
# Generated with starbug2-v%s
PARAM       =  STARBUGII PARAMETERS     //COMMENT


## GENERIC
VERBOSE     = 0          // (0:false 1:true)
OUTPUT      = .          //Directory or filename to output to 
HDUNAME     =            //If using a non standard HDU name, name it here (str or int)
FILTER      =            //Set a custom filter for the image

## DETECTETION 
FWHM        = -1         // Custom FWHM for image (-1 to use WEBBPSF)
SIGSKY      = 2.0        // (float>0) 
SIGSRC      = 5.0        // (float>0) Source value mininmum N sigma above background
DOBGD2D     = 1          // Run background2D step (usually finds more sources but takes time)
DOCONVL     = 1          // Run convolution step (usually finds more sources)
CLEANSRC    = 1          // Run source cleaning after detection (removes likely contaminants)
SHARP_LO    = 0.4        // cutoff in detection
SHARP_HI    = 0.9        // cutoff in detection
ROUND1_HI   = 1.0        // cutoff in detection
ROUND2_HI   = 1.0        // cutoff in detection
SMOOTH_LO   =            //.
SMOOTH_HI   =            //.
RICKER_R    = 1.0        // Radius (pix) of ricker wavelet 

## APERTURE PHOTOMOETRY
APPHOT_R    = 1.5        //Radius in number of pixels
ENCENERGY   = -1         //Fraction encircled energy (mutually exclusive with APPHOT_R)
SKY_RIN     = 3          //Sky annulus inner radius
SKY_ROUT    = 4.5        //Sky annulus outer radius
APCORR_FILE =            //Aperture correction file. See full manual for details

## BACKGROUND ESTIMATION
BGD_R       = 0          //.
BOX_SIZE    = 2          // (int>0) Background estimation kernal size (pix)

## PHOTOMETRY
AP_FILE     =            //Detection file to use instead of detecting
BGD_FILE    =            //Background estimation file
PSF_FILE    =            //Non default PSF file
USE_WCS     = 1          //When loading an AP_FILE, do you want to use WCS or xy values (if available)
ZP_MAG      = 8.9        //Zero point (mag) to add to the magnitude columns 

CRIT_SEP    = 8          //minimum distance for grouping (pixels) between two sources
FORCE_POS   = 0          //Force centroid position (1) or allow psf fitting to fit position too (0)
DPOS_THRESH = -1         //If allowed to fit position, max separation (arcsec) from source list centroid
MAX_XYDEV   = 3p         //.
PSF_SIZE    = -1         //Set fit size of psf (>0) or -1 to take PSF file dimensions
GEN_RESIDUAL= 0          //generate a residual image

## SOURCE STATS
CALC_CROWD  = 1          //Run crowding metric calculation (execution time scales N^2)

## CATALOGUE MATCHING
MATCH_THRESH= 0.1        // when combining background subtraction catalogue, minimum separation (arcsec) of centroids to be considered separate sources
MATCH_COLS  =            // EXTRA columns to include in output matched table i.e sharpness
NEXP_THRESH = -1         // Keep sources that appear in NUM >= NEXP_THRESH (if -1 keep everything)
SN_THRESH   = -1         // Remove sources with SN ratio < SN_THRESH before matching (default -1 to not apply this cut)
BRIDGE_COL  =            // Bridge --band matching NIRCam and MIRI catalogues by ensuring NIRCam catalogue has a match in BRIDGE_COL

## ARTIFICAL STAR TESTS
NTESTS      = 100        // Number of artificial star tests
NSTARS      = 1          // Number of stars per artifical test
SUBIMAGE    = 500        //number of pixels ? to crop around artificial star
MIN_FLUX    = 0.00001      //minimun flux for artificial star
MAX_FLUX    = 100.0      //maximum flux for artificial star

## MISC EXTRAS
REGION_COL  = green      //DS9 region colour
REGION_SCAL = 1          //Scale region to flux if possible
REGION_RAD  = 3          //Region radius default
REGION_XCOL = RA         //X column name to use for region
REGION_YCOL = DEC        //Y column name to use for region
REGION_WCS  = 1          //If X/Y column names correspind to WCS values
"""%get_version()

def parse_param(line):
    """
    Parse a parameter line
    """
    param={}
    if line and line[0] not in "# \t\n":
        if "//" in line: key,value,_=parse("{}={}//{}",line)
        else: key,value=parse("{}={}",line)
        key=key.strip().rstrip()
        value=value.strip().rstrip()
        try:
            if '.' in value: value=float(value)
            else: value=int(value)
        except:
            pass

        ## Special case values
        if key in ("OUTPUT", "AP_FILE","BGD_FILE","PSF_FILE"): value=os.path.expandvars(value)
        param[key]=value
    return param



def load_default_params():
    config={}
    for line in default.split('\n'):
        config.update(parse_param(line))
    return config

def load_params(fname):
    """
    Convert a parameter file into a dictionary of options
    INPUT:  fname=path/to/file.param
    RETURN: dictionary of options
    """
    config={}
    if fname is None:
        config=load_default_params()
    elif os.path.exists(fname):
        with open(fname, "r") as fp:
            for line in fp.readlines():
                config.update(parse_param(line))
    else:
        perror("config file \"%s\" does not exist\n"%fname)
    return config

def local_param():
    with open("starbug.param", "w") as fp:
        fp.write(default)

def update_paramfile(fname):
    """
    When the local parameter file is from an older version, add or remove the
    new or obselete keys
    INPUT: fname=local file to update
    """
    default_param=load_default_params()
    current_param=load_params(fname)

    if os.path.exists(fname):
        printf("Updating \"%s\"\n"%fname)
        fpi=open(fname, 'r')
        fpo=open("/tmp/starbug.param",'w')

        add_keys=set(default_param.keys())-set(current_param.keys())
        del_keys=set(current_param.keys())-set(default_param.keys())
        if add_keys: printf("-> adding: %s  \n"%(', '.join(add_keys)))
        if del_keys: printf("-> removing: %s\n"%(', '.join(del_keys)))
        
        if not len(add_keys|del_keys): 
            printf("-> No updates needed\n")
            return 

        for inline in default.split("\n"):
            if inline and inline[0] not in "# \t\n":

                key,value,comment=parse("{}={}//{}",inline)
                key=key.strip().rstrip()

                if key not in add_keys:
                    value=current_param[key]
                outline="%-24s"%("%-12s"%key+"= "+str(value))+" //"+comment
            else: outline=inline

            fpo.write("%s\n"%outline)
        fpi.close()
        fpo.close()
        os.system("mv /tmp/starbug.param %s"%fname)
    else: perror("local parameter file '%s' does not exist\n"%fname)

