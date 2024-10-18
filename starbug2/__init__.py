import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter("ignore",category=AstropyWarning)
warnings.simplefilter("ignore",category=RuntimeWarning) ## bit dodge that

logo="""
                   *          *  __  *  __   - * --   - 
 STARBUGII               *      / ___ /    \   --  -      - 
 ---------                 *___---.    .___/  -   --   -
 JWST photometry in       ./=== \  \.     \      * 
 complex crowded fields   | (O)  |  |     |           *
                           \._._/ ./    _(\)   *   
 conor.nally@ed.ac.uk     /   ~--\ ----~   \      *
                        ---      ___       ---      
 > %s
"""

motd="https://starbug2.readthedocs.io/en/latest/"

from os import getenv
_=getenv("STARBUG_DATDIR") 
DATDIR=_ if _ else "%s/.local/share/starbug"%(getenv("HOME"))


## HASHDEFS
MIRI=1
NIRCAM=2

NULL=0
LONG=1
SHORT=2

PIX=0
ARCSEC=1
ARCMIN=2
DEG=3

## SOURCE FLAGS
SRC_GOOD=0
SRC_BAD=0x01
SRC_JMP=0x02
SRC_VAR=0x04 ##source frame mean >5% differnet than median
SRC_FIX=0x08 ##psf fit with fixed centroid
SRC_UKN=0x10 ##source unknown


##DQ FLAGS
DQ_DO_NOT_USE=0x01
DQ_SATURATED =0x02
DQ_JUMP_DET  =0x04

## DEFAULT MATCHING COLS
match_cols=["RA","DEC","flag","flux","eflux", "NUM"]#"stdflux", "NUM"]

# ZERO POINT...
ZP={   "F070W"	:[3631,0],
       "F090W"	:[3631,0],
       "F115W"	:[3631,0],
       "F140M"	:[3631,0],
       "F150W"	:[3631,0],
       "F162M"	:[3631,0],
       "F164N"	:[3631,0],
       "F150W2" :[3631,0],
       "F182M"	:[3631,0],
       "F187N"	:[3631,0],
       "F200W"	:[3631,0],
       "F210M"	:[3631,0],
       "F212N"	:[3631,0],
       "F250M"	:[3631,0],
       "F277W"	:[3631,0],
       "F300M"	:[3631,0],
       "F322W2" :[3631,0],
       "F323N"	:[3631,0],
       "F335M"	:[3631,0],
       "F356W"	:[3631,0],
       "F360M"	:[3631,0],
       "F405N"	:[3631,0],
       "F410M"	:[3631,0],
       "F430M"	:[3631,0],
       "F444W"	:[3631,0],
       "F460M"	:[3631,0],
       "F466N"	:[3631,0],
       "F470N"	:[3631,0],
       "F480M"	:[3631,0],

       "F560W"  :[3631,0],
       "F770W"  :[3631,0],
       "F1000W" :[3631,0],
       "F1130W" :[3631,0],
       "F1280W" :[3631,0],
       "F1500W" :[3631,0],
       "F1800W" :[3631,0],
       "F2100W" :[3631,0],
       "F2550W" :[3631,0],
       }

class _F: #(struct) containing JWST filter info
    def __init__(self, wavelength, aFWHM, pFWHM, instr, length):
        self.wavelength=wavelength
        self.aFWHM=aFWHM
        self.pFWHM=pFWHM
        self.instr=instr
        self.length=length

filters={ #08/06/2023
"F070W":_F(0.704   ,0.023,0.742,NIRCAM,SHORT),
"F090W":_F(0.901   ,0.030,0.968,NIRCAM,SHORT),
"F115W":_F(1.154   ,0.037,1.194,NIRCAM,SHORT),
"F140M":_F(1.404   ,0.046,1.484,NIRCAM,SHORT),
"F150W":_F(1.501   ,0.049,1.581,NIRCAM,SHORT),
"F162M":_F(1.626   ,0.053,1.710,NIRCAM,SHORT),
"F164N":_F(1.644   ,0.054,1.742,NIRCAM,SHORT),
"F150W2":_F(1.671  ,0.045,1.452,NIRCAM,SHORT),
"F182M":_F(1.845   ,0.060,1.935,NIRCAM,SHORT),
"F187N":_F(1.874   ,0.061,1.968,NIRCAM,SHORT),
"F200W":_F(1.990   ,0.064,2.065,NIRCAM,SHORT),
"F210M":_F(2.093   ,0.068,2.194,NIRCAM,SHORT),
"F212N":_F(2.120   ,0.069,2.226,NIRCAM,SHORT),
"F250M":_F(2.503   ,0.082,1.302,NIRCAM,LONG),
"F277W":_F(2.786   ,0.088,1.397,NIRCAM,LONG),
"F300M":_F(2.996   ,0.097,1.540,NIRCAM,LONG),
"F322W2":_F(3.247  ,0.096,1.524,NIRCAM,LONG),
"F323N":_F(3.237   ,0.106,1.683,NIRCAM,LONG),
"F335M":_F(3.365   ,0.109,1.730,NIRCAM,LONG),
"F356W":_F(3.563   ,0.114,1.810,NIRCAM,LONG),
"F360M":_F(3.621   ,0.118,1.873,NIRCAM,LONG),
"F405N":_F(4.055   ,0.132,2.095,NIRCAM,LONG),
"F410M":_F(4.092   ,0.133,2.111,NIRCAM,LONG),
"F430M":_F(4.280   ,0.139,2.206,NIRCAM,LONG),
"F444W":_F(4.421   ,0.140,2.222,NIRCAM,LONG),
"F460M":_F(4.624   ,0.151,2.397,NIRCAM,LONG),
"F466N":_F(4.654   ,0.152,2.413,NIRCAM,LONG),
"F470N":_F(4.707   ,0.154,2.444,NIRCAM,LONG),
"F480M":_F(4.834   ,0.157,2.492,NIRCAM,LONG),

"F560W":_F(5.589   ,0.207,1.882,MIRI,NULL),
"F770W":_F(7.528   ,0.269,2.445,MIRI,NULL),
"F1000W":_F(9.883  ,0.328,2.982,MIRI,NULL),
"F1130W":_F(11.298 ,0.375,3.409,MIRI,NULL),
"F1280W":_F(12.712 ,0.420,3.818,MIRI,NULL),
"F1500W":_F(14.932 ,0.488,4.436,MIRI,NULL),
"F1800W":_F(17.875 ,0.591,5.373,MIRI,NULL),
"F2100W":_F(20.563 ,0.674,6.127,MIRI,NULL),
"F2550W":_F(25.147 ,0.803,7.300,MIRI,NULL),
}


helpstrings={
"DETECTION":"""
    Source Detection
    ----------------

    This routine locates point sources in an image. The input is a
    FITS image and the output is a FITS table, containing a list of
    point source locations, their geometric properties and flux/magnitude
    measurements as calculated by aperture photometry. The output file
    will have the suffix "-ap", note this is the same as the output
    for the aperture photometry routine.

    To run this routine, use the core command:

        $~ starbug2 -D image.fits

    Alter the parameter file options under "DETECTION" to tune the 
    performance of starbug2. Two of the key parameters are:
    
        - SIGSKY : Set the background level of the image
        - SIGSRC : Set the detection threshold of the sources

    Full documentation is at https://starbug2.readthedocs.io

""",

"BACKGROUND":"""
    Diffuse Background Estimations
    ------------------------------

    This routine estimates the "dusty" emissions in an image, given
    a source list. It is used to subtract from the image, thus removing
    the flux contribution on a source brightness from the dusty environment.
    
    The routine requires a list of sources to be generated (by source 
    detection) or loaded with [-d sourcelist.fits] and requires a FITS
    image to work on. The routine will ouput a FITS image, with the same
    dimensions and spatial coverage as the input image, with the suffix
    "-bgd". This background image can be used in the photometry later.

    To run the routine, use the core command:

        $~ starbug2 -B -d sourcelist.fits image.fits

    Alter the parameter file options under "BACKGROUND ESTIMATION" 
    to tune the performance of starbug2. Two key parameters are:

        - BGD_R    : Set a fixed aperture mask radius around each source
        - BOX_SIZE : Set the estimation resolution (larger will be more blurred)

    Full documentation is at https://starbug2.readthedocs.io

""",

"APPHOT":"""
    Aperture Photometry
    -------------------

    This routine conducts aperture photometry on an image given a list
    of sources. It requires a FITS image to run on and a FITS table source
    list with either RA/DEC columns, or x/y_centroid or x/y_0 columns. 
    The routine outputs a table with the suffix "-ap". Note this filename
    is the same as the source detection routine because aperture photometry 
    is automatically run at the end of the source detection step. The output
    table contains 2flux/magnitude information on every source

    To run this routine, use the core command:

        $~ starbug2 -A -d sourcelist.fits image.fits

    Alter the parameter file options under "APERTURE PHOTOMETRY" to tune
    the performance of starbug2. Three key parameters are:

        - APPHOT_R : Set the aperture radius for photometry (in pixels)
        - SKY_RIN  : Set the inner sky annulus radius (in pixels)
        - SKY_ROUT : Set the outer sky annulus radius (in pixels)

    Full documentation is at https://starbug2.readthedocs.io

""",

"PSFPHOT":"""
    PSF Photometry
    --------------

    This routine conducts PSF fitting photometry on an image given
    a list of sources. Its requires a FITS image to run on and a FITS table
    sourcelist with either RA/DEC columns, or x/y_centroid or x/y_0 columns.
    The routine outputs a table with the suffix "-psf". The output table 
    contains 2flux/magnitude information on every source
    
    To run this routine, use the core command:
        
        $~ starbug2 -P -d sourcelist image.fits 

    Alter the parameter file options under "PHOTOMETRY" to tune the
    performance of starbug2. Two key parameters are:

        - FORCE_POS    : Hold the cetroid positions of source fixed (forced photometry)
        - GEN_RESIDUAL : Generate a residual image from all the fit source

    Full documentation is at https://starbug2.readthedocs.io

""",

"MATCHOUTPUTS":"""
    Match Outputs
    -------------

    This option is set if the user wishes to combine all the output catalogues
    from starbug together. It would be used in the case that a routine is
    being ran on a list of images (either in series or parallel) and the
    final catalogues should all be combined into a single source list.
    It outputs two files, one with the suffix "full" and another with "match".
    The first is all columns from all table preserved into a single large
    catalogue, the second averages all the similar columns into a reduced
    table.

    To run this routine, use the core code:
        
        $~ starbug2 -DM image1.fits image2.fits image3.fits ...

    Alter the parameter file options under "CATALOGUE MATCHING" to tune the
    performance of starbug2. Two key parameters are:

        - MATCH_THRESH : Set the separation threshold (arcsec) to match two sources
        - NEXP_THRESH  : Set the minimum number of catalogues a source must be present in

""",
}


