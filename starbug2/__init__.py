import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter("ignore",category=AstropyWarning)

logo="""\
                   *          *  __  *  __   - * --   - 
 STARBUGII               *      / ___ /    \   --  -      - 
 ---------                 *___---.    .___/  -   --   -
 JWST photometry in       ./=== \  \.     \      * 
 complex crowded fields   | (O)  |  |     |           *
                           \._._/ ./    _(\)   *   
 conor.nally@ed.ac.uk     /   ~--\ ----~   \      *
                        ---      ___       ---      
 > %s"""

from os import getenv
_=getenv("STARBUG_DATDIR") 
DATDIR=_ if _ else "%s/.local/share/starbug"%(getenv("HOME"))


## HASHDEFS
MIRI=1
NIRCAM=2

NULL=0
LONG=1
SHORT=2

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
match_cols=["RA","DEC","flag","flux","eflux", "stdflux", "crowding"]

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

#            name    wavelen FWHM(arcsec,pix) instrument length
filters={   "F070W"	:_F(0.704,  0.030,0.987, NIRCAM,SHORT),
            "F090W"	:_F(0.902,  0.034,1.103, NIRCAM,SHORT),
            "F115W"	:_F(1.154,  0.040,1.298, NIRCAM,SHORT),
            "F140M"	:_F(1.405,  0.048,1.553, NIRCAM,SHORT),
            "F150W"	:_F(1.501,  0.050,1.628, NIRCAM,SHORT),
            "F162M"	:_F(1.627,  0.055,1.770, NIRCAM,SHORT),
            "F164N"	:_F(1.645,  0.056,1.801, NIRCAM,SHORT),
            "F150W2":_F(1.659,  0.046,1.494, NIRCAM,SHORT),
            "F182M"	:_F(1.845,  0.062,1.990, NIRCAM,SHORT),
            "F187N"	:_F(1.874,  0.064,2.060, NIRCAM,SHORT),
            "F200W"	:_F(1.989,  0.066,2.141, NIRCAM,SHORT),
            "F210M"	:_F(2.095,  0.071,2.304, NIRCAM,SHORT),
            "F212N"	:_F(2.121,  0.072,2.341, NIRCAM,SHORT),
            "F250M"	:_F(2.503,  0.084,1.340, NIRCAM,LONG),
            "F277W"	:_F(2.762,  0.091,1.444, NIRCAM,LONG),
            "F300M"	:_F(2.989,  0.100,1.585, NIRCAM,LONG),
            "F322W2":_F(3.232,  0.097,1.547, NIRCAM,LONG),
            "F323N"	:_F(3.237,  0.108,1.711, NIRCAM,LONG),
            "F335M"	:_F(3.362,  0.111,1.760, NIRCAM,LONG),
            "F356W"	:_F(3.568,  0.115,1.830, NIRCAM,LONG),
            "F360M"	:_F(3.624,  0.120,1.901, NIRCAM,LONG),
            "F405N"	:_F(4.052,  0.136,2.165, NIRCAM,LONG),
            "F410M"	:_F(4.082,  0.137,2.179, NIRCAM,LONG),
            "F430M"	:_F(4.281,  0.145,2.300, NIRCAM,LONG),
            "F444W"	:_F(4.408,  0.145,2.302, NIRCAM,LONG),
            "F460M"	:_F(4.630,  0.155,2.459, NIRCAM,LONG),
            "F466N"	:_F(4.654,  0.158,2.507, NIRCAM,LONG),
            "F470N"	:_F(4.708,  0.160,2.535, NIRCAM,LONG),
            "F480M"	:_F(4.874,  0.162,2.574, NIRCAM,LONG),

            "F560W" :_F(5.589,  0.182,1.636, MIRI,NULL),
            "F770W" :_F(7.528,  0.243,2.187, MIRI,NULL),
            "F1000W":_F(9.883,  0.321,2.888, MIRI,NULL),
            "F1130W":_F(11.298, 0.368,3.318, MIRI,NULL),
            "F1280W":_F(12.712, 0.412,3.713, MIRI,NULL),
            "F1500W":_F(14.932, 0.483,4.354, MIRI,NULL),
            "F1800W":_F(17.875, 0.580,5.224, MIRI,NULL),
            "F2100W":_F(20.563, 0.665,5.989, MIRI,NULL),
            "F2550W":_F(25.147, 0.812,7.312, MIRI,NULL),
            }

##need to convert to DN or something

GAIN={ SHORT:(2.05,0.4), LONG:(1.82,0.4), NULL:(4,0)}
