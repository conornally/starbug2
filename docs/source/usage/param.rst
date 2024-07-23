**************
Parameter File
**************

The parameter file is where any dataset specific parameters can be tweaked. Ideally the default values should be sufficient however if the diffuse dust emissions cause a very complex background or the field is very crowded, certain parameters may need tuned. Additionally differing sensitivities between photometric bands may require different parameters to detect on.
To generate a local parameter file, run::

    $~ starbug2 --local-param

This will create :code:`./starbug.param`, a file which will be loaded by default when *starbug2* is ran from the folder that contains it. However it can be named anything and explicitly loaded at runtime with the :code:`-p paramfile` OR :code:`--param paramfile` like::

    $~ starbug2 -p file.param  ...

It may be the case that you keep several parameter files on hand to conduct slightly different routines, this is how you would load them.

For quick testing or tweaking of certain parameters, *starbug2* can override a file setting with the addition of :code:`-s PARAM=VALUE` OR :code:`--set PARAM=VALUE` in the command. This will not change the parameter file but will prioritise this value at runtime. This feature may be useful for quick testing of different detection threshold cuts for example, without losing your current "safe" settings. This will look something like::

    $~ starbug2 -p file.param -s SHARP_HI=0.8 ...


If your parameter file was generated in an older version of *starbug2*, new parameters may have been added or old ones deprecated. If this happens, *starbug2* will abort the run with a warning and ask you to update the file. This can be done automatically with::
    
    $~ starbug2 --update-param


Settable Parameters
-------------------

As of the current version of *starbug2*.

These parameters will be discussed in more detail in their relevant sequences.

VERBOSE [INT 0:1]
    Run *starbug2* in verbose mode.

OUTPUT [STR]
    Output file or directory name. If using a specific filename, *starbug2* will append the relevant suffixes to the data products, i.e a detection list will go to --> OUTPUT-ap.fits 

HDUNAME [INT/STR]
    *starbug2* will look for "SCI" "BGD" "RES" fits extension names, if these are not present it will use the first extension in the file (HDUNAME=0). If *starbug2* is being run on a non standard fits image, set the index or extension name with this parameter.

FILTER [STR]
    If the keyword "FILTER" is not present in the fits image header file, *starbug2* will use "mag" as the filter name. This may have unintended consequences in for example the matching. Set the filter name explicitly with this parameter.

FWHM [FLOAT>0]
    If *starbug2* can't find the full-width half maximum of the image, set it explicitly here. This is required for many steps in the detection and photometry routines. When running *starbug2* on standard JWST images, it should be able to automatically set the FWHM, if it doesn't, it will warn you. Set the FWHM in pixel units with this parameter.

SIGSKY [FLOAT>0] default 3.0
    Number of sigma below the image median which pixels gets removed as sky. In images with bright diffuse emissions, drop this value slowly, in steps of ~0.1 and watch the number of sources detected increase, but be careful of false detections of bright spots.

SIGSRC [FLOAT>0] default:5.0
    Minimum number of sigma above the median that a source must be to be detected. This is often 5 sigma for a robust search or 3 sigma for faint detections. It is usually not encouraged to go below 3 sigma.

DOBGD2D [INT 0:1] default:1
    Do the BGD2D detection step. This increases the execution time but may find new sources.
    
DOCONVL [INT 0:1] default:1
    Do the CONVL detection step. This increases the execution time but may find new sources.

CLEANSRC [INT 0:1] default:1
    Do source cleaning at the end of the detection routine. This removes sources that lie outwith the geometric quality parameters (defined below), to remove bad detections and background galaxies. Turn this off with '0' but expect the catalogue to contain a lot of extra bad sources.

SHARP_LO [FLOAT] default:0.4
    Set the lower bound for how "sharp" a point source is.

SHARP_HI [FLOAT] default:0.9
    Set the upper bound for how "sharp" a point source is.

ROUND1_HI [FLOAT] default:1.0
    ROUNDNESS1 is a measure of symmetry in the source. It is a symmetric distribution centred on zero. Set the magnitude of the outer limit with this parameter.

ROUND2_HI [FLOAT] default:1.0
    ROUNDNESS2 is a measure of ratio between a 1D gaussian fit horizontally and vertically to the source. It is a symmetric distribution centred on zero. Set the magnitude of the outer limit with this parameter.

SMOOTH_LO [FLOAT] 
    Set the lower bound for how "smooth" a source is. This distribution is unbounded but clean sources normally shouldn't go below 0.

SMOOTH_HI [FLOAT] default:1
    Set the upper bound for how "smooth" a source is. This parameter should be tweaked very slowly and the results investigated in detail. It is very effective at removing spurious detections in dusty regions but may remove "good" sources in crowded regions.

RICKER_R [FLOAT>0] default:1
    Set the radius for the wavelet convolved with the image during the CONVL routine. In MIRI images, this will likely need to be increased to limit spurious detections in the dust structures. Set it high and bring it down slowly to see the effect.

APPHOT_R [FLOAT>0] default:1.5
    Aperture radius in pixel units. If left blank or set <0, *starbug2* will use **ENCENRGY** to calculate the aperture radius.

ENCENERGY [FLOAT 0-1] default:-1
    Calculate the aperture radius from the "percentage encircled energy". This requires an aperture correction file (APCORR_FILE) to be either explicitly set, or automatically loaded in the case of JWST images. **APPHOT_R** takes precedence over this parameter, but if it fails, it will use the value for **FWHM**.

SKY_RIN [FLOAT>0] default:3
    Set the inner radius for the sky annulus for aperture photometry in pixel units. This should be greater than **APPHOT_R**, in the case it is not, *starbug2* will automatically set it to **APPHOT_R** +1.

SKY_ROUT [FLOAT>0] default:4.5
    Set the outer radius for the sky annulus for aperture photometry in pixel units. This should be greater than **SKY_RIN**, in the case it is not, *starbug2* will automatically set it to **SKY_RIN** +1.

APCORR_FILE [STR]
    Set the filename for an explicit aperture correction file. This file should contain the columns "radius" "apcorr" and optionally "eefraction". This will be used to calculate the aperture correction for a given aperture radius.

BGD_R [FLOAT>0]
    Set a aperture source masking radius in pixel units to be used during the diffuse background estimation routine. By default *starbug2* will try to calculate good values for each source, but in cases where this doesnt act appropriately, a uniform radius can be set with this parameter.

PROF_SCALE [FLOAT>0] default:1.0
    The aperture masking in the background estimation steps scales the aperture radii with source flux. This parameter is used to change the scale factor *A* of the profile.

PROF_SLOPE [FLOAT>0] default 0.5
    The aperture masking in the background estimation steps scales the aperture radii with source flux. This parameter is used to change the slope exponent *B* of the profile.
    
BGD_CHECKFILE [STR]
    The scaled aperture masking is done under the hood of *starbug2*, set this value to a filename to output a ds9 region file containing the calculated aperture radii for each source.

BOX_SIZE [FLOAT>1] default:2.0
    Set the kernel size in pixel units to be used with estimating the background emission. When set small, the resulting model will contain a lot of detail at a high resolution but may be influenced by extraneous bright pixel or undetected sources. When set large, the resulting model will be low resolution but will be less influenced by bright pixels.

AP_FILE [STR]
    Set the filename for a source list to be explicitly used during any *starbug2* runs. The file must contain columns xcentroid,ycentroid or RA,DEC and the image must contain WCS information to convert that to pixel coordinates.

BGD_FILE [STR]
    Set the filename for a background estimation file to be loaded into any *starbug2* runs. This must be a fits image with the same dimensions as the image being actively worked on.

PSF_FILE [STR]
    Set the filename for a psf file to be loaded into any *starbug2* runs. This must be a fits image.

USE_WCS [INT 0:1] default:1
    If a loaded **AP_FILE** contains both RADEC and xycentroid columns, do you want *starbug2* to use the RADEC columns true (1) or false (0).

ZP_MAG [FLOAT] default:8.9
    Set a zeropoint to be added to the magnitudes. By default it is 8.9 for AB magnitudes.

CRIT_SEP [FLOAT>0] default 8.0
    Deprecated

FORCE_POS [INT 0:1] default:0
    Set whether you want the conduct "forced centre" PSF photometry.

DPOS_THRESH [FLOAT>0] 
    Deprecated

MAX_XYDEV [FLOAT>0] default:3p
    Set the maximum deviation of the central position of a source from its initial guess during PSF photometry. Sources that deviate by more than this value will be refit with forced centres. This value by default is in unit pixels but can be set in other units by appending (p:pixels, s:arcseconds, m:arcminutes, d:degrees) to the end of the values.

PSF_SIZE [INT>0]
    Explicitly set the size of the PSF to be used in pixel units. If not set, *starbug2* will use the whole psf image.

GEN_RESIDUAL [INT 0:1] default:0
    Generate a residual image at the end of PSF photometry. This will have the estimated background and the fitted sources subtracted.

CALC_CROWD [INT 0:1]
    Calculate the crowding parameter during source geometry routine. Scales exponentially with the number of sources in the list.


MATCH_THRESH [FLOAT>0] default:0.1
    The maximum separation of two sources matching together in arcsecond units.

MATCH_COLS [STR] 
    Extra columns to include in matching. Space separated list i.e. "sharpness smoothness".

NEXP_THRESH [INT>1]
    Threshold for the minimum number of detections a source must have while combining catalogues. If a source is only detected in two of four exposures and **NEXP_THRESH** =3, these sources will not make it into the finalised catalogue.

SN_THRESH [FLOAT]
    in prep.

BRIDGE_COL [STR]
    in prep.

NUMBER_ARTIFICIAL_STARS [INT>0]

SUBIMAGE_SIZE [INT>0]

MIN_FLUX [FLOAT>0]

MAX_FLUX [FLOAT>0]

SEPARATION_THRESH [FLOAT>0]

