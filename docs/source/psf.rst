**************
PSF Photometry
**************

The PSF photometry routine requires a source list to be generated with either the :doc:`detection routine <./detection>` or loaded with :code:`-d sourcelist.fits`. Additionally, an optional background estimation image can be generated with the :doc:`background estimation routine <./background>` or loaded with :code:`-b image-bgd.fits`, this will be subtracted from the raw image before the photometry is conducted.

.. tip::

    You can automatically load the corresponding source list and background file for a given image by adding the :code:`-f` or :code:`--find` flag to the command. 
    This will look for files with the same base filename but search for -ap and -bgd suffixed. For example :code:`$~ starbug2 -f image.fits` will locate "image-ap.fits" and "image-bgd.fits" if they exist in the working directory. 

    Note, you can also still load a global sourcelist which will take priority over the automatic one. For example :code:`$~ starbug2 -d sourcelist.fits -f image.fits` will load "sourcelist.fits" and "image-bgd.fits". This is useful in the situation where you are running PSF photometry on a large list of images simultaneously and don't want to manually include the appropriate background file every time but still want to use a carefully constructed sourcelist.

The PSF fitting routine is run with the :code:`-P` or :code:`--psf` flags. A simple run may look like either of the commands below. The first detects on the image and estimates the background, then conducts the PSF routine with these product. The second loads an existing source list and background image::
    
    $~ starbug2 -vDBP image.fits

    $~ starbug2 -vd sourcelist.fits -b background.fits -P image.fits

The routine will produce a binary fits table with the following columns:

   +------------------+--------------------------------------------------+
   | Name             | Description                                      |
   +------------------+--------------------------------------------------+
   | Catalogue_Number | Source index in table                            |
   +------------------+--------------------------------------------------+
   | x/y_init         | Initial centroid position of each source         |
   +------------------+--------------------------------------------------+
   | x/y_fit          | Fitted PSF centre for each source                |
   +------------------+--------------------------------------------------+
   | RA/DEC           | WCS measured from x/y_fit                        |
   +------------------+--------------------------------------------------+
   | xydev            | Distance between x/y_init and x/y_fit            |
   +------------------+--------------------------------------------------+
   | ap_FILTER        | Initial guess for photometry                     |
   +------------------+--------------------------------------------------+
   | flux/eflux       | Fitted flux and photometric error                |
   +------------------+--------------------------------------------------+
   | FILTER/eFILTER   | Magnitude and photometric error of each source   |
   +------------------+--------------------------------------------------+
   | qfit             | Measure of PSF fit quality (in prep.)            |
   +------------------+--------------------------------------------------+
   | flag             | Source quality flags                             |
   +------------------+--------------------------------------------------+




Relevant Parameters
-------------------

CRIT_SEP
    .

FORCE_POS
    The PSF routine can be run in two different modes, "Forced" or "Unforced". In Forced photometry, the central point of the point source is fixed at the input position and the photometry only fits for flux. This is much faster but can result in poorer flux measurements. Unforced photometry allows the central position to be fit as well as the flux, this is much more accurate but there are three times as many free parameters to fit. Set whether to conduct the routine in forced mode **FORCE_POS** = 0 or 1 (false or true).

MAX_XYDEV
    In "Unforced" mode, the central position of the point source can wander far from the initial position, normally in the case of a non point-like object at the centre. Setting **MAX_XYDEV** will set the threshold maximum deviation from the initial position allowed. Source that move further will be have their fluxes refit at the original position in forced mode. 

    **MAX_XYDEV** is in pixel units by default, but can be set in other units by append 'p': pixel, 's':arcsecond, 'm':arcminute, 'd':degree to the number. For example **MAX_XYDEV** =0.3s for 0.3 arcseconds.

PSF_SIZE
    Set the size of the PSF to fit, in pixel units. By default it will use the whole **PSF_FILE** image, but if you want to reduce the size of the PSF to speed up the process, set it here.

GEN_RESIDUAL
    Generate a residual image, this is the raw image with a background and source PSFs subtracted. It will be exported to "filename-res.fits". This is useful to inspect the quality of the PSF fits and may reveal hidden faint sources which weren't detected. You may wish to run a second round of photometry on this residual.

PSF_FILE 
    Load a specific fits image to be used as a PSF.

AP_FILE and BGD_FILE
    Load a source list and or a background image.


A Typical Run
-------------
