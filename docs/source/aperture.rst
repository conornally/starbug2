*******************
Aperture Photometry
*******************

Aperture photometry is automatically excecuted at the end of the :doc:`detection routine <./detection>` but in the situation where you want to run photometry on a list of sources, the two routines have been separated. First a source list needs to be loaded into *starbug2*, this is done with the :code:`-d` or :code:`--apfile` flags in the command. The sourcelist must contain some kind of positional columns, either "RA DEC" or "xcentroid ycentroid" and be a binary fits table.
The apeture photometry routine is called with :code:`-A` or :code:`--apphot`::

    $~ starbug2 -d sourcelist.fits -A image.fits

The routine will create a table with the same name as the input image and append :code:`-ap.fits` to the end. Note this is the same as the detection routine discussed in the :doc:`Source Detection <./detection>` section and will contain the same columns. Note this might overwrite an original source list produced by the detecion routine. 

Relevant Parameters
-------------------

APPHOT_R
    Set the aperture radius in pixel units to be used. If this value is blank and no value for **ENCENRGY** is set, *starbug2* will use the **FWHM** value in place of **APPHOT_R**. If *starbug2* is run in verbose mode, the routine will report that this has happened.

ENCENERGY
    Mutually exclusive with **APPHOT_R**.

    Calculate an aperture radius to be used base on a fraction of enricled energy of the PSF (float between 0 and 1). To do this, *starbug2* requires an **APCORR_FILE** to be loaded, either automatically in the case of JWST images, or manually using the parameter below. 

    If no value is added here, *starbug2* will use **APPHOT_R** instead.

SKY_RIN
    Set the inner sky annulus in pixel units to be used the photometry routine. If this value is smaller than **APPHOT_R**, *starbug2* will used **APPHOT_R** +1 instead.

SKY_ROUT
    Set the outer sky annulus in pixel units to be used the photometry routine. If this value is smaller than or equal to **SKY_RIN**, *starbug2* will used **SKY_RIN** +1 instead.

APCORR_FILE
    Set the file containing an aperture correction table to be used to correct the flux. This file must be a binary fits table containing the following:

    - a column named "radius" in pixel units

    - a column called "apcorr" corresponding to these radii

    - a column called "eefraction" containing the fraction of encircled energy corresponding to these radii

    - an optional column called "filter", used to parse a particular filter from a table containing aperture corrections for  multiple photometry bands.

    If this value is left blank, *starbug2* will attempt to find one in the :code:`STABRUG_DATDIR` directory, if this fails, no aperture correction will be applied.

AP_FILE
    Set the source list to use for the aperture photometry.


A Typical Run
-------------

Running aperture photometry separately from source detection is something likely to be done if you are trialling different values for aperture radius and sky annulus sizes. The process can be drastically sped up by not constantly repeating the source detection routine. As such, the run may look like::

    $~ starbug2 -D image.fits

    $~ starbug2 -d image-ap.fits -sAPPHOT_R=2 -o image-r2.fits -A image.fits
    $~ starbug2 -d image-ap.fits -sAPPHOT_R=3 -o image-r3.fits -A image.fits
    $~ starbug2 -d image-ap.fits -sAPPHOT_R=4 -o image-r4.fits -A image.fits

Where a master source list called "image-ap.fits" is created in the first step, and then three separation aperture photometry runs are conducted. The original source list is loaded with :code:`-d image-ap.fits` , the aperture radius is altered with :code:`-s APPHOT_R=X` and the tables are output to new filenames with :code:`-o output`.
