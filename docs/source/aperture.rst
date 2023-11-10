*******************
Aperture Photometry
*******************

Aperture photometry is automatically excecuted at the end of the :doc:`detection routine <./detection>` but in the situation where you want to run photometry on a list of sources, the two routines have been separated. First a source list needs to be loaded into *starbug2*, this is done with the :code:`-d` or :code:`--apfile` flags in the command. The sourcelist must contain some kind of positional columns, either "RA DEC" or "xcentroid ycentroid" and be a binary fits table.
The apeture photometry routine is called with :code:`-A` or :code:`--apphot`::

    $~ starbug2 -d sourcelist.fits -A image.fits

Relevant Parameters
-------------------

APPHOT_R

ENCENERGY

SKY_RIN

SKY_ROUT

APCORR_FILE

AP_FILE
