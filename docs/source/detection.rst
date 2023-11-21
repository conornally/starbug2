****************
Source Detection
****************

The first routine that is likely to be run with *starbug2* is the source ddetection. This takes an image and identifies the location and brightness of the point present point-like sources. The routine is called with :code:`--detect` or :code:`-D` for short and including an appropriate fits image file::
    
    $~ starbug2 -D image.fits

This will create a binary fits table called "image-ap.fits" in the current working directory. This table will contain the following columns:

    +------------------+------------------------------------------------------------+
    | NAME             | Description                                                |
    +------------------+------------------------------------------------------------+
    | Catalogue Number | Source index in catalogue                                  |
    +------------------+------------------------------------------------------------+
    | RA/DEC           | World coordinates of each source centroid                  |
    +------------------+------------------------------------------------------------+
    | x/ycentroid      | Pixel coordinates of each source centroid                  |
    +------------------+------------------------------------------------------------+
    | sharpness        | Source sharpness parameter                                 |
    +------------------+------------------------------------------------------------+
    | roundness1/2     | Source roundness parameters                                |
    +------------------+------------------------------------------------------------+
    | smoothness       | Source smoothness parameter                                |
    +------------------+------------------------------------------------------------+
    | flux/eflux/sky   | Flux, photometric error and sky value of source            |
    +------------------+------------------------------------------------------------+
    | FILTER/eFILTER   | Magnitude and magnitude error of source                    |
    +------------------+------------------------------------------------------------+
    | flag             | Source quality flag                                        |
    +------------------+------------------------------------------------------------+


At the end of the detection routine, aperture photometry is automatically executed to measure the flux of the sources. See :doc:`Aperture Photometry <./aperture>` for more details.

To output the source list to a different filename or folder, include :code:`-o outputfile.fits` or :code:`-o path/to/folder/` in the *starbug2* command.

.. tip::

    As with all *starbug2* routines, calling **DETECT** in verbose mode with the :code:`-v` flag will allow you to see the progression of the code as well as any useful outputs or warnings that have occurred. 


To inspect the output source list, you can open it with a fits viewer like `Topcat <https://www.star.bris.ac.uk/~mbt/topcat/>`_ or astropy `Tables <https://docs.astropy.org/en/stable/table/>`_. Alternatively it can be converted to a DS9 region file with::
    
    $~ starbug2 --generate-region image-ap.fits

You can customise the output region file by setting ...

Relevant Parameters
-------------------

The default *starbug2* detection parameters have been designed to work well with most images, however you may find it needs to be honed for your particular image or filter. Typically it will be worth creating a local parameter file with :code:`$~starbug2 --local-param`, see the :doc:`parameters <./param>` section for more information. Discussed below are some of the relevant parameters to hone the detection.

FWHM
    When running *starbug2* on "standard" JWST images, it will automatically detect the **FWHM** for that given **FILTER** , however in the case where it doesn;t manage (either the header file of the image is non-standard or the image is not from JWST) **FWHM** should be included explicitly in pixel units.

SIG_SRC
    You may find that clean and unambiguous sources may be being missed by the routine, this might be due to the source flux being below the detection threshold. Reduce this float to set the number of sigma above the background flux that a source must have to be detected. A robust detection is generally considered to be 5.0 sigma above the sky, but in some cases it may likely be reduced to 3.0. Setting **SIG_SRC** below this will begin to increase the number of spurious detections.

SIG_SKY
    Sigma clipping occurs on the image during the detection routine. **SIG_SKY** sets the number of sigma above the image median below which gets clipped. This reduces the number of spurious detections that may occur. 
    
    However, if the image contains very bright or complex diffuse emission from dust, there are likely obscured stars that are fainter than the brightest "dusty peaks". In this situation, reducing this parameter is a very productive way of increasing the detection rate of faint sources. It is recommended to reduce **SIG_SKY** slowly, in around 0.1 steps and carefully inspect the output as this will also increase the number of "dusty peaks" visible to the *starbug2*. Careful setting of the following parameters are required to mitigate this effect. 

    It would usually not be recommended to set **SIG_SKY** below 1.0 unless being done with a lot of care and precision. 

SHARP_LO / SHARP_HI
    :code:`Sharpness` is the measurement of the ratio the peak pixel in the source to the mean of the surrounding pixels. In other words, a sharp source has a high centre with narrow wings. Artifacts such as cosmic rays or bright pixels are likely to have high sharpness values whereas "dusty peaks" or resolved background galaxies may be less sharp.

    **SHARP_LO** and **SHARP_HI** set acceptable bounds of detected sources. Inspect the distribution of :code:`sharpness` values and widen or reduce the bounds as necessary to reduce the detection rate of non point-like sources. Note that this distribution is dependant on the PSF of the image and bounds for one image may not be appropriate for another.

SMOOTHNESS
    :code:`Smoothness` like :code:`sharpness` is a measure of how "pointy" a source is. It takes the ratio of mean pixel values as measured in two apertures around the source. A very smooth detection, such as a dusty peak or resolved background galaxy will have a value around 1.0, whereas a "good" star will have a lower value.
    **SMOOTH_LO** and **SMOOTH_HI** set the acceptable bounds of detected sources.

    This parameter is designed to do a similar job as sharpness but from the other direction. It is very effective at mitigating bright dust or "null" detections sometimes seen in empty areas of an image. As it relies on aperture photometry to measure, it is affected by crowding in really dense regions, and "smooth" sources may in fact be close optical binaries.


.. note::

    :code:`smoothness` is currently an experimental parameter and the exact definition may change in the future.

ROUNDNESS
    :code:`Roundness` is a measure of source eccentricity. There are two versions of this metric. :code:`roundness1` describes the 4-fold symmetry of a source and :code:`roundness2` is a ratio of two fitted 1D gaussians to the source, one vertical and horizontal. Both values are symmetric distributions centred on zero. **ROUND1_HI** and **ROUND2_HI** set the outer limits for their respective distributions. 

    Highly eccentric sources have roundness values further from zero. These are usually PSF fringes or resolved background galaxies. Inspecting the two :code:`roundness` distributions often reveals an underlying normal-like distribution with wings, these wings can be clipped to leave the cleaner point-like sources.

    Fundamentally, both :code:`roundness` values measure similar things but often they trace slightly different distributions and can be tweaked independently to remove outlying sources. 

RICKER_R
    This parameter sets the radius in pixel units of the wavelet convolved with the image during the *CONVL* stage of the detection routine. In noisy images, small values of **RICKER_R** can over detect spurious sources. In this case, try increasing the number to ~5/10 and then decreasing it in integer steps, while inspecting the result.

A Typical Run
-------------






Introducing Dithers
-------------------

All the *starbug2* routines can run on more than one image at a time. You can set this up by adding more images to the end of the command::
    
    $~ starbug2 -vD image1.fits image2.fits image3.fits image4.fits ...

Or make use of bash wildcarding and string formatting::

    $~ starbug2 -vD image*.fits

    $~ starbug2 -vD image{1..4}.fits

By default *starbug2* will execute the routine on these images in series but the process can be parallelised with addition of the :code:`-n` or :code:`--ncores` flags which take an integer option for the number of cores to parallelise over. If the number of cores is smaller than the number of images supplied, the extra images will get added to a queue and wait till the current process is complete. If the number of cores is greater than the number of images, there will be no ill side-effects. Finally there is no limit to the number of images that *starbug2* can run on, other than the available resources on the computer.
An example of this would look like::

    $~ starbug2 -vD -n4 image*.fits

.. only:: html

    .. warning::
        :class: sphx-glr-download-link-warning

        There is a current issue where a child process runs out of memory and is killed but this is not reported back to the main process, which is waiting for the child to finish. The result of this, is the main process getting stuck waiting and never completes. This all happens silently so you would need to watch the number of cores being actively used on the system. 

        Reduce the number of cores being used if this happens.

The source lists produced from each image can be matched together at the end of the detection run by including the :code:`-M` or :code:`--match` flag to the command. More information on this is available in the :doc:`Catalogue Matching <./matching>` section. An example of this command would be::

    $~ starbug2 -vDM image1.fits image2.fits 
