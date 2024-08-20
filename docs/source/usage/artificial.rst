***********************
Artificial Star Testing
***********************

Artificial Star Testing (AST) is the injection of known model PSFs into an image and determining the quality of the detection and flux recovery of the routine. *Starbug2* has a separate executable called :code:`starbug2-ast` to conduct artificial star testing.

The AST routine operates on an image and injects *NSTARS* with known brightnesses from the range *MAG_MIN* - *MAG_MAX* onto the array, runs the :doc:`detection <./detection>`, :doc:`background estimation <./background>` and :doc:`PSF photometry <./psf>` routine, and reports the resulting measured flux (if detected). This process runs *NTESTS* times, building a large catalogue of detection rate and uncertainty as a function of input brightness.

Output
------

The output is in the form as a multi extension fits file, with the following extensions:

AST : TABLE
    The first extension is the compiled user facing results. This is a table with detections rate "*rec*" and error "*err*" and a function of input magnitude "*mag*". The error column is calculated as the standard deviation of the flux measurements within the magnitude bin. 

RAW : TABLE
    The raw injected source data.

CMP : IMAGE (in prep.)
    An image with the same dimensions as the working input image, with the detection rate fraction as a function of pixel coordinate.


If the parameter *PLOTAST* is specificed, the AST table will be compiled into a readable figure, showing the detection rate decay as a function of magnitude and overlaying the 70% and 50% completeness magnitude. The parameterised s-curve coefficients are given in the AST header file.

.. image:: ../_static/images/ast.pdf
   :width: 325
   :alt: Artificial Star Test output image
    

Usage
-----

Usage of :code:`starbug2-ast` can be shown with::

    $~ starbug2-ast -vh

    StarbugII Artificial Star Testing
    usage: starbug2-ast [-vh] [-N ntests] [-n ncores] [-p file.param] [-S nstars] [-s opt=val] image.fits
        -h  --help          : show help screen
        -N  --ntests    num : number of tests to run
        -n  --ncores  cores : number of cores to split the tests over
        -o  --output output : output directory or filename to export results to
        -p  --param    file : load a parameter file
        -S  --nstars    num : number of stars to inject per test
        -s  --set    option : set parameter at runtime with syntax "-s KEY=VALUE"
        -v  --verbose       : show verbose stdout output


Relevant Parameters
-------------------

NTESTS : 
    .

NSTARS : 
    .

SUBIMAGE : 
    .

MAX_MAG : 
    .

MIN_MAG : 
    .

PLOTAST : 
    .



