***************
Using StarbugII
***************

Getting Started
---------------






Basic usage for *starbug2* is viewed by running::

    $~ starbug2 -vh

    Starbug II - JWST PSF photometry
    usage: starbug2 [-ABDfGhMPSv] [-b bgdfile] [-d apfile] [-n ncores] [-o ouput] [-p file.param] [-s opt=val] image.fits ...
       -A  --apphot          : run aperture photometry on a source list
       -B  --background      : run background estimation
       -b  --bgdfile         : load background (-bgd.fits) file
       -d  --apfile  ap.fits : load a source detection (-ap.fits) file to skip the source detection step
       -D  --detect          : run source detection
       -f  --find            : attempt to find associated -ap -bgd files
       -G  --geom            : calculate geometric stats on source list
       -h  --help            : display uasage information
       -M  --match           : match outputs from all input image files
       -n  --ncores      num : number of CPU cores to split process between
       -o  --output      dir : output directory
       -p  --param   a.param : load parameter file
       -P  --psf             : run psf photometry
       -s  --set      option : set value in parameter file at runtime (-s SIGSKY=3)
       -S  --subbgd          : subtract background from image
       -v  --verbose         : display verbose outputs

       --> Single run commands
           --init                     : Initialise Starbug (post install)
           --local-param              : Make a local copy of the default parameter file
           --update-param             : Update an out-of-date local parameter file
           --generate-psf             : Generate a single PSF. Set FILTER, DET_NAME, PSF_SIZE with -s
           --generate-region   a.fits : Make a ds9 region file with a detection file
           --generate-run      *.fits : Generate a simple run script
           --version                  : Print starbug2 version

           --apply-zeropint    a.fits : Apply a zeropoint (-s ZEROPOINT=1.0) to a.fits
           --calc-instr-zp     a.fits : Calculate and apply an instrumental zero point onto a.fits

       --> typical runs
          $~ starbug2 -vD -p file.param image.fits      //Source detect on image with a parameter file
          $~ starbug2 -vDM -n4 images*.fits             //Source detect and match outputs of a list of images
          $~ starbug2 -vd image-ap.fits -BP image.fits  //PSF photometry on an image with a source file (image-ap.fits)


.. toctree::
   :caption: Detailed Guides
   :maxdepth: 1
   :titlesonly:

   param
   detection
   aperture
   background
   psf
   matching


