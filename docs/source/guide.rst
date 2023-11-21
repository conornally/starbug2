**********
User Guide
**********

*StarbugII* has a simple command line interface that follows standard GNU Linux argument input and parsing standards.  There are two main programs associated with the tool: :code:`starbug2` the core photometry suite of routines and :code:`starbug2-match` a basic collection of matching routines.

.. toctree::
   :caption: Detailed Guides
   :maxdepth: 1
   :titlesonly:

   install.rst
   usage/param.rst
   usage/detection.rst
   usage/aperture.rst
   usage/background.rst
   usage/psf.rst
   usage/matching.rst
   usage/artificial.rst


Basic Usage
-----------


The two basic input and output products are photometric source lists and images. The various *starbug2* routines will take and/or export one or both of these products. The program uses standard fits formatted files and for most of the routines, the inputs can either have been created by *starbug2* at a previous stage or user defined. The underlying ethos of *starbug2* is that the tool is simple and powerful and allows the user to do whatever they want, it is however up to their discretion, whether what they are doing is a good idea.
For example, one may conduct forced photometry from *HST* source lists on *JWST* images, or run faint source detection on residual images.


To see the basic usage information for *starbug2*, run::

    $~ starbug2 -vh    

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

               --apply-zeropint    a.fits : Apply a zeropoint (-s ZP_MAG=1.0) to a.fits
               --calc-instr-zp     a.fits : Calculate and apply an instrumental zero point onto a.fits

           --> typical runs
              $~ starbug2 -vD -p file.param image.fits      //Source detect on image with a parameter file
              $~ starbug2 -vDM -n4 images*.fits             //Source detect and match outputs of a list of images
              $~ starbug2 -vd image-ap.fits -BP image.fits  //PSF photometry on an image with a source file (image-ap.fits)


For basic usage information for *starbug2-match*, run::

    $~ starbug2-match -vh

        usage: starbug2-match [-BGfhv] [-o output] [-p file.param] [-s KEY=VAL] table.fits ...
            -B  --band               : match in "BAND" mode (does not preserve a column for every frame)
            -C  --cascade            : match in "CASCADE" mode (left justify columns)
            -D  --dither             : match in "DITHER" mode (preserves a column for every frame)
            -f  --full               : export full catalogue
            -G  --generic            : match in "GENERIC" mode
            -h  --help               : show help message
            -o  --output  file.fits  : output matched catalogue
            -p  --param   file.param : load starbug parameter file
            -s  --set     option     : set value in parameter file at runtime (-s MATCH_THRESH=1)

            --> typical runs
               $~ starbug2-match -Gfo outfile.fits tab1.fits tab2.fits
               $~ starbug2-match -sMATCH_THRESH=0.2 -sBRIDGE_COL=F444W -Bo out.fits F*W.fits

