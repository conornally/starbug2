# StarBugII Manual

JWST PSF photometry in dusty crowded fields.
Last updated: v0.2.1

## Installation

```bash
pip install starbug2-X.X.X.tar.gz

-- OR --

pip install --user starbug2-X.X.X.tar.gz 
```

After the package is installed, there are a few steps required to initialise Starbug.

**WEBBPSF** Is a dependency of Starbug that has its own initialisation process. The full installation is documented on https://webbpsf.readthedocs.io/en/latest/installation.html however it requires two main steps. Download the data file on the website, named something like webbpsf-data-X.X.X.tar.gz and expand it into a directory, then add append to your .bashrc (or equivalent) `export "WEBBPSF=PATH/TO/DIRECTORY"`.

**PSFDIR** This is the folder where starbug stores its relevant data files. By default this is "${HOME}/.local/share/starbug". Make sure this folder exists, or if you wish to save them elsewhere, change the folder (permanently) in "starbug2/starbug2/param/default.param [PSFDIR]", or (temporarily) in a local starbug.param file.

**PSF FILES** Starbug requires PSF files to be generated for the filters you are using. To do so, run `starbug2 --generate-psf` and they will be generated into "PSFDIR"

## Usage

```bash
Starbug II - JWST PSF photometry
usage: starbug2 [-ABCDfhMPv] [-b bgdfile] [-d apfile] [-o directory] [-p file.param] [-s opt=val] image.fits ...
   -A  --artific         : run artificial star tests
   -B  --background      : run background estimation
   -b  --bgdfile         : load background (-bgd.fits) file
   -C  --clean           : run source cleaning before photometry 
   -d  --apfile  ap.fits : load a source detection (-ap.fits) file to skip the source detection step
   -D  --detect          : run source detection
   -f  --find            : attempt to find associated -ap -bgd files
   -h  --help            : display uasage information
   -M  --match           : match outputs from all input image files
   -n  --ncores      num : number of CPU cores to split process between
   -o  --output      dir : output directory
   -p  --param   a.param : load parameter file
   -P  --photom          : run psf photometry
   -s  --set      option : set value in parameter file at runtime (-s SIGSKY=3)
   -v  --verbose         : display verbose outputs

   --> Single run commands
       --generate-psf             : Generate ALL the PSF files to "PSFDIR"
       --local-param              : Make a local copy of the default parameter file
       --generate-region   a.fits : Make a ds9 region file with a detection file
       --clean-table       a.fits : Clean up an individual table
       --generate-run      *.fits : Generate a simple run script
       --version                  : Print starbug2 version
```

## Parameter File

The Parameter file is where any dataset specific parameters can be tweaked. Ideally the default values should be sufficient however if the diffuse dust emissions cause a very complex background or the field is very crowded, certain parameters may need tuned.

### General Parameters 

**VERBOSE** (0:false, 1:true) Set whether there will be logging outputs throughout the execution

[NULLVAL=999.999] Currently not being used

**PSFDIR** Set the directory where the PSF files are stored. By default this is `${HOME}/.local/share/starbug`. Note this path will expand any environment variables (${HOME} --> /home/dlister)

### Detection Parameters

**SIGSKY** (float>0) Set the number of sigma below the median pixel value which will get clipped out during background subtraction as "definitely sky". Note, for very bright dust emissions the median value may be comparable to some of the faintest sources, setting this value very low will allow the recovery of these sources at the detriment to picking up some bright dust regions.

**SIGSRC** (float>0) The minimum number of sigma that a source must be above the median pixel value. Set this high for only the brightest sources or low to detect the faint ones. Be careful setting this below 3sigma unless you have a good reason too.

**BOX_SIZE** (int>0) Kernel size in pixels during background2D subtraction. For complex dusty regions this should be set low.

**FILTER_SIZE** Emmm...look this up?

**MATCH_THRESH** (float>0 arcsec) Source catalogues from all the background subtraction methods are matched together to create a single source list. If a source has a separation larger than MATCH_THRESH (arcsec) it will be considered a "new source" that isnt already in the source list.

**SHARP_LO/SHARP_HI** (float) Set the bounds of sharpness, outside which the source detector will discard a potential source.

**ROUND_LO/ROUND_HI** (float) Same effect as SHARP_LO/SHARP_HI but for the source roundness (symmetrical around zero).

**APPHOT_R** (float>0) Set the aperture radius for aperture photometry following source detection.

**APPHOT_R0/APPHOT_R1** (float>0) Set the inner (R0) and outer (R1) radii for sky calculation during aperture photometry.

### Catalogue Cleaning

**ERROR_CUT** (float>0) Cut sources with a magnitude error greater than this value

**SHARP_HI_SIG/SHARP_LO_SIG** (float>0) The cleaning routine fits a gaussian distribution to the catalogue sharpness values. Cut values n sigma above and below the mean.

**ROUND_HI_SIG/ROUND_LO_SIG** (float>0) The same as SHARP_HI_SIG/SHARP_LO_SIG but for roundness values

### PSF Photometry

**AP_FILE** (file path) Set a source detection aperture photometry file (apfile) to be used to during photometry. If you are running photometry on a different starbug call than source detection then the apfile must be set either this way or with `starbug -d filename ...` 

**CRIT_SEP** [Cant remember]

**PSF_ITERATIONS** [Not yet implemented]

### Artificial Star Testing

**NUMBER_ARTIFICIAL_STARS** (int>0) Set the number of tests to run.

**SUBIMAGE_SIZE** (int>0) Starbug crops the image into smaller sub images for every inserted star. Set the (square) size in pixels of the sub images. Note, the larger this is the longer the tests will take, however setting it too small will likely effect the accuracy of the test.

**MIN_FLUX/MAX_FLUX** (float>0) Set the minimum and maximum inserted source flux. This should be just smaller than the minimum and larger than the maximum flux of the stars in the psf catalogue.

**SEPARATION_THRESH** (float>0 pix) The separation between the input and recovered source must be below this threshold value, otherwise the test fails that source. Large values are more lenient but less accurate. However the PSF size varies a lot with the photometric band and instrument, so this will need tweaked with every different filter.



## Detection

## Photometry

## CleanUp

## Artificial Star Testing

# Notes to Me

Source detection 1,2,3 are basically the same, i should replace them with other methods
