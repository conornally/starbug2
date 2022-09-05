# Starbug II Manual

Last edited during starbug version v0.1.6

## Install 

```bash
pip install starbug2-X.X.X.tar.gz
```

One of the required packages is STSCIs "webbpsf". This package is used to generate the up to date PSF files. It has its own install requirements which can be found at: https://webbpsf.readthedocs.io/en/stable/installation.html#data-install. Essentially some package data must be downloaded and the location of their storage must be set as an environment variable "WEBBPSF_PATH"

\vspace{0.5cm}\hrule 

## First Run

### PSF Files

The PSF files must be generated to run any PSF fitting routines. These files are going to be stored in the folder specified in the parameter file, by default that is PSFDIR=~/.local/share/starbug but can be set to anything. The generation process will take several minutes

To generate the PSFs run:

```bash
starbug2 --generate-psfs 
```
\vspace{0.5cm}\hrule 

### Local Parameter File

It is advisable to generate a local parameter file for each project, so parameters can be set specifically to each dataset and not effect the global defaults. To do this run:

```bash
starbug2 --local-param
```
\vspace{0.5cm}\hrule 

### BASH TAB Completeion

A bash tab completion script has been included but not install by default. To install it:
```bash
echo "source PATH/TO/INSTALL/bin/complete.sh" >> ~/.bashrc
```

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
```

\newpage

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
