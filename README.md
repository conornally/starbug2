# StarBugII

JWST PSF photometry in dusty crowded fields.
Last updated: v0.2.14

**NOTE** This is still under heavy development, PSFs from JWST have not been officially released and as such, only aperture photometry can be used robustly.

[![Python application](https://github.com/conornally/starbug2/actions/workflows/python-app.yml/badge.svg)](https://github.com/conornally/starbug2/actions/workflows/python-app.yml)
[![PyPI version fury.io](https://badge.fury.io/py/starbug2.svg)](https://pypi.python.org/pypi/starbug2/)
[![Latest release](https://badgen.net/github/release/conornally/starbug2)](https://github.com/conornally/starbug2/releases)


## Installation

```bash
pip install starbug2

--- OR ---

git clone https://github.com/conornally/starbug2.git
python -m build
pip install .
```

After the package is installed, there are a few steps required to initialise Starbug.

**WEBBPSF** Is a dependency of Starbug that has its own initialisation process. The full installation is documented on https://webbpsf.readthedocs.io/en/latest/installation.html however it requires two main steps. Download the data file on the website, named something like webbpsf-data-X.X.X.tar.gz and expand it into a directory, then add append to your .bashrc (or equivalent) `export "WEBBPSF=PATH/TO/DIRECTORY"`.

**PSFDIR** This is the folder where starbug stores its relevant data files. By default this is "${HOME}/.local/share/starbug". Make sure this folder exists, or if you wish to save them elsewhere, change the folder (permanently) in "starbug2/starbug2/param/default.param [PSFDIR]", or (temporarily) in a local starbug.param file.

**PSF FILES** Starbug requires PSF files to be generated for the filters you are using. To do so, run `starbug2 --generate-psf` and they will be generated into "PSFDIR"

## Usage

```bash
Starbug II - JWST PSF photometry
usage: starbug2 [-ABCDfhMPv] [-b bgdfile] [-d apfile] [-o directory] [-p file.param] [-s opt=val] image.fits ...
   -A  --apphot          : run aperture photometry on a source list
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
   -S  --subbgd          : subtract background from image
   -v  --verbose         : display verbose outputs

   --> Single run commands
       --init                     : Initialise Starbug (post install)
       --generate-psf             : Generate ALL the PSF files to "PSFDIR"
       --local-param              : Make a local copy of the default parameter file
       --generate-region   a.fits : Make a ds9 region file with a detection file
       --clean-table       a.fits : Clean up an individual table
       --generate-run      *.fits : Generate a simple run script
       --version                  : Print starbug2 version

   --> typical runs
      $~ starbug2 -vD -p file.param image.fits      //Source detect on image with a parameter file
      $~ starbug2 -vDM -n4 images*.fits             //Source detect and match outputs of a list of images
      $~ starbug2 -vd image-ap.fits -BP image.fits  //PSF photometry on an image with a source file (image-ap.fits)
```

See starbug-manual.pdf for more detailed instructions.

