# StarBugII

JWST PSF photometry in dusty crowded fields.
Last updated: v0.4.2

[![Python application](https://github.com/conornally/starbug2/actions/workflows/python-app.yml/badge.svg)](https://github.com/conornally/starbug2/actions/workflows/python-app.yml)
[![PyPI version fury.io](https://badge.fury.io/py/starbug2.svg)](https://pypi.python.org/pypi/starbug2/)
[![Latest release](https://badgen.net/github/release/conornally/starbug2)](https://github.com/conornally/starbug2/releases)


## Installation

```bash
$~ pip install starbug2

--- OR ---

$~ git clone https://github.com/conornally/starbug2.git
$~ cd starbug2
$~ python -m build
$~ pip install .
```

After the package is installed, there are a few steps required to initialise Starbug.

**WEBBPSF** Is a dependency of Starbug that has its own initialisation process. The full installation is documented on https://webbpsf.readthedocs.io/en/latest/installation.html however it requires two main steps. Download the data file on the website, named something like webbpsf-data-X.X.X.tar.gz and expand it into a directory, then add append to your .bashrc (or equivalent) `export "WEBBPSF_PATH=PATH/TO/DIRECTORY"`.

**DATA FILES** Starbug needs to generate the WEBBPSFs, and collect some CRDS, to do this run `starbug2 --init`. It will generate these files by default into "${HOME}/.local/share/starbug" however if you wish to use a different directory, set the environment variable "STARBUG_DATDIR" to the desired destination.

```bash
$~ echo "export 'WEBBPSF_PATH=PATH/TO/WEBBPSF/DIRECTORY'" >> ~/.bashrc
$~ echo "export 'STARBUG_DATDIR=PATH/TO/DESTINATION'" >> ~/.bashrc

$~ starbug2 --init
```

Finally verify the installation by running `starbug2 --version`

## Usage

```bash
Starbug II - JWST PSF photometry
usage: starbug2 [-ABDfGhMPSv] [-b bgdfile] [-d apfile] [-n ncores] [-o directory] [-p file.param] [-s opt=val] image.fits ...
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
       --generate-psf             : Generate ALL the PSF files to "STARBUG_DATDIR"
       --generate-region   a.fits : Make a ds9 region file with a detection file
       --generate-run      *.fits : Generate a simple run script
       --clean-table       a.fits : Clean up an individual table
       --version                  : Print starbug2 version

       --apply-zeropint    a.fits : Apply a zeropoint (-s ZEROPOINT=1.0) to a.fits
       --calc-instr-zp     a.fits : Calculate and apply an instrumental zero point onto a.fits

   --> typical runs
      $~ starbug2 -vD -p file.param image.fits      //Source detect on image with a parameter file
      $~ starbug2 -vDM -n4 images*.fits             //Source detect and match outputs of a list of images
      $~ starbug2 -vd image-ap.fits -BP image.fits  //PSF photometry on an image with a source file (image-ap.fits)
```
```bash
StarbugII Matching 
usage: starbug2-match [-BGfhv] [-o output] [-p file.param] [-s KEY=VAL] table.fits ...
    -B  --band               : match in "BAND" mode (does not preserve a column for every frame)
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
```

See [starbug-manual](https://github.com/conornally/starbug2/blob/main/docs/starbug-manual.md) for more detailed instructions.
