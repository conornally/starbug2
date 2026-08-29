# StarBugII

<img src="docs/source/_static/images/starbug.png" align="left" width="250px"/>


## PSF photometry in complex and crowded fields

&emsp; • Detect deeply dust-embedded sources in NIRCam and MIRI images

&emsp; • Full suite of photometry tools

&emsp; • Powerful modular and simple GNU Linux standard command line interface


<br clear="left">

[![Python application](https://github.com/conornally/starbug2/actions/workflows/python-app.yml/badge.svg)](https://github.com/conornally/starbug2/actions/workflows/python-app.yml)
[![PyPI version fury.io](https://badge.fury.io/py/starbug2.svg)](https://pypi.python.org/pypi/starbug2/)
[![Latest release](https://badgen.net/github/release/conornally/starbug2)](https://github.com/conornally/starbug2/releases)
[![Documentation Status](https://readthedocs.org/projects/starbug2/badge/?version=latest)](https://starbug2.readthedocs.io/en/latest/?badge=latest)
[![ASCL.net](https://img.shields.io/badge/ascl-2309.012-blue.svg?colorB=262255)](https://ascl.net/2309.012)


## Installation

*StarbugII* requires [webbpsf](https://webbpsf.readthedocs.io/en/latest/installation.html) installed on the system to function fully.

```bash
$~ pip install starbug2
```

NOTE::
At this point star bug is fully functional for telescopes which are not JWST.
To bring into line for use with JWST, the "JWST initialization" section
needs to be followed, which can be done using the following command line.
```
$~ starbug2 --init
```

NOTE:: I found the stpsf stuff to fail due to outdated dependnecy. needed to 
do the following:
```
 pip install --upgrade pysiaf

```

### From source

If installing starbugII from source, the following commands are used to ensure
current dependencies are required. The venv name is not essential to the install
but was set up for consistency.

- ```python3 -m venv star_bug_env```
- ```./star_bug_env/bin/pip install --upgrade stpsf==2.2.0```
- ```./star_bug_env/bin/pip install numpy==2.3.5 photutils==2.0.1 --force```
- ```./star_bug_env/bin/pip install parse==1.22.0```
- ```./star_bug_env/bin/pip install setuptools==82.0.1```
- ```./star_bug_env/bin/pip install scikit-image==0.26.0```
- ```./star_bug_env/bin/pip install webbpsf==2.0.0```
- ```./star_bug_env/bin/pip install pytest==9.0.3```
- ```./star_bug_env/bin/python -m pip install build```

instead of:
- ```pip install .```
I used the following command to ensure live development.
- ```./star_bug_env/bin/python -m pip install -e . --no-deps```

For generating the doc images. graphviz is needed. To isntall that , run 
the following commands:
``` pip isntall graphviz```

> [!IMPORTANT]
> If you make use of *StarbugII* in any published or presented work,
> please include a [citation](https://ui.adsabs.harvard.edu/abs/2023ascl.soft09012N/abstract).
> 
> *StarbugII* uses methods and datatypes from [astropy](https://docs.astropy.org/en/stable/) and 
> [photutils](https://photutils.readthedocs.io/en/stable), please acknowledge them accordingly.

## Documentation

See the [full documentation](https://starbug2.readthedocs.io/en/latest/?badge=latest) for the complete installation and detailed 
guides to using the photometric routines. Basic usage information is 
produced by running:

```txt
StarbugII - JWST PSF photometry
usage: starbug2 [-ABDfGhMPSv] [-b bgdfile] [-d apfile] [-n ncores] [-o ouput] [-p file.param] [-s opt=val] image.fits ...
   -A  --apphot          : run aperture photometry on a source list
   -B  --background      : run background estimation
   -D  --detect          : run source detection
   -G  --geom            : calculate geometric stats on source list
   -M  --match           : match outputs from all input image files
   -P  --psf             : run psf photometry
   -S  --subbgd          : subtract background from image

   -b  --bgdfile         : load background (-bgd.fits) file
   -d  --apfile  ap.fits : load a source detection (-ap.fits) file to skip the source detection step
   -f  --find            : attempt to find associated -ap -bgd files
   -h  --help            : display uasage information
   -n  --ncores      num : number of CPU cores to split process between
   -o  --output      dir : output directory
   -p  --param   a.param : load parameter file
   -s  --set      option : set value in parameter file at runtime (-s SIGSKY=3)
   -v  --verbose         : display verbose outputs

   --> Single run commands
       --init                     : Initialise Starbug (post install)
       --local-param              : Make a local copy of the default parameter file
       --update-param             : Update an out-of-date local parameter file
       --generate-psf             : Generate a single PSF. Set FILTER, DET_NAME, PSF_SIZE with -s
       --generate-region   a.fits : Make a ds9 region file with a detection file
       --generate-run      *.fits : Generate a simple run script
       --version                  : Print starbug2 version

   --> typical runs
      $~ starbug2 -vD -p file.param image.fits      //Source detect on image with a parameter file
      $~ starbug2 -vDM -n4 images*.fits             //Source detect and match outputs of a list of images
      $~ starbug2 -vd image-ap.fits -BP image.fits  //PSF photometry on an image with a source file (image-ap.fits)

To see more detailed information on an option, run [OPTION] --help:
    $~ starbug2 -D --help

See https://starbug2.readthedocs.io for full documentation.
```

```txt
StarbugII Matching 
usage: starbug2-match [-BCGfhvX] [-e column] [-m mask] [-o output] [-p file.param] [-s KEY=VAL] table.fits ...
    -B  --band               : match in "BAND" mode (does not preserve a column for every frame)
    -C  --cascade            : match in "CASCADE" mode (left justify columns)
    -G  --generic            : match in "GENERIC" mode
    -X  --exact              : match in "EXACTVALUE" mode

    -e  --error   column     : photometric error column ("eflux" or "stdflux")
    -f  --full               : export full catalogue
    -h  --help               : show help message
    -m  --mask    eval       : column evaluation to mask out of matching e.g. -m"~np.isnan(F444W)"
    -o  --output  file.fits  : output matched catalogue
    -p  --param   file.param : load starbug parameter file
    -s  --set     option     : set value in parameter file at runtime (-s MATCH_THRESH=1)
    -v  --verbose            : display verbose outputs

        --band-depr          : match in "old" band mode

    --> typical runs
       $~ starbug2-match -Gfo outfile.fits tab1.fits tab2.fits
       $~ starbug2-match -sMATCH_THRESH=0.2 -sBRIDGE_COL=F444W -Bo out.fits F*W.fits
```
