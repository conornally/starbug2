# StarBugII

JWST PSF photometry in dusty crowded fields.

## Installation

```bash
pip install starbug2-X.X.X.tar.gz

-- OR --

pip install --user starbug2-X.X.X.tar.gz 
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
See starbug-manual.pdf for more detailed instructions.

