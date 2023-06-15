# StarBugII Manual

JWST PSF photometry in dusty crowded fields.
Last updated: v0.4.3

## Installation

```bash
$~ pip install starbug2
```

After the package is installed, there are a few steps required to initialise Starbug.

**WEBBPSF** Is a dependency of Starbug that has its own initialisation process. The full installation is documented on [webbpsf homepage](https://webbpsf.readthedocs.io/en/latest/installation.html) however it requires two main steps. Download the data file on the website, named something like webbpsf-data-X.X.X.tar.gz and expand it into a directory, then append to your .bashrc (or equivalent) `export "WEBBPSF_PATH=PATH/TO/DIRECTORY"`.

**DATA FILES** Starbug needs to generate the WEBBPSFs, and collect some CRDS, to do this run `starbug2 --init`. It will generate these files by default into "${HOME}/.local/share/starbug" however if you wish to use a different directory, set the environment variable "STARBUG_DATDIR" to the desired destination.

```bash
$~ echo "export 'WEBBPSF_PATH=PATH/TO/WEBBPSF/DIRECTORY'" >> ~/.bashrc
$~ echo "export 'STARBUG_DATDIR=PATH/TO/DESTINATION'" >> ~/.bashrc
$~ starbug2 --init

//verify it starts up without issue
$~ starbug2 --version
```

### TAB Completion

StarbugII has a bash completion script `starbug2/extras/starbug2.completion`. This can be installed directly into `/etc/bash_completion.d/` or `"source /PATH/TO/COMPLETION/FILE"` can be place within your .bashrc. Unfortunately this completion script works only in bash shells.


## Usage

```bash
Starbug II - JWST PSF photometry
usage: starbug2 [-ABCDfhMPv] [-b bgdfile] [-d apfile] [-n ncores] [-o directory] [-p file.param] [-s opt=val] image.fits ...
   -A  --apphot          : run aperture photometry on a source list
   -B  --background      : run background estimation
   -b  --bgdfile         : load background (-bgd.fits) file
   -C  --clean           : run source cleaning before photometry 
   -d  --apfile  ap.fits : load a source detection (-ap.fits) file to skip the source detection step
   -D  --detect          : run source detection
   -f  --find            : attempt to find associated -ap -bgd files
   -G  --geom            : calculate geometric stats on source list
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
       --generate-psf             : Generate ALL the PSF files to "STARBUG_DATDIR"
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

## Parameter File

The parameter file is where any dataset specific parameters can be tweaked. Ideally the default values should be sufficient however if the diffuse dust emissions cause a very complex background or the field is very crowded, certain parameters may need tuned. Additionally differing sensitivities between photometric bands may require different parameters to detect on.
To generate a local parameter file, run:

```bash
$~ starbug2 --local-param
```

This will create `./starbug.param`, a file which will be loaded by default when starbug is ran from the folder that contains it. However it can be named anything and explicitly loaded at runtime with `-p file.param` in the command. It may be the case that you keep several parameter files on hand to conduct slightly different routines, this is how you would load them.

For quick testing or tweaking of certain parameters, starbug can override a file setting with the addition of `-s PARAM=VALUE` OR `--set PARAM=VALUE` in the command. This will not change the parameter file but will use this parameter value instead.

### Settable Parameters

As of the current version. If your parameter file doesnt fit the template of the current version of the default file, starbug will warn you but may crash later if you dont update the local parameter file.

| NAME             | DTYPE        | DESCRIPTION                                                           |
|------------------|--------------|---------------------------------------------------------------------------------|
| VERBOSE          | INT 0:1      | Include verbose outputs.                                              |
| OUTPUT           | STR          | Output file or folder to send all products to.                        |
|------------------|--------------|-----------------------------------------------------------------------|
| SIGSKY           | FLOAT >0     | Number of sigma below which pixels gets removed as sky. In images with bright diffuse emissions, drop this value slowly, in steps of ~0.1 and watch the number of sources detected increase, but be careful of false detections of bright spots. |
| SIGSRC           | FLOAT >0     | Minimum number of sigma above the median that a source must be to be detected. This is often 5sigma for a robust search or 3sigma for faint detections. |
| DOBGD2D          | INT 0:1      | Run complex background subtraction during source detection (This usually results in more sources detected, however it takes a long time). |
| SHARP_LO         | FLOAT >0     | Set bounds of source sharpness, outside which sources will be ignored.|
| SHARP_HI         | FLOAT >0     | Set bounds of source sharpness, outside which sources will be ignored.|
| ROUND_LO         | FLOAT >0     | Set bounds of source roundness, outside which sources will be ignored.|
| ROUND_HI         | FLOAT >0     | Set bounds of source roundness, outside which sources will be ignored.|
| CLEANSRC         | INT 0:1      | in prep.|
| RICKER_R         | FLAOT >0     | in prep.|
|------------------|--------------|-----------------------------------------------------------------------|
| FIT_APP_R        | INT 0:1      | Fit fraction encircled energy to aperture radius and use this as aperture radii (1), OR, use explicit aperture radius (0). |
| ENCENERGY        | FLOAT >0     | Fraction encircled energy to fit aperture radius to.                  |
| APPHOT_R         | FLOAT >0     | Aperture photometry radius.                                           |
| SKY_RIN          | FLOAT >0     | Sky annulus inner radius.                                             |
| SKY_ROUT         | FLOAT >0     | Sky annulus outer radius.                                             |
|------------------|--------------|-----------------------------------------------------------------------|
| BGD_R            | - | DESCRIPTION |
| BOX_SIZE         | INT >1 [pix] | Kernel size during background subtraction to estimate background within. For complex fields this is set low. |
|------------------|--------------|-----------------------------------------------------------------------|
| AP_FILE          | STR          | Load source list file (-ap.fits) into starbug. This is equivalent to `-d file-ap.fits`.|
| BGD_FILE         | STR          | Load background estimation file (-bgd.fits) into starbug. This is equivalent to `-b file-bgd.fits`.|
| PSF_FILE         | STR          | DESCRIPTION |
| CRIT_SEP         | INT >0       | DESCRIPTION |
| FORCE_POS        | INT 0:1      | Conduct forced centroid photometry, if set no (0) then starbug will also fit centroid positions |
| DPOS_THRESH      | FLOAT >0 [arcsec] | If PSF photometry fits centroid positions that deviate from the original positions by a threshold greater than this value (in units arcsec), these sources will have PSFs refit with forced centroids. |
| PSF_SIZE         | INT >0       | Set a custom PSF to fit, by default it will take the dimensions of the WEBBPSF file. |
| GEN_RESIDUAL     | INT 0:1      | Generate a residual images with the fitted PSFs removed. |
|------------------|--------------|-----------------------------------------------------------------------|
| MATCH_THRESH     | FLOAT >0 [arcsec] | Separation threshold between coordinate during astrometric matching. Set low to avoid mismatching.|
| MATCH_COLS       | STR,STR,..   | Comma separated list of columns to include in outputs during matching.  |
| RM_MATCH         | INT          | (prep.) |
| NEXP_THRESH      | INT          | (prep.) |
| SN_THRESH        | FLOAT >0     | (prep.) |
| BRIDGE_COL       | FILTER       | (prep.) |
|------------------|--------------|-----------------------------------------------------------------------|
| NUMBER_ARTIFICIAL_STARS | - | DESCRIPTION |
| SUBIMAGE_SIZE    | - | DESCRIPTION |
| MIN_FLUX         | - | DESCRIPTION |
| MAX_FLUX         | - | DESCRIPTION |
| SEPARATION_THRESH| - | DESCRIPTION |
|------------------|--------------|-----------------------------------------------------------------------|
| REGION_COL       | STR          | DS9 region colour.                                                    |
| REGION_SCAL      | INT 0:1      | Scale region radius with flux?                                        |
| REGION_RAD       | FLOAT >0     | If REGION_SCAL=0, what is the region radius.                          |
| REGION_XCOL      | STR          | Table column name for X coordinate of region.                         |
| REGION_YCOL      | STR          | Table column name for Y coordinate of region.                         |
| REGION_WCS       | INT 0:1      | If the REGION_(X/Y)COL values are WCS. If not it will set them as pixel coordinates.|

\newpage

## A Typical Run

A typical run begins with a MAST download, I recommend using [jwst_mast_query](https://github.com/spacetelescope/jwst_mast_query), and it is worth splitting the data by photometric filter. 
For each filter, it is worth creating a unique parameter file with:
```bash
$~ starbug2 --local-param
```
This will place `starbug.param` in the current folder. If this file exists in the folder at runtime, starbug will automatically load it. However it can be renamed to anything as long as it is explicitly loaded with `-p file.param` during any subsequent commands.

### Source detection

Begin honing the source detection parameters by running starbug on a single exposure. This will be an iterative process of detecting and tweaking parameters. If the background of the image is fairly flat or uniform, the parameter `DOBGD2D` can be set to 0 to turn off the slower background subtraction pass during the routine. 

```bash
$~ starbug2 -vD [-p param.file] exposure1.fits

//To create a DS9 region file to look at the resulting source file
$~ starbug2 --generate-region exposure1-ap.fits
$~ ds9 exposure1.fits -regions exposure1-ap.fits
```
If the resulting source list seems to miss obvious sources, begin slowly lowering `SIG_SRC` and `SIG_SKY`. The former is usually left at 5 for a robust detections, below this we begin to get into the uncertainty limit of the data, however in certain cases this may be desirable, although be careful dropping it below 3. 
It situations where the background contains more complex dust emissions, or is has some non-uniformity, `SIG_SKY` can be a productive parameter to lower. Drop it slowly, by no more than 0.1 at a time, and watch the output source file carefully. This will drastically increase the detection of fainter sources however it will begin to make false detections of bright dusty peaks.

To remove (a significant number of) spurious sources and or resolved background galaxies, the geometric parameters `SHARP_HI/LO` and `ROUND_HI/LO` are available. Sharpness is a measure of how the peak of the source compares to the median within the FWHM. Cosmic rays are often very sharp and can be limited by lowering the upper limit. Whereas other artifacts often appear less sharp. It is worth opening the output source list `*-ap.fits` and plotting the distribution of `SHARPNESS`, usually there is a clear roughly normal distribution with wings - set `SHARP_HI` and `SHARP_LO` to cut off these wings. 
Roundness is a measure of eccentricity of a source. The distribution should be symmetric (`ROUND_HI` and `ROUND_LO` and measures of the same thing in orthogonal directions). Sources centred on 0 are round or pointlike, whereas resolved galaxies are often towards the edges of the distribution.

### Introducing Dithers

You should have now arrived on appropriate parameters but may still be finding that certain extraneous bad sources are squeaking through the cuts. The next method we can deploy to remove these is rolling in the set of dithers.
We can detect simultaneously on a set of dithers and match the outputs together, removing any sources that don't appear a threshold number of times. 

If you give a starbug instance a list of files to work on, it will by default run them in series and produce a separate source list for each. However we can parallelise and match with the `-n NUMBER_CORES` and `-M` options:

```bash

$~ starbug2 -vD -M exposure1.fits exposure2.fits exposure3.fits exposure4.fits

// parallelised 
$~ starbug2 -vD -M -n 4 exposure1.fits exposure2.fits exposure3.fits exposure4.fits

//General bash wildcarding etc. can also be used
$~ starbug2 -vDMn4 exposure{1..4}.fits
$~ starbug2 -vDMn4 exposure*.fits
```
The matched output will export into two catalogues: `exposure(1234)-apfull.fits` and `exposure(1234)-apmatch.fits`, the former being the complete catalogue containing every column, the latter being a condensed form that cuts sources without a threshold "NUM" value. If parameter **RM_MATCH** is set then the `-apmatch.fits` catalogue will remove any sources that haven't been detected in at least that many exposures.

### Loading a Source List

Once a source has been created, it can be loaded back into starbug. Source lists don't have to be in a specific shape or form or have been created by starbug however they must include positional columns. These can be pixel coordinates **x_0 y_0** or **xcentroid ycentroid** which will be used before any world coordinates **RA DEC**.
There are several methods to loading a source list into starbug to work on. It will sometimes be referred to as an AP file.
Most commonly it will occur in the command line using the argument `-d` or `--apfile` which will take require the name of the fits table after.

```bash
$~ starbug2 -d sourcelist.fits ....
```

If the source list is being used for many starbug runs, it can be added to the parameter file **AP_FILE**=/PATH/TO/FILE. This will take lower priority than loading it from the command line, so it can be overloaded at runtime.

Finally, starbug can automagically find a source list associated with an image file with the `-f` or `--find`  command line option. This will look for a fits table with "-ap.fits" added to the image filename. For example "image.fits" will try and locate "image-ap.fits". 

```bash
$~ starbug2 -vf image.fits
>>>	loaded AP_FILE='./image-ap.fits'
```

Source lists are required by several starbug routines; Aperture Photometry, Background Estimation, PSF Photometry. Either generated before, or all the routines can be rolled in and conducted in the same starbug run.

### Aperture Photometry

Aperture photometry is conducted as part of the source detection step, however it can be ran in isolation is required with `-A` or `--apphot`:

```bash
$~ starbug2 -vA -d sources.fits image.fits

// or on several files
$~ starbug2 -vA -d sources.fits image1.fits image2.fits ...
```

There are two modes that aperture photometry can be ran in. Setting an aperture radius or scaling the radius with percentage encircle energy. The former allows for constant radii between all the photometric bands within a dataset but introduces errors on the fit of the aperture correction which will deteriorate at very small, or large radii. Instead we can scale the aperture radius with encircled energy, this will change the aperture radius for every photometric band but will have a more solid aperture correction solution. To switch between modes, toggle **FIT_APP_R** in the parameter file, (0 ignores aperture radius and uses encircled energy, 1 ignores encircled energy and uses aperture radius).

If using the fixed aperture radius method, then **APPHOT_R** sets this value in pixel units. Alternatively set the fraction encircled energy with **ENCENERGY** using a number between zero and 1.
Regardless of the photometric method, set the sky annulus radii with **SKY_IN SKY_OUT** in pixel units.

### Background Estimation

To run the background estimation routine, first load or generate a source list and use the `-B` or `--background` argument. The resulting background file "(-bgd)" can be reloaded into starbug similarly to a source list with `-b` or `--bgdfile` or `-f` in the command line or **BGD_FILE** in the parameter file.

```bash
$~ starbug2 -d sources.fits -B image.fits

//run together to create a sourcelist and bgd file
$~ starbug2 -vBD image.fits

//load it into starbug later
$~ starbug2 -b image-bgd.fits ...
```


### PSF Photometry

PSF photometry can be ran in two modes, free or fixed centroids, which allows or otherwise the source position to be refit as well as the flux. To change between the two modes, toggle the parameter `FORCE_POS`= 0 (free positions) 1 (forced photometry). 
Running the photometry requires a source list of initial positions and a diffuse background estimation, both discussed above. A simple run of the photometry will look like:

```bash
$~ starbug2 -d sourcelist.fits -b background.fits -P image.fits
 OR
$~ starbug2 -vDBP image.fits
```

Starbug will automatically use the simulated WEBBPSF file associated with the image detector and module, however custom PSF files can be loaded instead by supplying the filename to a fits file with `PSF_FILE` in the parameter file. The routine will use the entire PSF array to fit each source in the sourcelist, if the PSF is very large this will slow the program significantly. Reduce the size of the fitting by setting `PSF_SIZE`.

There are two phases to free centroid PSF photometry: first all the sources are fit, then the routine inspects the new calculated source positions. If the resulting centroid has moved from the initial position by a distance greater than `DPOS_THRESH` (arcsec), a second round of fixed centroid photometry occurs using the initial positions, these sources will be flagged `SRC_FIX` in the resulting catalogue.
The quality of the PSF fits can be checked by inspecting the residual images. To generate them set `GEN_RESIDUAL=1` and the image `image-res.fits` will be produced with the background estimation and PSFs removed from the original image. 

Finally, a long standing crash in the routine occurs when many sources are being fit at the same time, causing a recursion error. This happens more often in compact fields in short wavelength images. It can be mitigated by reducing the size of `CRIT_SEP`.


### Instrumental Zero Points

PSF photometry is much more accurate but is not output into physical units. The instrumental zero point must be calculated from the aperture photometry. To do this we run:

```bash
$~ starbug2 -d APFILE.fits --calc-instr-zp table.fits
```
Where `APFILE.fits` is a very clean catalogue in the units you wish to move to. This can just be the same source list that has been used throughout however it may be prudent to clean the poorer signal to noise sources from this, as every source is weighted equally.
The resulting catalogue will not overwrite the input one but create a new `*-zp.fits` file.






\newpage

## A Full Run Through (in command line)

```bash
$~ ls
>	F115W/ F200W/ F444W

$~ cd F444W; ls
>	image1_1.fits image1_2.fits image1_3.fits image1_4.fits
>	image2_1.fits image2_2.fits image2_3.fits image2_4.fits

//determine good parameters
$~ starbug2 --local-param
$~ starbug2 -vD image1_1.fits
$~ starbug2 --generate-run image*_{1..4}.fits

>>> set CMDS to desired routines (something like '-DMn4')
$~ source run.sh
$~ starbug2-match -oF444W-ap.fits *-apmatch.fits

//now with a global sourcelist, we can do some background estimation and PSF photometry
$~ starbug2 -vd F444W-apmatch.fits -B image1_{1..4}.fits 
$~ starbug2 -vd F444W-apmatch.fits -fPMn4 image1_{1..4}.fits

>>> OR using run.sh set CMDS again (something like '-vd F444W-ap.fits -BPMn4')
$~ starbug2-match -oF444W-psf.fits *-psfmatch.fits
```

Run a similar set of commands for each photometric band, then match the outputs

```bash
$~ starbug2-match -vB -o out-psf.fits F115W/F115W-psf.fits F200W/F200W-psf.fits F444W/F444W-psf.fits
```




\newpage


## Source Flags

Sources are given quality flags at various points in starbug routines. These flags are applied as bitwise masks, so that the source retains information from all stages during data reduction. 

| BIT  | NAME      | DESCRIPTION                                                                                        |
|------|-----------|----------------------------------------------------------------------------------------------------|
|`0x00`| SRC_GOOD  | Source OK                                                                                          |
|`0x01`| SRC_BAD   | Source aperture contains a pixel marked as saturate or bad                                         |
|`0x02`| SRC_JMP   | Source aperture contains a pixel marked with a jump during integration (possible cosmic ray)       |
|`0x04`| SRC_VAR   | Source has an asymmetric flux distribution between matches (mean and median more than 5% different)|
|`0x08`| SRC_FIX   | Source has PSF photometry with a forced position                                                   |
|`0x10`| SRC_UKN   | Something was wrong..                                                                              |

## FAQ


