v0.5.3:
    Background:
        >   BGD_R added functionality
        >   Changed automatic radius estimation
    Param:
        >   Added ZP_MAG: which is added to magnitude columns ## This needs some thinking
        >   Added a differentiating roundness 1 and 2
        >   Removing ROUND_LO
        >   Adding SMOOTH bounds
    APPHOT:
        >   Catching another error in the aperture masking
        >   load_apfile now cuts extra sources outside fov
        >   load_apfile gives errors in the correct place now
    MATCHING:
        >   Fixing no default matching routine
    PSF:
        >   Updated PFS base to photutils 1.9 classes
    General:
        >   Starting to imprse minimum package versio requirments onto some dependencies 
    BIN:
        >   Moved the main starbug script to a console entry point
v0.5.2:
    APPHOT:
        >   error calculation takes same form as stsci example (+sumsq)
        >   sky measurement now done arraywise, speeds up execution ten fold
v0.5.1:
	Matching:
		>	Takes stdflux as the error instead of eflux (at least for now)


v0.5.0:
    Generalisation
        Detection and APhot
            > APCORR_FILE added, so that you can use custom values
            > Many catches for non existant values
            > Removing FIT_AP_R
            > Removing assumption of AREA and DQ flags
        PSF
            > started on generalising that too
            > Gnenerated PSFs are now PrimaryHDUs
            > DPOS_THRESH now supports pixel or arcsec units
            > MAX_XYDEV taking the place of DPOS_THRESH but with default units in pixels
            > utils.parse_unit handy for this
            > BGD files dont need to be generated, if left off it will take a sigma clipped median
            > fixed the bad removal of good stars during the forced photometry
        MATCHING:
            > Fixes to generic_match so that it as able to nan fill tables better
        PARAMS:
            > Added APCORR_FILE, FILTER and FWHM

    MATCH:
        > Default match changed from cascade to generic
    
    
    

            



v0.4.9:
    Bug fixes:
        >   webbpsf update has caused problems with F150W2 generation
            --> Ive just added a try catch on it for now
	Photometry:
		>	Updated CRDS MIRI apcorr download
		>	Switched APCORR calculation to interpolation rather than curve fitting
v0.4.8:
	Bug fixes:
		>	generic match has a filled(np.nan) bodge applied
v0.4.7:
	Bug fixes:
		>	bandmatch and calc-instr-zp had had bugs introduced with small changes to generic_match earlier

v0.4.2:
    Matching:
        >   Adding SN_THRESH to band matching
        >   Formally adding BRIDGE_COL to band matching
        >   Band match now exports -nircam -miri catalogues mid way through the matching
    General:
        >   Added "import_table()" but not implemented throughout

v0.4.1:
	Adding NIRCam to MIRI matching
v0.4.0:
	Matching:
		>	Starting work on the miri matching properly
		>	right now there is a hard coded increasing matching threshold wrt fwhm
			>>> nircam short 0.06
			>>>	nircam long  0.1
			>>> <=F560W		 0.15
			>>>	<=F1000W	 0.2
			>>>	<=F2100W	 0.5
v0.3.15:
	Match
		>	Genericmatch column names now drop the suffix if possible
v0.3.14:
	Param:
		>	Added NEXP_THRESH
        >   Beginning to remove RM_MATCH
    MATCH:
        >   Allowing NEXP_THRESH to be used during the normal matching
        >   Bandmatch just does magnitudes now
    General:
        >   --apply-zeropoint
        >   including a motd at the start now

v0.3.13:
	General:
		>	change argument -P --photom --> -P --psf
        >   *-apmatch etc files should now also go to the right outdir
        >   --calc-instr-zp added
        
v0.3.12:
    Param:
        >   Moved Boxsize to bgd estimation
        >   Changed OUTDIR to OUTPUT and can now change basename of output files too
    Detection and PSF:
        >   Fixed output column order
        >   Moved export table to ends of routines rather than in a separate export function
        
v0.3.11:
    General:
        >   Printing stdout for first file of -n set, hoping this will help diagnose the silent crash when memory runs out
	PSF:
		>	fitting now happens on the array in Jy not DN
v0.3.10:
	PSF:
		>	non conforming sourcelists now allowed
	General:
		>	If a apfile doesnt have sources or there isnt any overlap with the image, then
			starbug.verify will fail that run

v0.3.8:
	PSF:
		>	Bug(fix), recursion error during fitting source groups with more members than the
			system recursion limit. For now the recursion limit has just been overrided.
v0.3.7:
    Detection:
        >   changing convolution detection method to find peaks
		>	Added RICKER_R to params
    Utils
        >   I think update-param should be more robust..
v0.3.6:
	PSF:
		>	Trying to correct the sporadic PSF fit fail that occurs when -fP is given
v0.3.4:
    MISC:
        >   adding --update-param
v0.3.2/3:
	Mainly bug fixes 
    Source Detection:
        >   Added the convolution method
        >   Removing sourceExtractor method

v0.3.1:
	Detection:
		problem: some sources are being discarded because they have nan rnd/shp values. I dont know why they have calculated badly
		>	Used SourceProperties to fit shapes at the end of the routine
        >   Tring to remove the need for WCS in the detection
            >>> Need to finalise match threshold
	SourceProperties:
		>	Plugging in source properties routine
			>>> Crowding
			>>>	DAOPHOT stats
v0.3.0:
    PSF:
        >   load_psf: can now load a non default PSF_FILE
		>	PSF loading should now account or "MULTIPLE" DETECTOR key in mosaics
		>	added PSF_SIZE so can reduce the size fo the psf fitting, this speeds everything up a lot
        >   changed bgd BUNIT to DN
		>	Removed init flux guess for now
    General:
        >   Trying to edit MANIFEST.in to include files properly
		>	Loading in a sourcelist will ignore flag column if it doesnt exist
		>	Removing PSFDIR from param file and exchanging it for an optional environment variable
		>	Little bit of drafted work in the docs
v0.2.19:
    PSF:
        >   I fixed the fitting!
    General:
        >   fixed some verbose outputs
v0.2.18:
	bugs:
		>	SourceProperties verbose flag
v0.2.17:
	Fixing some readme information
v0.2.16:
    PSF:
        >   generate residual between forced and free xy fits
		>	new PSF param options 
			FORCE_POS
			DPOS_THRESH
			GEN_RESIDUAL
	GENERAL:
		>	can set HDUNAME in param file
		>	residuals are now just a single hdu 
v0.2.15:
    PSF:
        >   split forced and free xy position

v0.2.14:
	General: 
		>	I broke export_table earlier, I have now fixed it
		>	Fixing scipy mode future warning
	DOCS:
		>	Beginning to rework the manual
v0.2.13:
	General:
		Removed the := operator
v0.2.12:

v0.2.11:
    APPHOT:
        >   Can now ask to use closest encircled energy instead of fitting a raidus
    General:
        >   Added some configuration to region generation
v0.2.10:
    Starbug:
        >   Adding caliblevel into headers
        >   INIT now downloads (crudely) the CRDS apcorr files for nircam and miri
	APPHOT:
		>	BIG Simplifying error calculation
        >   Using error array instead of sqrt data
        >   picks nircam or miri apcor file
	MATCHING:
		>   changing name "Generic match" to cascade match
        >   first pass stage matching
    PSF!
        >   Started work on PSF routine again
        >   For a single exposure, it can load in the ap and bgd file and run RELATIVELY well
v0.2.9:
    General:
        >   Fixing a runtime warning on sqrt(nan) in apphot error array
        >   Turning off verbose on starbug.__setstate__
        >   match.generic_match adds magnitude properly
        >   removed peak from apphot catalogue
        >   region generation looks at flux not peak
v0.2.8:
    StarBug:
        >   Adding self.stage to track where the processes differ
        >   Can now handle stage two and stage three
    APPhot:
        >   Takes the arrays now, not the full hdulist
    Detection:
        >   Adding option to not do background 2d (needs streamlined)
    Match:
        >   Added Dither match

v0.2.7:
	Param:
		>	Changed the default values for APPHOT_R APPHOT_SKY* 
v0.2.6:
	General
		>Fixing the log10 warning in flux2mag
		>Adding a few tests
		>adding requirements.txt...yet more packaging files -.-
v0.2.5: I keep forgetting to update this file
		
v0.2.4:
    Minor tweaks

v0.2.3:
	Previously it successully installed with pip, but the param file failed, this hopefully fixes that
v0.2.2:
	Setting things up for PIP
v0.2.1:
    General
        >   Fix some MACOS compatibility errors
        >   Added OUTDIR to param file
v0.2.0:
	>	first pass aperture phot package..with some quirks
	>	--generate-run -> run.sh makes a simple run modifiable run script to aid the user
	>	matching stdev now (SHOULD) take the largest rather than combining

v0.1.7:
	>	APPHot local background estimated using MODE not median
	>	Can pass sigma sky into apphot to sigma clip background backgorund annulis
	>	band match first draft
	

v0.1.6:
	Aperture Photometry things:
		>	apcorr fits aperture radius
		>	errors
		>	source flags saturation or JUMP
		>	averaging now takes mean,std and numsq errors
        >   flag raised if the mean and median are more than 5% different for one source over frames (SRC_VAR)

v0.1.5:
	Matching has been put on hold as i need to add the ability to be able to run all the dithers on pipeline 2 images...

    matching.dither_match():
        >   this matches together all files input into bin/starbug on one run
        >   it outputs two files for now, -apfull (all matched columns) -apmatch (a stripped down version)

    bin/starbug multiprocessing:
        >   it will now multiprocess (not by default) all the starbug input file processes
        >   __setstate__ __getstate__ added to starbug to enable pickling

    matching.generic_match():
        >   match tables for different OBSRVN VISIT DETECTOR etc
    General:
        >   hcascade: left aligned hstack of tables
        >   find_colnames: look for substring in table colnames

v0.1.4:
	Matching:
        >   bin/starbug2-match initial draft
		>	added MATCH_COL to param file, to choose which columns to include in output
		>	added RM_MATCH to param file, to remove matches with less than RM_MATCH resolutions from table
		>	averaging now takes the median of the values in the table
		>	removed NULLVAL from the tables and just using nan for now
    AperturePhot:
        >   Done inside three radii, n=1..3 FWHM returning three different flux values
        >   Sky value now calculated as median in a much larger radius (default 25-30)
    General:
        >   move load_params to utils
        >   loading now has a print resolution
		>	made a first attempt at a completion script

v0.1.3:
	Background Estimation:
		>	The aperture size of the local sky estimation now scales with peak pixel value of the source
	Known issues:
		>	Need to make sure PSFPhot loads the background in properly
			>>>	sometimes it needs to load it, other times it needs to generate it
    AperturePhot:
        >   APPhot now takes HUDList not 2d array, in time a lot of routines should take this
        >   Flux multiplied by the area modifier in pipeline2
            >>> Need to do the same for PSFPhot
	General:
		>	region generation also scales aperture on peak source pixel

v0.1.2:
	--> The PSF photometry over subtracts the sources
		>	For every source in the given sourcelist, a mean sky value is calculated with an annuli (Rin=APHOT_R0, Rout=APHOT_R1 for now)
		>	any pixel less than 2FWHM from the centre of the source is then set to this value (i know this has some flaws)
		>	Background2D is then run over the image.
			>>> The idea is not to "remove" the star from the image, its to mitigate its effect during background2D calculations
        >   generate-region scales the circle to the source size
	Photometry:
		>	it now spits out the calculated background
		>	a new loading bar is visible during the source masking
    Background Estimation:
        >   The apertures to estimate the sky now scale with the peak pixel value of the source. This allows it to "hug" to the star, which will reduce the effect of crowding
		
v0.1.1:
	Photometry:
		>   saves residual images into a multiextension hdu
	Source Detection:
		>   fixed Background2D filtersize and boxsize inputs
		>   removed MedianBackground subtraction (no different from plain)
	Cleanup:
		>   now cutting when ap_flux<0    
		>	and then I removed that ^^^
	General:
		>   moved StarbugBase to its own file
		>   separated the filter info
		>   added some parameter file descriptions
		>	first draft on manual
		>   moved param file detection out of StarbuBase to bin file

v0.1.0 - a working version..

NOTES and TODO:

    >   During export step, add starbug parameters to the header file. That way the file will be completely reproducible
    >   Separate the background subtraction annuli radii from the apphot radii?
    >   Artificial stars needs to make sure it is using all the correct routines 
	>	Make sure you have the correct PSFs. I think the ones im currently using are too narrow? The residuals dont look fantastic
    >   Fix when verbose happens, or dont have it as an option at all, or make verbose mode a lot more verbose?

	>	Cleanup: Sharpness basically is not gaussian so the fit isnt very accurate

