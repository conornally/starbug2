import starbug2
from starbug2.param import load_params, load_default_params
from starbug2.utils import *
from starbug2.misc import *
from starbug2.routines import *
from starbug2.artificialstars import Artificial_Stars

from astropy.wcs import WCS
from astropy.table import hstack, vstack
from photutils.psf import EPSFModel, subtract_psf, FittableImageModel


class StarbugBase(object):
    """
    StarbugBase is the overall container for the photometry package. It holds the active image,
    the parameter file and the output images/tables.
    It is self contained enough to simply run "photometry" and everything should just take care
    of itself from there on.
    """
    filter=None
    stage=0
    outdir=""
    fname=None
    bname=""

    detections=None
    psfcatalogue=None
    residuals=None
    background=None
    source_stats=None
    psf=None

    _image=None
    _nHDU=-1
    _unit=None
    wcs=None
    def __init__(self, fname, pfile=None, options={}):
        """
        fname : FITS image file name
        pfile : parameter file name
        options : extra options to load into starbug
        """
        if not pfile:
            if os.path.exists("starbug.param"): pfile="starbug.param"
            else: pfile=None
        self.options=load_params(pfile)
        self.options.update(options)
        self.load_image(fname)   ## Load the fits image

        if self.options["AP_FILE"]: self.load_apfile() ## Load the source list if given
        if self.options["BGD_FILE"]: self.load_bgdfile()

        #_=self.image ## Force self._nHDU

    @property
    def header(self):
        """
        Construct relevant base header information for routine products

        Returns
        -------
        `fits.Header` : Header file containing a series of relevvant information
        """
        head={}
        head["STARBUG"]=get_version()
        head["CALIBLEVEL"]=self.stage
        if self.filter: head["FILTER"]=self.filter
        head.update(self.options)
        head.update(self.info)
        return utils.collapse_header(head)


    @property
    def info(self):
        """
        Get some useful information from the image header file
        """
        out={}
        keys=("FILTER","DETECTOR","TELESCOP","INSTRUME",
              "BUNIT","PIXAR_A2", "PIXAR_SR")
        if self._image:
            for hdu in self._image:
                out.update( { (key,hdu.header[key]) for key in keys if key in hdu.header})

        return out

    @property
    def image(self):
        """
        automagically find the main image array to use
        Order of importance is:
        > self._nHDU (if set)
        > param[ HDUNAME ]
        > SCI, BGD, RES
        > first ImageHDU
        > image[0]
        """
        if self._nHDU >=0: return self._image[self._nHDU]
        enames=extnames(self._image)

        ## HDUNAME in param file
        n=self.options["HDUNAME"]
        if n and n in enames:
            self._nHDU=enames.index(n)
            return self._image[n]

        ##index?
        if type(n) in (int,float):
            self._nHDU=int(n)
            return self._image[self._nHDU]

        ## SCI, BGD, RES (common names)
        for n in ("SCI","BGD","RES"):
            if n in enames:
                self._nHDU=enames.index(n)
                return self._image[n]

        ## First ImageHDU
        for n,hdu in enumerate(self._image):
            if type[hdu]==fits.ImageHDU:
                self._nHDU=enames.index(n)
                return hdu

        self._nHDU=0
        return self._image[0]

    def log(self, msg):
        """
        Print message if in verbose mode 

        Parameters:
        -----------
        msg : str
            Message to print out
        """
        if self.options["VERBOSE"]:
            printf(msg)
            sys.stdout.flush()


    def load_image(self, fname):
        """
        Given fname, load the image into starbug to be worked on.

        Parameters:
        -----------
        fname : str
            Filename of fits image (with any number of extensions. If using
            a non standard HDU index, set the name or index of the extension
            with "HDUNAME=XXX" in the parameter file
        """
        self.fname=fname
        if fname:
            #########################################
            # Sorting out the file names and what not
            #########################################
            self.outdir,self.bname,extension=split_fname(fname)
            if (tmp_outname:=self.options.get("OUTPUT")) and tmp_outname !='.':
                outdir,bname,_=split_fname(tmp_outname)
                if os.path.exists(outdir) and os.path.isdir(outdir): self.outdir=outdir
                else: perror("unable to locate output directory \"%s\"\n"%outdir)
                if bname: self.bname=bname

            if extension==".fits":
                if os.path.exists(fname):
                    self.log("loaded: \"%s\"\n"%fname)
                    self._image=fits.open(fname)
                    _=self.image ## Force assigning _nHDU
                    self.log("-> using image HDU: %d (%s)\n"%(self._nHDU,self.image.name))
                    if self.image.data is None:
                        warn("Image seems to be empty.\n")

                    if (val:=self.header.get("TELESCOP")) is None or (val.find("JWST")<0):
                        warn("Telescope not JWST, there may be undefined behaviour.\n")

                    self.filter=self.options.get("FILTER")
                    if ("FILTER" in self.header) and (self.header["FILTER"] in starbug2.filters.keys()):
                        self.filter=self.header["FILTER"]
                        if self.options["FWHM"]<0: self.options["FWHM"]=starbug2.filters[self.filter].pFWHM
                    if self.filter:
                        self.log("-> photometric band: %s\n"%self.filter)
                    else:
                        warn("Unable to determine image filter\n")

                    if "DETECTOR" in self.info.keys():
                        self.log("-> detector module: %s\n"%self.info["DETECTOR"])
                    else: warn("Unable to determine Telescope DETECTOR.\n")

                    if "BUNIT" in self.image.header:
                        self._unit=self.image.header["BUNIT"]
                    else: warn("Unable to determine image BUNIT.\n")

                    self.wcs=WCS(self.image.header)

                    ## I NEED TO DETERMINE BETTER WHAT STAGE IT IS IN
                    exts=extnames(self._image)
                    if "DQ" in exts:
                        if "AREA" in exts: self.stage=2
                        else: self.stage=2.5
                    elif "WHT" in exts: self.stage=3
                    elif "CALIBLEVEL" in self.image.header: self.stage=self.image.header["CALIBLEVEL"]
                    else:
                        warn("Unable to determine calibration level, assuming stage 3\n")
                        self.stage=3


                    #self.log("loaded: \"%s\"\n"%fname)
                    self.log("-> pipeline stage: %d\n"%self.stage)

                else: warn("fits file \"%s\" does not exist\n"%fname)
            else: warn("included file must be FITS format\n")

    def load_apfile(self,fname=None):
        """
        Load a AP_FILE to be used during photometry

        Parameters:
        -----------
        fname : str
            Filename for fits table containg source coordinates. These coorindates can be
            xcentroid/ycentroid, x_init/y_init, x_0,y_0 or RA/DEC. The latter is used if 
            starbug gets "USE_WCS=1" in the parameter file.
        """
        if not fname: fname=self.options["AP_FILE"]
        if os.path.exists(fname):
            self.detections=import_table(fname)
            cn=set(self.detections.colnames)

            self.log("loaded AP_FILE='%s'\n"%fname)

            if self.options.get("USE_WCS"):
                if len(cn & set(("RA","DEC")))==2:
                    self.log("-> using RADEC coordinates\n")
                    try:
                        xy=self.wcs.all_world2pix(self.detections["RA"], self.detections["DEC"],0)
                    except:
                        warn("Something went wrong converting WCS to pixels, trying wcs_world2pix next.\n")
                        xy=self.wcs.wcs_world2pix(self.detections["RA"], self.detections["DEC"],0)
                    if "xcentroid" in cn: self.detections.remove_column("xcentroid")
                    if "ycentroid" in cn: self.detections.remove_column("ycentroid")
                    self.detections.add_columns(xy,names=("xcentroid","ycentroid"),indexes=[0,0])
                else:
                    warn("No 'RA' or 'DEC' found in AP_FILE\n")
                    #self.options["USE_WCS"]=0

            elif len( set(("x_0","y_0"))&cn)==2:
                self.detections.rename_columns(("x_0","y_0"),("xcentroid","ycentroid"))
            elif len( set(("x_init","y_init"))&cn)==2:
                self.detections.rename_columns(("x_init","y_init"),("xcentroid","ycentroid"))

            if len( set(("xcentroid","ycentroid"))&cn)==2:
                mask=(self.detections["xcentroid"]>=0) & (self.detections["xcentroid"]<self.image.shape[0]) & (self.detections["ycentroid"]>=0) & (self.detections["ycentroid"]<self.image.shape[1])
                self.detections.remove_rows(~mask)
                self.log("-> loaded %d sources from AP_FILE\n"%len(self.detections))
            else:
                warn("Unable to determine physical coordinates from detections table\n")
        else: perror("AP_FILE='%s' does not exists\n"%fname)

    def load_bgdfile(self,fname=None):
        """
        Load a BGD_FILE to be used during photometry
        
        Parameters:
        -----------
        fname : str
            Filename of fits image the same dimensions as the main image
        """
        if not fname: fname=self.options["BGD_FILE"]
        if os.path.exists(fname):
            self.background=fits.open(fname)[1]
            self.log("loaded BGD_FILE='%s'\n"%fname)
        else: perror("BGD_FILE='%s' does not exist\n"%fname)

    def load_psf(self,fname=None):
        """
        Load a PSF_FILE to be used during photometry

        Parameters:
        -----------
        fname : str
            Filename of a PSF fits image
        """
        status=0
        if not fname:
            fltr=starbug2.filters.get(self.filter)
            if fltr:
                dtname=self.info["DETECTOR"]
                if dtname=="NRCALONG": dtname="NRCA5"
                if dtname=="NRCBLONG": dtname="NRCB5"
                if dtname=="MULTIPLE":
                    if   fltr.instr==starbug2.NIRCAM and fltr.length==starbug2.SHORT: dtname="NRCA1"
                    elif fltr.instr==starbug2.NIRCAM and fltr.length==starbug2.LONG:  dtname="NRCA5"
                    elif fltr.instr==starbug2.MIRI:  dtname=""
                if dtname=="MIRIMAGE": dtname=""
                fname="%s/%s%s.fits"%(starbug2.DATDIR,self.filter,dtname)
            else: status=1
        if os.path.exists(fname):
            fp=fits.open(fname)

            if fp[0].data is None: 
                perror("There is a version mismatch between starbug and webbpsf. Please reinitialise with: starbug2 --init.\n")
                quit("Fatal error, quitting\n")


            self.psf=fp[0].data ####hmm
            fp.close()
            self.log("loaded PSF_FILE='%s'\n"%(fname))
        else:
            perror("PSF_FILE='%s' does not exist\n"%fname)
            status=1
        return status

    def prepare_image_arrays(self):
        """
        Make a copy of the original image, and prepare the other image arrays

        Returns:
        --------
        image : 
            A copy of the image array, scale into Jy if appropriate

        error : np.array
            An error array in the same unit as image

        bgd : np.array or None
            A copy of the loaded background image in the same units as the image
            
        mask : ?
        """

        # Collect scale factor
        if self.header.get("BUNIT")=="MJy/sr":
            scalefactor=get_MJysr2Jy_scalefactor(self.image)
            self.log("-> converting unit from MJy/sr to Jr with factor: %e\n"%scalefactor)
        else: scalefactor=1

        image=self.image.data.copy() * scalefactor

        # scale by area
        if "AREA" in extnames(self._image):
            image*= self._image["AREA"].data ## AREA distortion correction

        # collect and scale error
        if "ERR" in extnames(self._image) and np.shape(self._image["ERR"]):
            error=self._image["ERR"].data.copy() * scalefactor
        else: error=np.sqrt(np.abs(image))
        
        # create mask
        if "DQ" in extnames(self._image):
            mask=self._image["DQ"].data & (DQ_DO_NOT_USE|DQ_SATURATED) #|DQ_JUMP_DET)
            mask=mask.astype(bool)
        else:mask=np.isnan(image)

        # collect and scale background array
        if self.background is not None:
            bgd=self.background.data.copy() * scalefactor
        else: bgd=None

        return image, error, bgd, mask

    def detect(self):
        """
        Full source detection routine
        Saves the result as a table self.detections
        """
        self.log("Detecting Sources\n")
        status=0
        if self.image:# and self.filter:
            _f=starbug2.filters.get(self.filter)
            if self.options["FWHM"]>0: FWHM=self.options["FWHM"]
            elif _f: FWHM=_f.pFWHM
            else: FWHM=2
            #FWHM=_f.pFWHM if _f else self.options["FWHM"]
            #FWHM=starbug2.filters.get(self.filter).pFWHM

            detector=Detection_Routine( sig_src=self.options["SIGSRC"],
                                        sig_sky=self.options["SIGSKY"],
                                        fwhm=FWHM,
                                        sharplo=self.options["SHARP_LO"],
                                        sharphi=self.options["SHARP_HI"],
                                        round1hi=self.options["ROUND1_HI"],
                                        round2hi=self.options["ROUND2_HI"],
                                        smoothlo=self.options["SMOOTH_LO"], 
                                        smoothhi=self.options["SMOOTH_HI"],
                                        ricker_r=self.options["RICKER_R"],
                                        dobgd2d=self.options["DOBGD2D"],
                                        doconvl=self.options["DOCONVL"],
                                        boxsize=int(self.options["BOX_SIZE"]),
                                        cleansrc=self.options["CLEANSRC"],
                                        verbose=self.options["VERBOSE"])

            self.detections=detector(self.image.data.copy())["xcentroid","ycentroid","sharpness","roundness1","roundness2"]

            ra,dec=self.wcs.all_pix2world(self.detections["xcentroid"], self.detections["ycentroid"],0)
            self.detections.add_column( Column(ra, name="RA"), index=2)
            self.detections.add_column( Column(dec, name="DEC"), index=3)
            self.detections.meta=dict(self.header.items())
            self.detections.meta.update({"ROUNTINE":"DETECT"})
            self.aperture_photometry()

        else:
            perror("Something went wrong.\n")
            status=1
        return status



    def aperture_photometry(self):

        if self.detections is None:
            perror("No detection source file loaded (-d file-ap.fits)\n")
            return
        if len(set(("x_0","y_0","x_init","y_init","xcentroid","ycentroid")) & set(self.detections.colnames))<2:
            perror("No pixel coordinates in source file\n")
            return

        new_columns=("smoothness","flux","eflux","sky", "flag", self.filter,"e%s"%self.filter)
        self.detections.remove_columns( set(new_columns)&set(self.detections.colnames) )


        #######################
        # APERTURE PHOTOMETRY #
        #######################
        self.log("\nRunning Aperture Photometry\n")

        image, error, _, mask = self.prepare_image_arrays()

        """
        image=self.image.data.copy() ##dont work on the real image!

        #########################
        # Unit Conversion to Jy #
        #########################
        error=None
        if self.header.get("BUNIT")=="MJy/sr":
            scalefactor=get_MJysr2Jy_scalefactor(self.image)
            self.log("-> converting unit from MJy/sr to Jr with factor: %e\n"%scalefactor)
        else: scalefactor=1

        image*=scalefactor
        if "ERR" in extnames(self._image) and np.shape(self._image["ERR"]):
            error=self._image["ERR"].data.copy()
            error*=scalefactor
        else: 
            error=np.sqrt(image)
        """

        #######################
        # Aperture Correction #
        #######################
        apcorr=1
        apcorr_fname=None
        if (_apcorr_fname:=self.options.get("APCORR_FILE")): apcorr_fname=_apcorr_fname
        elif   self.info.get("INSTRUME")=="NIRCAM": apcorr_fname="%s/apcorr_nircam.fits"%starbug2.DATDIR
        elif self.info.get("INSTRUME")=="MIRI":   apcorr_fname="%s/apcorr_miri.fits"%starbug2.DATDIR

        if apcorr_fname: self.log("-> apcorr file: %s\n"%apcorr_fname)
        else: 
            warn("No apcorr file available for instrument\n")

        radius=self.options["APPHOT_R"]
        eefrac=self.options["ENCENERGY"]
        skyin= self.options["SKY_RIN"]
        skyout=self.options["SKY_ROUT"]

        if eefrac >=0:
            radius=APPhot_Routine.radius_from_encenrgy(self.filter, eefrac, apcorr_fname)
            if radius >0: self.log("-> calculating aperture radius from encirlced energy\n")

        if radius <=0: 
            if (radius:=self.options["FWHM"])>0:
                self.log("-> using FWHM as aprture radius\n")
            else:
                radius=2

        apcorr=APPhot_Routine.calc_apcorr(self.filter, radius, table_fname=apcorr_fname, verbose=self.options["VERBOSE"])

        ##################
        # Run Photometry #
        ##################
        apphot=APPhot_Routine( radius, skyin, skyout, verbose=self.options["VERBOSE"])

        if "DQ" in extnames(self._image):
            dqflags=self._image["DQ"].data.copy()
        else: dqflags=None
        """
        if self.stage==2:
            if "AREA" in extnames(self._image):
                image*= self._image["AREA"].data ## AREA distortion correction
            if "DQ" in extnames(self._image):
                mask=self._image["DQ"].data & (DQ_DO_NOT_USE|DQ_SATURATED) #|DQ_JUMP_DET)
                image[mask]=np.nan
                error[mask]=np.nan
                dqflags=self._image["DQ"].data
            else: dqflags=None

            ap_cat=apphot(image, self.detections, error=error, dqflags=dqflags, apcorr=apcorr, sig_sky=self.options["SIGSKY"])

        else: ##stage 3 version
            ap_cat=apphot(image, self.detections, error=error, apcorr=apcorr, sig_sky=self.options["SIGSKY"])
        """
        ap_cat=apphot(image, self.detections, error=error, dqflags=dqflags, apcorr=apcorr, sig_sky=self.options["SIGSKY"])


        fltr=self.filter if self.filter else "mag"
        mag,magerr=flux2mag( ap_cat["flux"], ap_cat["eflux"])
        ap_cat.add_column(Column(mag+self.options.get("ZP_MAG"),fltr))
        ap_cat.add_column(Column(magerr,"e%s"%fltr))
        self.detections=hstack((self.detections,ap_cat))

        if self.options.get("CLEANSRC"):
            N=len(self.detections)
            if (smoothlo:=self.options.get("SMOOTH_LO")) != "": 
                self.detections.remove_rows( self.detections["smoothness"]<smoothlo)
            if (smoothhi:=self.options.get("SMOOTH_HI")) != "": 
                self.detections.remove_rows( self.detections["smoothness"]>smoothhi)
            if len(self.detections)!=N:
                self.log("-> removing %d sources outside SMOOTH range\n"%(N-len(self.detections)))

        reindex(self.detections)
        self.detections.meta["FILTER"]=self.filter

        if not self.options.get("QUIETMODE"):
            _fname="%s/%s-ap.fits"%(self.outdir, self.bname)
            self.log("--> %s\n"%_fname)
            export_table(self.detections, _fname, header=self.header)


    def bgd_estimate(self):
        """
        Estimate the background of the active image
        Saves the result as an ImageHDU self.background
        """
        self.log("\nEstimating Background\n")
        if self.detections:

            _f=starbug2.filters.get(self.filter)
            if self.options["FWHM"]>0: FWHM=self.options["FWHM"]
            elif _f: FWHM=_f.pFWHM
            else: FWHM=2

            bgd=BackGround_Estimate_Routine(self.detections,
                                            boxsize=int(self.options["BOX_SIZE"]),
                                            fwhm=FWHM,
                                            sigsky=self.options["SIGSKY"],
                                            bgd_r=self.options["BGD_R"],
                                            verbose=self.options["VERBOSE"])
            header=self.header
            header.update(self.wcs.to_header())
            self.background=fits.ImageHDU(data=bgd(self.image.data.copy()), header=header)
            _fname="%s/%s-bgd.fits"%(self.outdir, self.bname)
            self.log("--> %s\n"%_fname)
            self.background.writeto(_fname,overwrite=True)
        else:
            perror("unable to estimate background, no source list loaded\n")



    def bgd_subtraction(self):
        """
        Internally subtract a background array from an image array
        """
        self.log("Subtracting Background\n")

        if self.background is None:
            perror("No background array loaded (-b file-bgd.fits)\n")
            return
        array= self.image.data - self.background.data
        self.residuals = array
        self._image[self._nHDU].data=array
        header=self.header
        header.update(self.wcs.to_header())
        fits.ImageHDU(data=self.residuals, name="RES", header=header).writeto("%s/%s-res.fits"%(self.outdir,self.bname), overwrite=True)

    def photometry(self):
        """
        Full photometry routine
        Saves the result as a table self.psfcatalogue
        // Additionally it appends a residual Image onto the self.residuals HDUList
        """
        if self.image:
            self.log("\nRunning PSF Photometry\n")


            image, error, bgd, mask = self.prepare_image_arrays()

            if bgd is None:
                _,median,_=sigma_clipped_stats(image,sigma=self.options["SIGSKY"])
                bgd=np.ones(self.image.shape)*median
                self.log("-> no background file loaded, measuring sigma clipped median\n")






            ###################################
            # Collect relevent files and data #
            ###################################
            #image=self.image.data.copy()

            if self.detections is None:
                perror("unable to run photometry: no source list loaded\n")
                return

            if self.psf is None and self.load_psf(os.path.expandvars(self.options["PSF_FILE"])):
                perror("unable to run photometry: no PSF loaded\n")
                return

            """
            if self.background is None:
                _,median,_=sigma_clipped_stats(self.image.data,sigma=self.options["SIGSKY"])
                bgd=np.ones(self.image.shape)*median
                self.log("-> no background file loaded, measuring sigma clipped median\n")
            else:
                bgd = self.background.data.copy()

            if "ERR" in extnames(self._image) and np.shape(self._image["ERR"]):
                error=self._image["ERR"].data
            else: 
                error=None
            """


            #_scalefactor=self.image.header.get("PHOTMJSR")#https://spacetelescope.github.io/jdat_notebooks/notebooks/psf_photometry/NIRCam_PSF_Photometry_Example.html
            #_bunit=self.image.header.get("BUNIT")
                #self.log("-> PHOTMJSR: %f\n"%_scalefactor)
            """
            if self.header.get("BUNIT")=="MJy/sr":
                scalefactor=get_MJysr2Jy_scalefactor(self.image)
                self.log("-> converting unit from MJy/sr to Jr with factor: %e\n"%scalefactor)
                image*=scalefactor
                bgd*=scalefactor
                if error is not None: error*=scalefactor
            else: scalefactor=1
            """

            psfmask= ~np.isfinite(self.psf)
            if psfmask.sum():
                self.psf[psfmask]=0
                self.log("-> masking INF pixels in PSF_FILE\n")

            psf_model=FittableImageModel(self.psf)
            if self.options["PSF_SIZE"]>0: size=int(self.options["PSF_SIZE"])
            else: size=psf_model.shape[0]
            if not size%2: size-=1
            self.log("-> psf size: %d\n"%size)

            #########################
            # Sort out Init guesses #
            #########################

            init_guesses=self.detections.copy()
            if "xcentroid" in init_guesses.colnames: init_guesses.rename_column("xcentroid", "x_init")
            if "ycentroid" in init_guesses.colnames: init_guesses.rename_column("ycentroid", "y_init")

            init_guesses=init_guesses[ init_guesses["x_init"]>=0 ]
            init_guesses=init_guesses[ init_guesses["y_init"]>=0 ]
            init_guesses=init_guesses[ init_guesses["x_init"]<self.image.header["NAXIS1"]]
            init_guesses=init_guesses[ init_guesses["y_init"]<self.image.header["NAXIS2"]]

            ######
            # Allow tables that dont have the correct columns through
            ######
            required=["x_init","y_init","flux",self.filter, "flag"]
            for notfound in  set(required)-set(init_guesses.colnames):
                dtype=np.uint16 if notfound=="flag" else float
                init_guesses.add_column( Column( np.zeros(len(init_guesses)), name=notfound, dtype=dtype) )

            init_guesses=init_guesses[required]
            init_guesses.remove_column("flux")
            init_guesses.rename_column(self.filter,"ap_%s"%self.filter)
            #init_guesses.rename_column("flux","flux_0")
            #init_guesses=init_guesses[init_guesses["flux_0"]>0]
            #init_guesses.remove_column("flux_0")

            
            ###########
            # Run Fit #
            ###########
            apphot_r=self.options.get("APPHOT_R")
            if not apphot_r or apphot_r <=0: apphot_r=3

            min_separation=self.options.get("CRIT_SEP")
            if min_separation is None:
                min_separation = min(5, 2.5*self.options.get("FWHM"))

            if self.options["FORCE_POS"]:
                phot=PSFPhot_Routine(psf_model, size, min_separation=min_separation, apphot_r=apphot_r, background=bgd, force_fit=1, verbose=self.options["VERBOSE"])
                psf_cat=phot(image,init_params=init_guesses, error=error, mask=mask)
                psf_cat["flag"] |= starbug2.SRC_FIX

            else:
                phot=PSFPhot_Routine(psf_model, size, min_separation=min_separation, apphot_r=apphot_r, background=bgd, force_fit=0, verbose=self.options["VERBOSE"])
                psf_cat=phot(image,init_params=init_guesses, error=error, mask=mask)



                ##################################
                # Setting position max variation #
                ##################################
                maxydev,unit=utils.parse_unit(self.options["MAX_XYDEV"])
                if unit is not None:
                    if unit==starbug2.DEG: 
                        maxydev*=60
                        unit=starbug2.ARCMIN
                    if unit==starbug2.ARCMIN:
                        maxydev*=60
                        unit=starbug2.ARCSEC
                    if unit==starbug2.ARCSEC:
                        if not self.header.get("PIXAR_A2"):
                            warn("MAX_XYDEV is units arcseconds, but starbug cannot locate a pixel scale in the header. Please use syntax MAX_XYDEV=%sp to set change to pixels\n"%maxydev)
                        else: maxydev /= np.sqrt(self.header.get("PIXAR_A2"))

                if maxydev>0:
                    self.log("-> position fit threshold: %.2gpix\n"%maxydev)
                    phot=PSFPhot_Routine(psf_model, size, min_separation=min_separation, apphot_r=apphot_r, background=bgd, force_fit=1, verbose=self.options["VERBOSE"])
                    ii=np.where( psf_cat["xydev"]>maxydev)
                    fixed_centres= psf_cat[ii][["x_init","y_init","ap_%s"%self.filter,"flag"]]
                    if len(fixed_centres):
                        self.log("-> forcing positions for deviant sources\n")
                        fixed_cat=phot(image, init_params=fixed_centres, error=error, mask=mask)
                        fixed_cat["flag"]|=starbug2.SRC_FIX
                        psf_cat.remove_rows(ii)
                        psf_cat=vstack((psf_cat, fixed_cat))
                    else: self.log("-> no deviant sources\n")


            ra,dec=self.wcs.all_pix2world(psf_cat["x_fit"], psf_cat["y_fit"],0)
            psf_cat.add_column( Column(ra, name="RA"), index=2)
            psf_cat.add_column( Column(dec, name="DEC"), index=3)

            ##################
            # Residual Image #
            ##################

            if self.options["GEN_RESIDUAL"]:
                self.log("-> generating residual\n")
                _tmp=psf_cat["x_fit","y_fit","flux"].copy()
                _tmp.rename_column("flux","flux_fit")
                residual = subtract_psf(image-bgd, psf_model, _tmp, subshape=(size,size))
                self.residuals=residual/get_MJysr2Jy_scalefactor(self.image)
                header=self.header
                header.update(self.wcs.to_header())
                fits.ImageHDU(data=self.residuals, name="RES", header=header).writeto("%s/%s-res.fits"%(self.outdir,self.bname), overwrite=True)

            #psf_cat.rename_column("flux_fit","flux")
            mag,magerr=flux2mag(psf_cat["flux"],psf_cat["eflux"])

            fltr= self.filter if self.filter else "mag"
            psf_cat.add_column(mag+self.options.get("ZP_MAG"),name=fltr)
            psf_cat.add_column(magerr,name="e%s"%fltr)
            self.psfcatalogue=tabppend(self.psfcatalogue, psf_cat)
            self.psfcatalogue.meta=dict(self.header.items())
            self.psfcatalogue.meta["AP_FILE"]=self.options["AP_FILE"]
            self.psfcatalogue.meta["BGD_FILE"]=self.options["BGD_FILE"]

            reindex(self.psfcatalogue)
            if not self.options.get("QUIETMODE"):
                _fname="%s/%s-psf.fits"%(self.outdir, self.bname)
                self.log("--> %s\n"%_fname)
                fits.BinTableHDU(data=self.psfcatalogue, header=self.header).writeto(_fname,overwrite=True)

    def artificial_stars(self):
        """
        Run artificial star testing

        >>> This needs to get the background loaded into it somewhere!!
        >>> Need to make sure all the appropriate parameters are being passed into the functions
        """
        status=0
        self.log("\nArtificial Star Testing (n=%d)\n"%(self.options["NTESTS"]))
        
        ################################
        # Collect files and sort units #
        ################################
        image=self.image.data.copy()
        bgd=None

        if self.background is not None:
            bgd=self.background.data.copy()

        if self.header.get("BUNIT")=="MJy/sr-1":
            scalefactor=get_MJysr2Jy_scalefactor(self.image)
            image/=scalefactor
            if bgd is not None: bgd/=scalefactor

        self.load_psf(self.options.get("PSF_FILE"))
        psf_model=FittableImageModel(self.psf)

        #############################
        # Build the Routine Classes #
        #############################

        detector=Detection_Routine( sig_src=self.options["SIGSRC"],
                                    sig_sky=self.options["SIGSKY"],
                                    fwhm=starbug2.filters[self.filter].pFWHM,
                                    sharplo=self.options["SHARP_LO"],
                                    sharphi=self.options["SHARP_HI"],
                                    round1hi=self.options["ROUND1_HI"],
                                    verbose=0)
        phot=APPhot_Routine ( self.options["APPHOT_R"],
                              self.options["SKY_RIN"],
                              self.options["SKY_ROUT"])
        
        phot=PSFPhot_Routine( psf_model, psf_model.shape,
                apphot_r=self.options["APPHOT_R"], force_fit=False, 
                background=bgd, verbose=0)
        export_table(phot(image,detector(image)),fname="/tmp/out.fits")



        art=Artificial_Stars(detector=detector, photometry=phot, psf=psf_model)

        ###########
        # Execute #
        ###########
        """
        I can parallelise here
        """

        MINMAG=27
        MAXMAG=18

        ZP = self.options.get("ZP_MAG") if self.options.get("ZP_MAG") else 0
        min_flux=np.exp( (np.log(10)/2.5)*(ZP-MINMAG) )
        max_flux=np.exp( (np.log(10)/2.5)*(ZP-MAXMAG) )

        #print("calc", min_flux, max_flux, MINMAG, MAXMAG)
        #print("cat", np.nanmin( self.detections["flux"]), np.nanmax( self.detections["flux"]), np.nanmax(self.detections["F444W"]), np.nanmin(self.detections["F444W"]))


        result= art.run_auto( image, background=bgd, ntests=self.options.get("NTESTS"), stars_per_test=self.options.get("NSTARS"), 
                            subimage_size=self.options.get("SUBIMAGE"), flux_range=( min_flux, max_flux))

        _fname="%s/%s-afs.fits"%(self.outdir, self.bname)
        export_table(result, fname=_fname)

        return status





    def source_geometry(self):
        """
        Calculate source geometry stats for a given image and source list
        """
        if self.detections is None: perror("No source file loaded\n")
        else:
            self.log("Running Source Geometry\n")
            slist=self.detections[["xcentroid","ycentroid"]].copy()
            slist=slist[ slist["xcentroid"]>=0 ]
            slist=slist[ slist["ycentroid"]>=0 ]
            slist=slist[ slist["xcentroid"]<self.image.header["NAXIS1"]]
            slist=slist[ slist["ycentroid"]<self.image.header["NAXIS2"]]

            sp=SourceProperties(self.image.data, slist, verbose=self.options["VERBOSE"])
            stat=sp(fwhm=starbug2.filters[self.filter].pFWHM, do_crowd=self.options["CALC_CROWD"])
            #geom=sp.calculate_geometry(fwhm=starbug2.filters[self.filter].pFWHM)
            
            self.source_stats=hstack((slist,stat))
            _fname="%s/%s-stat.fits"%(self.outdir, self.bname)
            self.log("--> %s\n"%_fname)
            reindex(self.source_stats)
            fits.BinTableHDU(data=self.source_stats,header=self.header).writeto(_fname, overwrite=True)
    
    def calculate_psf(self):
        """
        """
        pass



    def verify(self):
        """
        This simple function verifies that everything necessary has been loaded properly
        RETURN: 
            0 - on success
            1 - on fail
        """
        status=0
        #warn=lambda :perror(sbold("WARNING: "))
        
        self.log("Checking internal systems..\n")

        if not self.filter:
            warn("No FILTER set, please set in parameter file or use \"-s FILTER=XXX\"\n")
            status=1

        dname = os.path.expandvars(starbug2.DATDIR)
        if not os.path.exists(dname):
            warn("Unable to locate STARBUG_DATDIR='%s'\n"%dname)
            #status=1

        if not os.path.exists(self.outdir):
            warn("Unable to locate OUTPUT='%s'\n"%self.outdir)
            status=1

        tmp=load_default_params()
        if set(tmp.keys()) - set(self.options.keys()):
            warn("Parameter file version mismatch. Run starbug2 --update-param to update\n")
            status=1
        
        if self._image is None or self.image.data is None:
            warn("Image did not load correctly\n")
            status=1

        if self.options["AP_FILE"] and self.detections is not None:
            test=self.detections[["xcentroid","ycentroid"]]
            test=test[ test["xcentroid"]>=0 ]
            test=test[ test["ycentroid"]>=0 ]
            test=test[ test["xcentroid"]<self.image.header["NAXIS1"]]
            test=test[ test["ycentroid"]<self.image.header["NAXIS2"]]
            if not len(test):
                warn("Detection file empty or no sources overlap the image.\n")
                status=1

        return status

    def __getstate__(self):
        self._image.close()
        #if self.background: self.background.close()
        state=self.__dict__.copy()
        if "_image" in state:
            del state["_image"] ##Sorry but we cant have that
        if "background" in state: ## This currently doesnt get reloaded
            del state["background"]
            
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        v=self.options["VERBOSE"]
        self.options["VERBOSE"]=0
        self.load_image(self.fname)
        self.options["VERBOSE"]=v
