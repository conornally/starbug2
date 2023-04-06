import starbug2
from starbug2.utils import *
from starbug2.misc import *
from starbug2.routines import *

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
    fname=None
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
        if not pfile: pfile="%s/default.param"%pkg_resources.resource_filename("starbug2","param/")
        self.options=load_params(pfile)
        self.options.update(options)

        self.load_image(fname)   ## Load the fits image
        if self.options["AP_FILE"]: self.load_apfile() ## Load the source list if given
        if self.options["BGD_FILE"]: self.load_bgdfile()

        #_=self.image ## Force self._nHDU

    @property
    def header(self):
        """
        """
        head=fits.Header({**self.options,**self.info})
        head["CALIBLEVEL"]=self.stage
        head["STARBUG"]=pkg_resources.get_distribution("starbug2").version
        head["FILTER"]=self.filter
        #head["BUNIT"]=self.image.header["BUNIT"]
        head.update(self.info)
        #head.update(self.wcs.to_header())
        return head


    @property
    def info(self):
        """
        Get some useful information from the image header file
        """
        out={}
        keys=("FILTER","DETECTOR","TELESCOP","INSTRUME",
              "BUNIT","PIXAR_A2")
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
        Print message if in verbose mode (just a macro really)
        INPUT:  msg=message to print out
        """
        if self.options["VERBOSE"]:
            printf(msg)
            sys.stdout.flush()


    def load_image(self, fname):
        """
        Given fname, load the image into starbug to be worked on.
        Will return None if file isnt readable
        """
        self.fname=fname
        if fname:
            dname,bname,extension=split_fname(fname)
            if extension==".fits":
                if os.path.exists(fname):
                    self.log("loaded: \"%s\"\n"%fname)
                    self._image=fits.open(fname)
                    _=self.image ## Force assigning _nHDU
                    self.log("-> using image HDU: %d (%s)\n"%(self._nHDU,self.image.name))

                    if ("FILTER" in self.header) and (self.header["FILTER"] in starbug2.filters.keys()):
                        self.filter=self.header["FILTER"]
                        self.log("-> photometric band: %s\n"%self.filter)
                    else: warn();perror("Unable to determine image filter\n")

                    if "DETECTOR" in self.info.keys():
                        self.log("-> detector module: %s\n"%self.info["DETECTOR"])
                    if "BUNIT" in self.image.header:
                        self._unit=self.image.header["BUNIT"]
                    self.wcs=WCS(self.image.header)

                    ## I NEED TO DETERMINE BETTER WHAT STAGE IT IS IN
                    exts=extnames(self._image)
                    if "DQ" in exts:
                        if "AREA" in exts: self.stage=2
                        else: self.stage=2.5
                    elif "WHT" in exts: self.stage=3
                    elif "CALIBLEVEL" in self.image.header: self.stage=self.image.header["CALIBLEVEL"]
                    else:
                        warn();
                        perror("Unable to determine jwst pipeline level, assuming 3\n")
                        self.stage=3


                    #self.log("loaded: \"%s\"\n"%fname)
                    self.log("-> pipeline stage: %d\n"%self.stage)

                else: warn();perror("fits file \"%s\" does not exist\n"%fname)
            else: warn();perror("included file must be FITS format\n")

    def load_apfile(self,fname=None):
        """
        Load a AP_FILE to be used during photometry
        INPUT:
            fname : file-ap.fits (this file is exported during source detection step
        """
        if not fname: fname=self.options["AP_FILE"]
        if os.path.exists(fname):
            self.detections=Table().read(fname,format="fits")#data=fp[1].data._get_raw_data())
            cn=self.detections.colnames

            if "flag" in cn:
                self.detections["flag"]=Column(self.detections["flag"], dtype=np.uint16)

            self.log("loaded AP_FILE='%s'\n"%fname)

            if not any( _ in cn for _ in ("xcentroid","ycentroid","x_0","y_0")):
                if all( _ in cn for _ in ("RA","DEC")):
                    xy=self.wcs.all_world2pix(self.detections["RA"], self.detections["DEC"],0)
                    self.detections.add_columns(xy,names=("xcentroid","ycentroid"),indexes=[0,0])
                    self.log("-> using RADEC coordinates\n")
                else: perror("WARNING, unable to determine physical coordinates from detections table\n")
            if len( set(("x_0","y_0"))&set(self.detections.colnames))==2: self.detections.rename_columns(("x_0","y_0"),("xcentroid","ycentroid"))
        else: perror("AP_FILE='%s' does not exists\n"%fname)

    def load_bgdfile(self,fname=None):
        """
        Load a BGD_FILE to be used during photometry
        INPUT:
            fname : file-bgd.fits (this file is exported during background estimation step
        """
        if not fname: fname=self.options["BGD_FILE"]
        if os.path.exists(fname):
            self.background=fits.open(fname)[1]
            self.log("loaded BGD_FILE='%s'\n"%fname)
        else: perror("BGD_FILE='%s' does not exist\n"%fname)

    def load_psf(self,fname=None):
        """
        Load a PSF_FILE to be used during photometry
        INPUT:
            fname : psf.fits
        """
        status=0
        if not fname:
            fltr=starbug2.filters[self.filter]
            dtname=self.info["DETECTOR"]
            #print(dtname)
            if dtname=="MULTIPLE":
                if   fltr.instr==starbug2.NIRCAM and fltr.length==starbug2.SHORT: dtname="NRCA1"
                elif fltr.instr==starbug2.NIRCAM and fltr.length==starbug2.LONG:  dtname="NRCALONG"
                elif fltr.instr==starbug2.MIRI:  dtname=""
            if dtname=="MIRIMAGE": dtname=""
            fname="%s/%s%s.fits"%(starbug2.DATDIR,self.filter,dtname)
        if os.path.exists(fname):
            fp=fits.open(fname)
            self.psf=fp[1].data ####hmm
            fp.close()
            self.log("loaded PSF_FILE='%s'\n"%(fname))
        else:
            perror("PSF_FILE='%s' does not exist\n"%fname)
            status=1
        return status


    def detect(self):
        """
        Full source detection routine
        Saves the result as a table self.detections
        """
        self.log("Detecting Sources\n")
        if self.image and self.filter:
            FWHM=starbug2.filters[self.filter].pFWHM
            detector=Detection_Routine( sig_src=self.options["SIGSRC"],
                                        sig_sky=self.options["SIGSKY"],
                                        fwhm=FWHM,
                                        sharplo=self.options["SHARP_LO"],
                                        sharphi=self.options["SHARP_HI"],
                                        roundlo=self.options["ROUND_LO"],
                                        roundhi=self.options["ROUND_HI"],
                                        bgd2d=self.options["DOBGD2D"],
                                        boxsize=int(self.options["BOX_SIZE"]),
                                        cleansrc=self.options["CLEANSRC"],
                                        verbose=self.options["VERBOSE"])

            self.detections=detector(self.image.data.copy())["xcentroid","ycentroid","sharpness","roundness1","roundness2"]

            ra,dec=self.wcs.all_pix2world(self.detections["xcentroid"], self.detections["ycentroid"],0)
            self.detections.add_column( Column(ra, name="RA"), index=1)
            self.detections.add_column( Column(dec, name="DEC"), index=2)
            self.aperture_photometry()

            self.detections.meta=dict(self.header.items())
            self.detections.meta.update({"ROUNTINE":"DETECT"})


    def aperture_photometry(self):

        if self.detections is None:
            perror("No detection source file loaded (-d file-ap.fits)\n")
            return
        if len(set(("x_0","y_0","xcentroid","ycentroid")) & set(self.detections.colnames))<2:
            perror("No pixel coordinates in source file\n")
            return

        new_columns=("flux","eflux","sky", "flag", self.filter,"e%s"%self.filter)
        self.detections.remove_columns( set(new_columns)&set(self.detections.colnames) )


        #######################
        # APERTURE PHOTOMETRY #
        #######################
        self.log("Running Aperture Photometry\n")
        image=self.image.data.copy() ##dont work on the real image!

        #########################
        # Unit Conversion to Jy #
        #########################
        error=None

        scalefactor=get_MJysr2Jy_scalefactor(self.image)
        self.log("-> converting unit from MJy/sr to Jr with factor: %e\n"%scalefactor)
        #if self._image["SCI"].header["BUNIT"]=="MJy/sr":
            #scalefactor=1e6*float(self._image["SCI"].header["PIXAR_SR"])

        image*=scalefactor
        if "ERR" in extnames(self._image) and np.shape(self._image["ERR"]):
            error=self._image["ERR"].data
            error*=scalefactor
        else: error=np.sqrt(image)

        #######################
        # Aperture Correction #
        #######################
        apcorr=1
        fname=None
        radius=self.options["APPHOT_R"]
        skyin= self.options["SKY_RIN"]
        skyout=self.options["SKY_ROUT"]

        if   self.info["INSTRUME"]=="NIRCAM": fname="%s/apcorr_nircam.fits"%starbug2.DATDIR
        elif self.info["INSTRUME"]=="MIRI":   fname="%s/apcorr_miri.fits"%starbug2.DATDIR
        else: perror("No apcorr file available for instrument\n")

        if self.options["FIT_APP_R"]:
            apcorr=APPhot_Routine.calc_apcorr(self.filter, self.options["APPHOT_R"], table_fname=fname, verbose=self.options["VERBOSE"])
        else:
            apcorr,radius=APPhot_Routine.apcorr_from_encenergy(self.filter,self.options["ENCENERGY"],table_fname=fname, verbose=self.options["VERBOSE"])

        apphot=APPhot_Routine( radius, skyin, skyout, encircled_energy=self.options["ENCENERGY"], fit_radius=self.options["FIT_APP_R"], verbose=self.options["VERBOSE"])

        if self.stage==2:
            image*= self._image["AREA"].data ## AREA distortion correction
            mask=self._image["DQ"].data & (DQ_DO_NOT_USE|DQ_SATURATED) #|DQ_JUMP_DET)
            image[mask]=np.nan
            error[mask]=np.nan
            ap_cat=apphot(self.detections, image, error=error, dqflags=self._image["DQ"].data, apcorr=apcorr, sig_sky=self.options["SIGSKY"])

        else: ##stage 3 version
            ap_cat=apphot(self.detections, image, error=error, apcorr=apcorr, sig_sky=self.options["SIGSKY"])

        mag,magerr=flux2ABmag( ap_cat["flux"], ap_cat["eflux"], filter=self.filter)
        ap_cat.add_column(Column(mag,self.filter))
        ap_cat.add_column(Column(magerr,"e%s"%self.filter))
        self.detections=hstack((self.detections,ap_cat))

    def bgd_estimate(self):
        """
        Estimate the background of the active image
        Saves the result as an ImageHDU self.background
        """
        self.log("Estimating Background\n")
        if self.detections:
            xname="xcentroid" if "xcentroid" in self.detections.colnames else "x_0"
            yname="ycentroid" if "ycentroid" in self.detections.colnames else "y_0"

            sources=self.detections[[xname,yname]]
            sources=sources[ sources[xname]>=0 ]
            sources=sources[ sources[yname]>=0 ]
            sources=sources[ sources[xname]<self.image.header["NAXIS1"]]
            sources=sources[ sources[yname]<self.image.header["NAXIS2"]]
            bgd=BackGround_Estimate_Routine(sources,
                                            boxsize=int(self.options["BOX_SIZE"]),
                                            fwhm=starbug2.filters[self.filter].pFWHM,
                                            verbose=self.options["VERBOSE"])
            self.background=fits.ImageHDU(data=bgd(self.image.data.copy()), header=self.wcs.to_header())
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

    def photometry(self):
        """
        Full photometry routine
        Saves the result as a table self.psfcatalogue
        // Additionally it appends a residual Image onto the self.residuals HDUList
        """
        if self.detections is None:
            perror("unable to run photometry: no source list loaded\n")
            return

        if self.background is None:
            perror("unable to run photometry: no background estimation loaded\n")
            return

        if self.psf is None and self.load_psf(self.options["PSF_FILE"]):
            perror("unable to run photometry: no PSF loaded\n")
            return

        if self.image:
            self.log("Running PSF Photometry\n")

            ###################################
            # Collect relevent files and data #
            ###################################

            image=self.image.data.copy()
            bgd = self.background.data.copy()

            _bunit=self.image.header.get("BUNIT")
            _scalefactor=self.image.header.get("PHOTMJSR")
            if _scalefactor:#https://spacetelescope.github.io/jdat_notebooks/notebooks/psf_photometry/NIRCam_PSF_Photometry_Example.html
                self.log("-> PHOTMJSR: %f\n"%_scalefactor)
                image/=_scalefactor
                bgd/=_scalefactor

            psf_model=FittableImageModel(self.psf)
            #psf_model=EPSFModel(fp[1].data)
            if self.options["PSF_SIZE"]>0: size=int(self.options["PSF_SIZE"])
            else: size=psf_model.shape[0]
            if not size%2: size-=1
            self.log("-> psf size: %d\n"%size)

            #########################
            # Sort out Init guesses #
            #########################

            init_guesses=self.detections.copy()
            if "xcentroid" in init_guesses.colnames: init_guesses.rename_column("xcentroid", "x_0")
            if "ycentroid" in init_guesses.colnames: init_guesses.rename_column("ycentroid", "y_0")

            init_guesses=init_guesses[ init_guesses["x_0"]>=0 ]
            init_guesses=init_guesses[ init_guesses["y_0"]>=0 ]
            init_guesses=init_guesses[ init_guesses["x_0"]<self.image.header["NAXIS1"]]
            init_guesses=init_guesses[ init_guesses["y_0"]<self.image.header["NAXIS2"]]
            init_guesses=init_guesses[["x_0","y_0","flux",self.filter, "flag"]]
            init_guesses.remove_column("flux")
            #init_guesses.rename_column("flux","flux_0")
            init_guesses.rename_column(self.filter,"ap_%s"%self.filter)
            #init_guesses=init_guesses[init_guesses["flux_0"]>0]
            #init_guesses.remove_column("flux_0")

            ###########
            # Run Fit #
            ###########

            _psf_cat=None
            _fixpsf_cat=None

            if not self.options["FORCE_POS"]:
                dpos= self.options["DPOS_THRESH"] / np.sqrt( self.image.header["PIXAR_A2"])
                self.log("-> position fit threshold [pix]: %.2g\n"%dpos)
                phot=PSFPhot_Routine(self.options["CRIT_SEP"], psf_model, size, background=bgd, force_fit=0, verbose=self.options["VERBOSE"])
                _psf_cat=phot(image,init_guesses=init_guesses)

                d = (_psf_cat["x_0"]-_psf_cat["x_fit"])**2.0 + (_psf_cat["y_0"]-_psf_cat["y_fit"])**2.0
                ii=np.where(d>=dpos**2.0)
                init_guesses=init_guesses[ii]
                _psf_cat.remove_rows(ii)
                if len(init_guesses): self.log("-> number bad position fits: %d\n"%len(init_guesses))

            if len(init_guesses):
                phot=PSFPhot_Routine(self.options["CRIT_SEP"], psf_model, size, background=bgd, force_fit=1, verbose=self.options["VERBOSE"])
                _fixpsf_cat=phot(image,init_guesses=init_guesses)
                _fixpsf_cat["flag"] |= starbug2.SRC_FIX

            if _psf_cat is not None and _fixpsf_cat is not None: psf_cat=vstack((_psf_cat,_fixpsf_cat))
            elif _psf_cat is None: psf_cat=_fixpsf_cat
            else: psf_cat=_psf_cat

            ra,dec=self.wcs.all_pix2world(psf_cat["x_fit"], psf_cat["y_fit"],0)
            psf_cat.add_column( Column(ra, name="RA"), index=1)
            psf_cat.add_column( Column(dec, name="DEC"), index=2)

            ##################
            # Residual Image #
            ##################

            if self.options["GEN_RESIDUAL"]:
                self.log("-> generating residual\n")
                residual = subtract_psf(image-bgd, psf_model, psf_cat[["x_fit","y_fit","flux_fit"]], subshape=(size,size))
                self.residuals=residual

            ######################
            # Photometric offset # takes the top 50% least crowded sources
            ######################

            #crowd=SourceProperties(self._image["SCI"].data, psf_cat[["RA","DEC"]]).calculate_crowding()
            #ii=np.argsort(crowd)[len(crowd)//2:]
            #apmag,_=flux2ABmag(psf_cat["apflux"],None,filter=self.filter)
            psf_cat.rename_column("flux_fit","flux")
            mag,magerr=flux2ABmag(psf_cat["flux"],psf_cat["eflux"],filter=self.filter)
            #dmag= np.nanmean( mag[ii]-apmag[ii] )
            #mag-=dmag
            #self.log("Photometric offset: %f\n"%dmag)

            psf_cat.add_column(mag,name=self.filter)
            psf_cat.add_column(magerr,name="e%s"%self.filter)
            self.psfcatalogue=tabppend(self.psfcatalogue, psf_cat)
            self.psfcatalogue.meta=dict(self.header.items())
            #self.background=fits.ImageHDU(data=phot.bkg_estimator.bgd, name="BACKGROUND", header=self.wcs.to_header()) ##So is it supposed to be a fits image or a numpy array?!

    def cleanup(self):
        """
        """
        self.log("Cleaning up..\n")

        if self.detections:
            cln=Cleaning_Routine(self.detections, verbose=self.options["VERBOSE"])
            self.detections=cln.run( mag_unc=self.options["ERROR_CUT"],
                                    sharp_sig_hi= self.options["SHARP_HI_SIG"],
                                    sharp_sig_lo= self.options["SHARP_LO_SIG"],
                                    round_sig_hi= self.options["ROUND_HI_SIG"],
                                    round_sig_lo= self.options["ROUND_LO_SIG"])
        if self.psfcatalogue:
            cln=Cleaning_Routine(self.psfcatalogue, verbose=self.options["VERBOSE"])
            self.psfcatalogue=cln.run( mag_unc=self.options["ERROR_CUT"],
                                    sharp_sig_hi= self.options["SHARP_HI_SIG"],
                                    sharp_sig_lo= self.options["SHARP_LO_SIG"],
                                    round_sig_hi= self.options["ROUND_HI_SIG"],
                                    round_sig_lo= self.options["ROUND_LO_SIG"])

    def artificial_stars(self):
        """
        Run artificial star testing

        >>> This needs to get the background loaded into it somewhere!!
        """

        fname=os.path.expandvars("%s/%s.fits"%(starbug2.DATDIR, self.filter))
        with fits.open(fname) as fp:
            psf_model=DiscretePRF(fp[0].data)

        detector=Detection_Routine( sig_src=self.options["SIGSRC"],
                                    sig_sky=self.options["SIGSKY"],
                                    fwhm=starbug2.filters[self.filter].pFWHM,
                                    sharplo=self.options["SHARP_LO"],
                                    sharphi=self.options["SHARP_HI"],
                                    roundlo=self.options["ROUND_LO"],
                                    roundhi=self.options["ROUND_HI"],
                                    wcs=WCS(self.image.header),
                                    verbose=0)

        phot=PSFPhot_Routine(   self.options["CRIT_SEP"],
                                starbug2.filters[self.filter].pFWHM,
                                psf_model,
                                psf_model.shape,
                                sig_sky=self.options["SIGSKY"],
                                sig_src=self.options["SIGSRC"],
                                sharplo=self.options["SHARP_LO"],
                                sharphi=self.options["SHARP_HI"],
                                roundlo=self.options["ROUND_LO"],
                                roundhi=self.options["ROUND_HI"],
                                wcs=WCS(self.image.header),
                                verbose=0)

        art=ArtificialStar_Routine(detector, phot, psf_model)
        self.log("Artificial Star Testing (n=%d)\n"%(self.options["NUMBER_ARTIFICIAL_STARS"]))
        result=art.run(self.image.data, ntests=self.options["NUMBER_ARTIFICIAL_STARS"], flux_range=(self.options["MIN_FLUX"], self.options["MAX_FLUX"]),
                subimage_size=self.options["SUBIMAGE_SIZE"], separation_thresh=self.options["SEPARATION_THRESH"], fwhm=starbug2.filters[self.filter].pFWHM)
        export_table(result, "/tmp/artificialstars.fits")

    def export(self, outdir=None):
        """
        Export all the current catalogues
        """
        if not outdir: outdir=self.options["OUTDIR"]
        if not os.path.exists("%s/"%outdir):
            perror("output directory '%s' does not exist, using /tmp instead\n"%outdir)
            outdir="/tmp"

        dname,fname,ext=split_fname(self.fname)
        if self.detections:
            self.detections.meta["FILTER"]=self.filter
            reindex(self.detections)
            hdulist=[fits.PrimaryHDU(header=self.header),fits.BinTableHDU(data=self.detections)]
            fits.HDUList(hdulist).writeto("%s/%s-ap.fits"%(outdir,fname), overwrite=True)
            #export_table(self.detections, fname="%s/%s-ap.fits"%(outdir,fname))
        if self.psfcatalogue:
            reindex(self.psfcatalogue)
            hdulist=[fits.PrimaryHDU(header=self.header),fits.BinTableHDU(data=self.psfcatalogue)]
            fits.HDUList(hdulist).writeto("%s/%s-psf.fits"%(outdir,fname), overwrite=True)
            #export_table(self.psfcatalogue, fname="%s/%s-psf.fits"%(outdir,fname))
        if self.background: 
            #self.background.header.update(header)
            self.background.writeto("%s/%s-bgd.fits"%(outdir,fname), overwrite=True)
        if self.residuals is not None:
            im=fits.ImageHDU(data=self.residuals, name="RES", header=self.header)
            im.header.update(self.wcs.to_header())
            im.writeto("%s/%s-res.fits"%(outdir,fname), overwrite=True)
        if self.source_stats is not None:
            reindex(self.source_stats)
            hdulist=[fits.PrimaryHDU(header=self.header),fits.BinTableHDU(data=self.source_stats)]
            fits.HDUList(hdulist).writeto("%s/%s-stat.fits"%(outdir,fname), overwrite=True)

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



    def verify(self):
        """
        This simple function verifies that everything necessary has been loaded properly
        RETURN: 
            0 - on success
            1 - on fail
        """
        status=0
        #warn=lambda :perror(sbold("WARNING: "))

        if self.filter not in starbug2.filters.keys():
            warn()
            perror("Unknown filter '%s'\n"%self.filter)
            status=1

        dname = os.path.expandvars(starbug2.DATDIR)
        if not os.path.exists(dname):
            warn()
            perror("Unable to locate STARBUG_DATDIR='%s'\n"%dname)
            status=1

        else:
            pass
            #if not os.path.exists("%s/%s%s.fits"%(dname, self.filter, self.info["DETECTOR"])):
            #        warn()
            #        perror("Unable to locate filter PSF for '%s'\n"%self.filter)
            #        status=1
        
        if not os.path.exists((dname:=os.path.expandvars(self.options["OUTDIR"]))):
            warn()
            perror("Unable to locate OUTDIR='%s'\n"%dname)
            status=1

        tmp=load_params("%sdefault.param"%pkg_resources.resource_filename("starbug2","param/"))
        if set(tmp.keys()) - set(self.options.keys()):
            warn()
            perror("parameter file version mismatch. Run starbug --update-param to update\n")
            status=1

        return status

    def __getstate__(self):
        state=self.__dict__.copy()
        if "_image" in state:
            del state["_image"] ##Sorry but we cant have that
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        v=self.options["VERBOSE"]
        self.options["VERBOSE"]=0
        self.load_image(self.fname)
        self.options["VERBOSE"]=v
