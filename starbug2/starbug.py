import starbug2
from starbug2.utils import *
from starbug2.misc import *
from starbug2.routines import *

from astropy.table import hstack

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
    residuals=[]
    background=None
    image=None
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
        #if self.options["TEST"]: print(self.options["TEST"])

        self.image=self.load_image(fname)   ## Load the fits image
        if self.options["AP_FILE"]: self.load_apfile() ## Load the source list if given
        if self.options["BGD_FILE"]: self.load_bgdfile()

    @property
    def header(self):
        """
        """
        head=fits.Header({**self.options,**self.info})
        head["CALIBLEVEL"]=self.stage
        head["STARBUG"]=pkg_resources.get_distribution("starbug2").version
        head["FILTER"]=self.filter
        return head


    @property
    def info(self):
        """
        Get some useful information from the image header file
        """
        out={}
        keys=("FILTER","DETECTOR","TELESCOP","INSTRUME")
        if self.image:
            out.update( { (key,self.image[0].header[key]) for key in keys if key in self.image[0].header})
        return out

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
        image=None
        self.fname=fname
        if fname:
            dname,bname,extension=split_fname(fname)
            if extension==".fits":
                if os.path.exists(fname):
                    image=fits.open(fname)
                    #self.output=dname

                    if ("FILTER" in image[0].header) and (image[0].header["FILTER"] in starbug2.filters.keys()):
                        self.filter=image[0].header["FILTER"]
                    else: perror("Unable to dtermine image filter\n")
                        
                    ## I NEED TO DETERMINE BETTER WHAT STAGE IT IS IN
                    exts=extnames(image)
                    hdr=image[0].header

                    if "SCI"  in exts: self.wcs=WCS(image["SCI"].header)
                    else: self.wcs=WCS(hdr)


                    if "DQ" in exts:
                        if "AREA" in exts: self.stage=2
                        else: self.stage=2.5
                    elif "WHT" in exts: self.stage=3 
                    else: perror("unable to determine jwst pipeline level\n")


                    self.log("loaded: \"%s\"\n"%fname)
                    self.log("pipeline stage: %d\n"%self.stage)

                else: perror("fits file \"%s\" does not exist\n"%fname)
            else: perror("included file must be FITS format\n")
        return image

    def load_apfile(self,fname=None):
        """
        Load a AP_FILE to be used during photometry
        INPUT:  
            fname : file-ap.fits (this file is exported during source detection step
        """
        if not fname: fname=self.options["AP_FILE"]
        if os.path.exists(fname):
            with fits.open(fname) as fp:
                self.detections=Table(data=fp[1].data._get_raw_data())
                self.log("loaded AP_FILE='%s'\n"%fname)

                cn=self.detections.colnames
                if not any( _ in cn for _ in ("xcentroid","ycentroid","x_0","y_0")):
                    if all( _ in cn for _ in ("RA","DEC")):
                        xy=self.wcs.all_world2pix(self.detections["RA"], self.detections["DEC"],0)
                        print(xy)
                        self.detections.add_columns(xy,names=("xcentroid","ycentroid"),indexes=[0,0])
                    else: perror("WARNING, unable to determine physical coordinates from detections table\n")
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
        else: perror("BGD_FILE='%s' does not exists\n"%fname)

    def detect(self):
        """
        Full source detection routine
        Saves the result as a table self.detections
        """
        self.log("Detecting Sources\n")
        if self.image and self.filter:
            FWHM=starbug2.filters[self.filter][2]
            detector=Detection_Routine( sig_src=self.options["SIGSRC"],
                                        sig_sky=self.options["SIGSKY"],
                                        fwhm=FWHM,
                                        sharplo=self.options["SHARP_LO"],
                                        sharphi=self.options["SHARP_HI"],
                                        roundlo=self.options["ROUND_LO"],
                                        roundhi=self.options["ROUND_HI"],
                                        wcs=WCS(self.image[1].header),
                                        bgd2d=self.options["DOBGD2D"],
                                        boxsize=int(self.options["BOX_SIZE"]),
                                        filtersize=int(self.options["FILTER_SIZE"]),
                                        verbose=self.options["VERBOSE"])

            dat=detector(self.image["SCI"].data)
            colnames=("RA","DEC","xcentroid","ycentroid","sharpness","roundness1","roundness2", "peak")
            self.detections=dat[colnames]
            self.aperture_photometry()


    def aperture_photometry(self):

        if self.detections is None:
            perror("No detection source file loaded (-d file-ap.fits)\n")
            return
        colnames=list( name for name in ("RA","DEC","xcentroid","ycentroid","sharpness","roundness1","roundness2", "peak") if name in self.detections.colnames)
        dat=self.detections[colnames]
        #dat=self.detections[("RA","DEC","xcentroid","ycentroid","sharpness","roundness1","roundness2", "peak")]
        #######################
        # APERTURE PHOTOMETRY #
        #######################
        self.log("Running Aperture Photometry\n")
        image=self.image["SCI"].data.copy() ##dont work on the real image!

        ######################### 
        # Unit Conversion to Jy #
        ######################### 
        error=None
        scalefactor=1
        if self.image["SCI"].header["BUNIT"]=="MJy/sr":
            scalefactor=1e6*float(self.image["SCI"].header["PIXAR_SR"])
            self.log("-> converting unit from MJy/sr to Jr with factor: %e\n"%scalefactor)

        image*=scalefactor
        if "ERR" in extnames(self.image) and np.shape(self.image["ERR"]):
            error=self.image["ERR"].data
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

        if   self.info["INSTRUME"]=="NIRCAM": fname="%s/apcorr_nircam.fits"%self.options["PSFDIR"]
        elif self.info["INSTRUME"]=="MIRI":   fname="%s/apcorr_miri.fits"%self.options["PSFDIR"]
        else: perror("No apcorr file available for instrument\n")

        if self.options["FIT_APP_R"]:
            apcorr=APPhot_Routine.calc_apcorr(self.filter, self.options["APPHOT_R"], table_fname=fname, verbose=self.options["VERBOSE"])
        else:
            apcorr,radius=APPhot_Routine.apcorr_from_encenergy(self.filter,self.options["ENCENERGY"],table_fname=fname, verbose=self.options["VERBOSE"])

        apphot=APPhot_Routine( radius, skyin, skyout, encircled_energy=self.options["ENCENERGY"], fit_radius=self.options["FIT_APP_R"], verbose=self.options["VERBOSE"])

        if self.stage==2:
            image*= self.image["AREA"].data ## AREA distortion correction
            mask=self.image["DQ"].data & (DQ_DO_NOT_USE|DQ_SATURATED) #|DQ_JUMP_DET)
            image[mask]=np.nan
            error[mask]=np.nan
            ap_cat=apphot(dat, image, error=error, dqflags=self.image["DQ"].data, apcorr=apcorr, sig_sky=self.options["SIGSKY"])

        else: ##stage 3 version
            ap_cat=apphot(dat, image, error=error, apcorr=apcorr, sig_sky=self.options["SIGSKY"])

        mag,magerr=flux2ABmag( ap_cat["flux"], ap_cat["eflux"], filter=self.filter)
        ap_cat.add_column(Column(mag,self.filter))
        ap_cat.add_column(Column(magerr,"e%s"%self.filter))

        self.detections=hstack((dat,ap_cat))
        self.detections.meta=dict(self.header.items())
        self.detections.meta.update({"ROUNTINE":"DETECT"})

        #else: perror("Failed to run aperture photometry")

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
            sources=sources[ sources[xname]<self.image["SCI"].header["NAXIS1"]]
            sources=sources[ sources[yname]<self.image["SCI"].header["NAXIS2"]]
            bgd=BackGround_Estimate_Routine(sources, 
                                            boxsize=int(self.options["BOX_SIZE"]),
                                            fwhm=starbug2.filters[self.filter][2],
                                            verbose=self.options["VERBOSE"])
            self.background=fits.ImageHDU(data=bgd(self.image[1].data))
        else:
            perror("unable to estimate background, no source list loaded\n")

    def bgd_subtraction(self):
        
        self.log("Subtracting Background\n")

        if self.background is None:
            perror("No background array loaded (-b file-bgd.fits)\n")
            return 
        array= self.image["SCI"].data - self.background.data
        self.residuals.append(array)
        self.image["SCI"].data=array

    def photometry(self):
        """
        Full photometry routine
        Saves the result as a table self.psfcatalogue
        Additionally it appends a residual Image onto the self.residuals HDUList
        """
        if self.detections is None:
            perror("unable to run photometry: no source list loaded\n")
            return

        if self.background is None:
            perror("unable to run photometry: no background estimation loaded\n")
            return

        if self.image and self.filter:
            self.log("Running PSF Photometry")
            fname=os.path.expandvars("%s/%s.fits"%(self.options["PSFDIR"], self.filter))
            self.log(" <-- %s\n"%fname)
            with fits.open(fname) as fp:
                psf_model=DiscretePRF(fp[0].data)

                print(psf_model.fixed)
                psf_model.fixed["x_0"]=False
                psf_model.fixed["y_0"]=False
            phot=PSFPhot_Routine(   self.options["CRIT_SEP"],
                                    starbug2.filters[self.filter][2],
                                    psf_model,
                                    psf_model.prf_shape,
                                    sig_sky=self.options["SIGSKY"],
                                    sig_src=self.options["SIGSRC"],
                                    sharplo=self.options["SHARP_LO"],
                                    sharphi=self.options["SHARP_HI"],
                                    roundlo=self.options["ROUND_LO"],
                                    roundhi=self.options["ROUND_HI"],
                                    background=self.background.data,
                                    wcs=self.wcs)

            init_guesses=self.detections.copy()

            if "xcentroid" in init_guesses.colnames: init_guesses.rename_column("xcentroid", "x_0")
            if "ycentroid" in init_guesses.colnames: init_guesses.rename_column("ycentroid", "y_0")

            init_guesses=init_guesses[ init_guesses["x_0"]>=0 ]
            init_guesses=init_guesses[ init_guesses["y_0"]>=0 ]
            init_guesses=init_guesses[ init_guesses["x_0"]<self.image["SCI"].header["NAXIS1"]]
            init_guesses=init_guesses[ init_guesses["y_0"]<self.image["SCI"].header["NAXIS2"]]


            init_guesses=init_guesses[["x_0","y_0"]]

            with open("out.reg","w") as fp:
                for line in init_guesses:
                    fp.write("circle %f %f 3\n"%( line["x_0"]+1, line["y_0"]+1))
            self.psfcatalogue=tabppend(self.psfcatalogue, phot(self.image[1].data, init_guesses=init_guesses))
            self.residuals.append(phot.get_residual_image())
            self.background=fits.ImageHDU(data=phot.bkg_estimator.bgd, name="BACKGROUND") ##So is it supposed to be a fits image or a numpy array?!

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

        fname=os.path.expandvars("%s/%s.fits"%(self.options["PSFDIR"], self.filter))
        with fits.open(fname) as fp:
            psf_model=DiscretePRF(fp[0].data)

        detector=Detection_Routine( sig_src=self.options["SIGSRC"],
                                    sig_sky=self.options["SIGSKY"],
                                    fwhm=starbug2.filters[self.filter][2],
                                    sharplo=self.options["SHARP_LO"],
                                    sharphi=self.options["SHARP_HI"],
                                    roundlo=self.options["ROUND_LO"],
                                    roundhi=self.options["ROUND_HI"],
                                    wcs=WCS(self.image[1].header),
                                    verbose=0)

        phot=PSFPhot_Routine(   self.options["CRIT_SEP"],
                                starbug2.filters[self.filter][2],
                                psf_model,
                                psf_model.prf_shape,
                                sig_sky=self.options["SIGSKY"],
                                sig_src=self.options["SIGSRC"],
                                sharplo=self.options["SHARP_LO"],
                                sharphi=self.options["SHARP_HI"],
                                roundlo=self.options["ROUND_LO"],
                                roundhi=self.options["ROUND_HI"],
                                wcs=WCS(self.image[1].header),
                                verbose=0)

        art=ArtificialStar_Routine(detector, phot, psf_model)
        self.log("Artificial Star Testing (n=%d)\n"%(self.options["NUMBER_ARTIFICIAL_STARS"]))
        result=art.run(self.image[1].data, ntests=self.options["NUMBER_ARTIFICIAL_STARS"], flux_range=(self.options["MIN_FLUX"], self.options["MAX_FLUX"]),
                subimage_size=self.options["SUBIMAGE_SIZE"], separation_thresh=self.options["SEPARATION_THRESH"], fwhm=starbug2.filters[self.filter][2])
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

        #hdulist=[fits.PrimaryHDU(header=header)]
        hdulist=[fits.PrimaryHDU()]
        for n,res in enumerate(self.residuals,1):
            im=fits.ImageHDU(data=res,name="RESIDUAL%d"%n)
            #im.header.update(header)
            hdulist.append(im)

        if len(hdulist)>1: 
            fits.HDUList(hdulist).writeto("%s/%s-res.fits"%(outdir,fname), overwrite=True)



    def verify(self):
        """
        This simple function verifies that everything necessary has been loaded properly
        RETURN: 
            0 - on success
            1 - on fail
        """
        status=0
        warn=lambda :perror(sbold("WARNING: "))

        if self.filter not in starbug2.filters.keys():
            warn()
            perror("Unknown filter '%s'\n"%self.filter)
            status=1

        dname = os.path.expandvars(self.options["PSFDIR"])
        if not os.path.exists(dname):
            warn()
            perror("Unable to locate PSFDIR='%s'\n"%dname)
            status=1

        else:
            if not os.path.exists("%s/%s.fits"%(dname, self.filter)):
                warn()
                perror("Unable to locate filter PSF for '%s'\n"%self.filter)
                status=1
        
        if not os.path.exists((dname:=os.path.expandvars(self.options["OUTDIR"]))):
            warn()
            perror("Unable to locate OUTDIR='%s'\n"%dname)
            status=1

        tmp=load_params("%sdefault.param"%pkg_resources.resource_filename("starbug2","param/"))
        if set(tmp.keys()) - set(self.options.keys()):
            warn()
            perror("parameter file version mismatch. Run starbug --local-param to update\n")
            status=1

        return status

    def __getstate__(self):
        state=self.__dict__.copy()
        if "image" in state:
            del state["image"] ##Sorry but we cant have that
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        v=self.options["VERBOSE"]
        self.options["VERBOSE"]=0
        self.load_image(self.fname)
        self.options["VERBOSE"]=v
