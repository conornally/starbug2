import starbug2
from starbug2.utils import *
from starbug2.misc import *
from starbug2.routines import *

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter("ignore",category=AstropyWarning)

from astropy.table import hstack

class StarbugBase(object):
    """
    StarbugBase is the overall container for the photometry package. It holds the active image,
    the parameter file and the output images/tables.
    It is self contained enough to simply run "photometry" and everything should just take care 
    of itself from there on.
    """
    filter=None
    fname=None
    detections=None
    psfcatalogue=None
    residuals=[]
    background=None
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

        if self.options["AP_FILE"]: self.load_apfile() ## Load the source list if given
        if self.options["BGD_FILE"]: self.load_bgdfile()
        self.image=self.load_image(fname)   ## Load the fits image

    def info(self):
        """
        Get some useful information from the image header file
        """
        out={}
        keys=("FILTER","DETECTOR")
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
        else: perror("AP_FILE='%s' does not exists\n"%fname)

    def load_bgdfile(self,fname=None):
        """
        Load a BGD_FILE to be used during photometry
        INPUT:  
            fname : file-bgd.fits (this file is exported during background estimation step
        """
        if not fname: fname=self.options["BGD_FILE"]
        if os.path.exists(fname):
            self.background=fits.open(fname)
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
                                        boxsize=int(self.options["BOX_SIZE"]),
                                        filtersize=int(self.options["FILTER_SIZE"]),
                                        verbose=self.options["VERBOSE"])

            dat=detector(self.image["SCI"].data)
            colnames=("RA","DEC","xcentroid","ycentroid","sharpness","roundness1","roundness2","peak")
            dat=dat[colnames]

            self.log("Running Aperture Photometry\n")

            ######################### 
            # Unit Conversion to Jy #
            ######################### 

            if self.image["SCI"].header["BUNIT"]=="MJy/sr":
                factor=1e6*float(self.image["SCI"].header["PIXAR_SR"])
                self.image["SCI"].data*=factor
                self.image["SCI"].header.update({"BUNIT":"Jy"})

                if "ERR" in extnames(self.image) and np.shape(self.image["ERR"]):
                    self.image["ERR"].data*=factor
                    self.image["ERR"].header.update({"BUNIT":"Jy"})

                self.log("-> converting unit from MJy/sr to Jr with factor: %e\n"%factor)

            apphot=APPhot_Routine( self.options["APPHOT_R"], self.options["SKY_RIN"], self.options["SKY_ROUT"], verbose=self.options["VERBOSE"])
            apcorr=apphot.calc_apcorr( self.filter, self.options["APPHOT_R"],"%s/jwst_nircam_apcorr_0004.fits"%(self.options["PSFDIR"] ))
            ap_cat=apphot(dat, self.image, apcorr=apcorr, sig_sky=self.options["SIGSKY"])

            ap_cat.add_column(Column( -2.5*np.log10(ap_cat["flux"] / ZP[self.filter][0]), name="%s_mag"%self.filter))
            ap_cat.add_column(Column( 2.5*np.log10(1+(ap_cat["eflux"]/ap_cat["flux"])), name="%s_emag"%self.filter ))

            self.detections=hstack((dat,ap_cat))
            self.detections.meta=self.info()
            self.detections.meta.update({"STARBUGROUNTINE":"DETECT"})

        else: perror("Failed to run aperture photometry")


    def bgd_estimate(self):
        """
        Estimate the background of the active image
        Saves the result as an ImageHDU self.background
        """
        self.log("Estimating Background\n")
        if self.detections:
            bgd=BackGround_Estimate_Routine(self.detections, 
                                            boxsize=int(self.options["BOX_SIZE"]),
                                            fwhm=starbug2.filters[self.filter][2],
                                            verbose=self.options["VERBOSE"])
            self.background=fits.ImageHDU(data=bgd(self.image[1].data))
        else:
            perror("unable to estimate background, no source list loaded\n")


    def photometry(self):
        """
        Full photometry routine
        Saves the result as a table self.psfcatalogue
        Additionally it appends a residual Image onto the self.residuals HDUList
        """
        if not self.detections:
            perror("unable to run photometry: no source list loaded\n")
            return

        if not self.background:
            perror("unable to run photometry: no background estimation loaded\n")
            return

        if self.image and self.filter:
            self.log("Running PSF Photometry")
            fname=os.path.expandvars("%s/%s.fits"%(self.options["PSFDIR"], self.filter))
            self.log(" <-- %s\n"%fname)
            with fits.open(fname) as fp:
                psf_model=DiscretePRF(fp[0].data)

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
                                    #boxsize=self.options["BOX_SIZE"],
                                    background=self.background[1].data,
                                    #filtersize=self.options["FILTER_SIZE"],
                                    wcs=WCS(self.image[1].header))

            init_guesses=self.detections.copy()
            if "xcentroid" in init_guesses.colnames: init_guesses.rename_column("xcentroid", "x_0")
            if "ycentroid" in init_guesses.colnames: init_guesses.rename_column("ycentroid", "y_0")

            init_guesses.remove_columns(("sharpness", "roundness1", "roundness2", "npix", "sky", "flux", "mag", "ap_flux", "ap_sky"))

            self.psfcatalogue=tabppend(self.psfcatalogue, phot(self.image[1].data, init_guesses=init_guesses))
            self.residuals.append(phot.get_residual_image())
            self.background=fits.ImageHDU(data=phot.bkg_estimator.bgd, name="BACKGROUND")

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

        head=fits.Header(self.options)
        dname,fname,ext=split_fname(self.fname)
        if self.detections: 
            hdulist=[fits.PrimaryHDU(header=head),fits.BinTableHDU(data=self.detections)]
            fits.HDUList(hdulist).writeto("%s/%s-ap.fits"%(outdir,fname), overwrite=True)
            #export_table(self.detections, fname="%s/%s-ap.fits"%(outdir,fname))
        if self.psfcatalogue: 
            hdulist=[fits.PrimaryHDU(header=head),fits.BinTableHDU(data=self.psfcatalogue)]
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

        if self.filter not in starbug2.filters.keys():
            perror("WARNING: Unknown filter '%s'\n"%self.filter)
            status=1

        dname = os.path.expandvars(self.options["PSFDIR"])
        if not os.path.exists(dname):
            perror("WARNING: Unable to locate PSFDIR='%s'\n"%dname)
            status=1

        else:
            if not os.path.exists("%s/%s.fits"%(dname, self.filter)):
                perror("WARNING: Unable to locate filter PSF for '%s'\n"%self.filter)
                status=1
        
        if not os.path.exists((dname:=os.path.expandvars(self.options["OUTDIR"]))):
            perror("WARNING: Unable to locate OUTDIR='%s'\n"%dname)
            status=1

        return status

    def __getstate__(self):
        state=self.__dict__.copy()
        if "image" in state:
            del state["image"] ##Sorry but we cant have that
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.load_image(self.fname)
