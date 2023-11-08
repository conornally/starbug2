"""
Routines for Starbug
"""
import os
import sys
import time
import numpy as np
import pkg_resources
from scipy.stats import norm#, mode
from scipy.optimize import curve_fit
from scipy.ndimage import convolve
from skimage.feature import match_template

from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.coordinates import SkyCoord
from astropy.table import Column, Table, QTable, hstack, vstack
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.convolution import RickerWavelet2DKernel

from photutils.background import Background2D, BackgroundBase
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.detection import StarFinderBase, DAOStarFinder, find_peaks
from photutils.psf.groupstars import DAOGroup
from photutils.psf import PSFPhotometry, IntegratedGaussianPRF, SourceGrouper

from photutils.datasets import make_model_sources_image, make_random_models_table
from starbug2.utils import loading, printf, perror, warn
from starbug2 import *

class Detection_Routine(StarFinderBase):
    """
    Detection routine- called by starbug
    A standalone detection that runs on a 2D image.
    It uses DAOStarFinder as the base for peak detection but run
    several times on a series of background subtracted images.
    Each run the background subtraction is differemt, bringing out a
    different set of sources
    """
    def __init__(self,  sig_src=5, sig_sky=3, fwhm=1,
                        sharplo=0.2, sharphi=1, round1hi=1, round2hi=1,
                        smoothlo=-np.inf, smoothhi=np.inf, ricker_r=1.0,
                        verbose=0, cleansrc=1, dobgd2d=1, boxsize=2, doconvl=1):
        self.sig_src=sig_src
        self.sig_sky=sig_sky
        self.fwhm = fwhm
        self.sharphi=sharphi
        self.sharplo=sharplo
        self.round1hi= round1hi if round1hi is not None else np.inf
        self.round2hi= round2hi if round2hi is not None else np.inf
        self.smoothlo= smoothlo if smoothlo is not None else -np.inf
        self.smoothhi= smoothhi if smoothhi is not None else np.inf

        self.ricker_r=ricker_r
        self.cleansrc=cleansrc

        self.catalogue=Table()
        self.verbose=verbose

        self.dobgd2d=dobgd2d
        self.boxsize=boxsize
        self.doconvl=doconvl

    def detect(self, data, bkg_estimator=None, xycoords=None, method=None):
        """
        The core detection step (DAOStarFinder
        INPUT:  
                data=array to detect on
                bkg_estimator=background array the same shape as data array
                xycorrds= table of initial guesses (xcentroid,ycentroid)
        RETURNS:
                Sourcelist Table
        """
        bkg=np.zeros(data.shape)
        if bkg_estimator:
            bkg=bkg_estimator(data)

        _,median,std=sigma_clipped_stats(data,sigma=self.sig_sky)
        if method=="findpeaks":
            return find_peaks(data-bkg, median+std*self.sig_src, box_size=11)

        else:
            roundhi=max((self.round1hi,self.round2hi))
            find=DAOStarFinder(std*self.sig_src, self.fwhm, sharplo=self.sharplo, sharphi=self.sharphi,
                roundlo=-roundhi, roundhi=roundhi, peakmax=np.inf, xycoords=xycoords)
            return find(data - bkg)


    def _bkg2d(self, data):
        return Background2D(data, self.boxsize, filter_size=3).background

    def match(self, base, cat):
        """
        Internal function to class
        Used to match detenctoins from separate background subtracted images
        into the main catalogue. This will append a source if its matched separation
        is above the threshold self.match_threshold
        """

        added=0
        base_sky=SkyCoord(x=base["xcentroid"], y=base["ycentroid"], z=np.zeros(len(base)), representation_type="cartesian")
        cat_sky=SkyCoord(x=cat["xcentroid"], y=cat["ycentroid"], z=np.zeros(len(cat)), representation_type="cartesian")
        _,separation,_=cat_sky.match_to_catalog_3d(base_sky)
        mask= separation.to_value() >self.fwhm
        return vstack((base,cat[mask]))


    def find_stars(self, data, mask=None):
        """
        Main function of Routine
        FUNC:
            This routine runs source detection several times, but on a different form
            of the data array each time. Each form has been "skewed" somehow to brighten the
            most faint sources and flatten the differential background.
            1- Plain detections
            2- Subtract Background estimation
            3- RickerWave convolution

        INPUT:
            data=array to detect on
            mask=pixels to mask out
        """
        if data is None: return None
        if mask is None: mask=np.where(np.isnan(data))
        _,median,_=sigma_clipped_stats(data,sigma=self.sig_sky)
        data[mask]=median

        self.catalogue=self.detect(data)
        if self.verbose: printf("-> [PLAIN] pass: %d sources\n"%len(self.catalogue))

        if self.dobgd2d:
            self.catalogue=self.match(self.catalogue, self.detect(data, self._bkg2d))
            if self.verbose: printf("-> [BGD2D] pass: %d sources\n"%len(self.catalogue))

        ## 2nd order differential detection
        if self.doconvl:
            kernel=RickerWavelet2DKernel(self.ricker_r)
            conv=convolve(data, kernel)
            corr=match_template(conv/np.amax(conv), kernel.array)
            _detections=self.detect(corr, method="findpeaks")
            if _detections:
                _detections["x_peak"]+=kernel.shape[0]//2
                _detections["y_peak"]+=kernel.shape[0]//2
                _detections.rename_columns( ("x_peak","y_peak"),("xcentroid","ycentroid"))
                self.catalogue=self.match(self.catalogue, _detections)
            if self.verbose: printf("-> [CONVL] pass: %d sources\n"%len(self.catalogue))

        ## Now with xycoords DAOStarfinder will refit the sharp and round values at the detected locations
        #self.catalogue=self.detect(data, xycoords=np.array([self.catalogue["xcentroid"],self.catalogue["ycentroid"]]).T)#, clean=0)
        tmp=SourceProperties(data,self.catalogue, verbose=self.verbose).calculate_geometry(self.fwhm)
        if tmp: self.catalogue=tmp

        mask=(~np.isnan(self.catalogue["xcentroid"]) & ~np.isnan(self.catalogue["ycentroid"]))
        #self.catalogue.remove_rows(~mask)

        if self.cleansrc: 
            mask &=((self.catalogue["sharpness"]>self.sharplo)
                   &(self.catalogue["sharpness"]<self.sharphi)
                   &(self.catalogue["roundness1"]> -self.round1hi)
                   &(self.catalogue["roundness1"]<  self.round1hi)
                   &(self.catalogue["roundness2"]> -self.round2hi)
                   &(self.catalogue["roundness2"]<  self.round2hi))
        if self.verbose: printf("-> cleaning %d unlikley point sources\n"%sum(~mask))
        self.catalogue.remove_rows(~mask)
        
        if self.verbose: printf("Total: %d sources\n"%len(self.catalogue))

        self.catalogue.replace_column("id", Column(range(1,1+len(self.catalogue))))

        return self.catalogue

class APPhot_Routine():
    """
    Aperture photometry called by starbug
    Given photometry radius, sky annuli radii rad_inner rad_outer
    """
    def __init__(self, radius, sky_in, sky_out, verbose=0):
        if sky_in < radius:
            warn()
            perror("Sky annulus radii must be larger than aperture radius.\n")
            sky_in=radius+1

        if sky_in >= sky_out: 
            warn()
            perror("Sky annulus outer radii must be larger than the inner.\n")
            sky_out=sky_in+1

        self.radius=radius
        self.sky_in=sky_in
        self.sky_out=sky_out
        self.catalogue=Table(None)#, names=["ap_flux_r%d"%n for n in range(len(radii))]+["sky_median"])
        self.verbose=verbose

    def __call__(self, detections, image, **kwargs):
        return self.run(detections, image, **kwargs)

    def run(self, detections, image, error=1, dqflags=None, apcorr=1.0, sig_sky=3):
        """
        Forced aperture photometry on a list of detections
        detections are a astropy.table.Table with columns xcentroid ycentroid or x_0 y_0
        This will add extra columns into this table ap_flux ap_sky
        INPUT:  detections - Astropy.Table containing source list
                image       -2D image array to run photometry on
                error       -2D image array OR scalar to act as photometric error
                dqflags     -2D array of data quality flags where appropriate

        RETURN: Photometry catalogue
        """
        if len( set(("xcentroid","ycentroid")) & set(detections.colnames))==2:
            pos=[(line["xcentroid"],line["ycentroid"]) for line in detections]
        elif len( set(("x_0","y_0")) & set(detections.colnames))==2:
            pos=[(line["x_0"],line["y_0"]) for line in detections]
        else:
            perror("Cannot identify position in detection catalogue (x_0/xcentroid)\n");
            return None

        mask=np.isnan(image)

        apertures=CircularAperture(pos,self.radius)
        smooth_apertures=CircularAperture(pos, min(1.5*self.radius,self.sky_in))
        annulus_aperture=CircularAnnulus(pos, r_in=self.sky_in, r_out=self.sky_out)
        
        self.log("-> apertures: %.2g (%.2g - %.2g)\n"%(self.radius, self.sky_in, self.sky_out))
        phot=aperture_photometry(image, (apertures,smooth_apertures), error=error, mask=mask)
        self.catalogue=Table(np.full((len(pos),4),np.nan),names=("smoothness","flux","eflux","sky"))

        self.log("-> calculating sky values\n")
        masks=annulus_aperture.to_mask(method="center")
        dat=list(map(lambda a:a.multiply(image),masks))

        try: dat=np.array(dat).astype(float)
        except:
            ## Cases where the array is inhomegenoeus
            ## If annulus reaches the edge of the image, it will create a mask the wrong shape
            ## If for whatever reason the point lies outside the image, it will have None
            ## in the list, this needs to be caught too
            warn()
            perror("Ran into issues with the sky annuli, trying to fix them..\n")
            size=np.max( [np.shape(d) for d in dat if d is not None ])
            for i,d in enumerate(dat):
                if d is None: dat[i]=np.zeros((size,size))
                elif (shape:=np.shape(d))!=(size,size):
                    dat[i]=np.zeros((size,size))
                    dat[i][:shape[0],:shape[1]]+=d
            dat=np.array(dat)

        mask=(dat>0 & np.isfinite(dat))
        dat[~mask]=np.nan
        dat=sigma_clip(dat.reshape(dat.shape[0],-1), sigma=sig_sky,axis=1)
        self.catalogue["sky"]=np.ma.median(dat,axis=1)
        std=np.ma.std(dat,axis=1)

        epoisson=phot["aperture_sum_err_0"]
        esky_scatter= apertures.area*std**2
        esky_mean=  (std**2 * apertures.area**2) / annulus_aperture.area

        self.catalogue["eflux"]=np.sqrt( epoisson**2 +esky_scatter**2 +esky_mean**2)
        self.catalogue["flux"]=apcorr*(phot["aperture_sum_0"] - (self.catalogue["sky"]*apertures.area))
        
        ######################
        # Source "smoothness", the gradient of median pixel values within the two test apertures
        ######################
        self.catalogue["smoothness"] = (phot["aperture_sum_1"]/smooth_apertures.area) / (phot["aperture_sum_0"]/apertures.area)

        col=Column(np.full(len(apertures),SRC_GOOD), dtype=np.uint16, name="flag")
        if dqflags is not None:
            self.log("-> flagging unlikely sources\n")
            for i, mask in enumerate(apertures.to_mask(method="center")):
                _tmp=mask.multiply(dqflags)
                if _tmp is not None:
                    dat=np.array(_tmp,dtype=np.uint32)
                    if np.sum( dat & (DQ_DO_NOT_USE|DQ_SATURATED)): col[i]|=SRC_BAD
                    if np.sum( dat & DQ_JUMP_DET): col[i]|=SRC_JMP
                #else: col[i]|=SRC_UKN
        self.catalogue.add_column(col)
        return self.catalogue


    @staticmethod
    def calc_apcorr(filter, radius, table_fname=None, verbose=0):
        """
        Using CRDS apcorr table, fit a curve to the radius vs apcorr
        columns and then return aporr to respective input radius
        """
        if not table_fname or not os.path.exists(table_fname): return 1
        tmp=Table.read(table_fname, format="fits")
        
        if "filter" in tmp.colnames:
            t_apcorr=tmp[(tmp["filter"]==filter)]
        else: t_apcorr=tmp


        if "pupil" in t_apcorr.colnames:
            t_apcorr=t_apcorr[ t_apcorr["pupil"]=="CLEAR"]
        
        apcorr= np.interp(radius, t_apcorr["radius"], t_apcorr["apcorr"])
        if verbose: printf("-> estimating aperture correction: %.3g\n"%apcorr)

        #eefrac= np.interp(radius, t_apcorr["radius"], t_apcorr["eefraction"])
        #if verbose: printf("-> effective encircled energy fraction: %.3g\n"%eefrac)
        return apcorr


    @staticmethod
    def apcorr_from_encenergy(filter, encircled_energy, table_fname=None, verbose=0):
        """
        Rather than fitting radius to the APCORR CRDS, use the closes Encircled energy value
        """
        if not table_fname or not os.path.exists(table_fname): return 1
        tmp=Table.read(table_fname, format="fits")

        if "filter" in tmp.colnames:
            t_apcorr=tmp[(tmp["filter"]==filter)]
        else: t_apcorr=tmp

        line=t_apcorr[(np.abs(t_apcorr["eefraction"]-encircled_energy)).argmin()]
        if verbose:
            printf("-> best matching encircled energy %.1f, with radius %g pixels\n"%(line["eefraction"],line["radius"]))
            printf("-> using aperture correction: %f\n"%line["apcorr"])

        return line["apcorr"], line["radius"]
    
    def radius_from_encenrgy(filter, eefrac, table_fname):
        """
        """
        if not table_fname or not os.path.exists(table_fname): return -1
        t_apcorr=Table.read(table_fname, format="fits")

        if len( set(["eefraction","radius"])&set(t_apcorr.colnames))!=2: return -1

        if "filter" in t_apcorr.colnames: # Crop down table
            t_apcorr=t_apcorr[(t_apcorr["filter"]==filter)]

        if "pupil" in t_apcorr.colnames: # Crop down table
            t_apcorr=t_apcorr[ t_apcorr["pupil"]=="CLEAR"]

        return np.interp( eefrac, t_apcorr["eefraction"], t_apcorr["radius"])

    def log(self,msg):
        """
        log message if in verbose mode
        """
        if self.verbose: 
            printf(msg)
            sys.stdout.flush()


class BackGround_Estimate_Routine(BackgroundBase):
    """
    """
    def __init__(self, sourcelist, boxsize=2, fwhm=2, sigsky=2, bgd_r=-1, verbose=0, bgd=None):#mask_r0=7, mask_r1=9
        self.sourcelist=sourcelist
        self.boxsize=boxsize
        self.fwhm=fwhm
        self.sigsky=sigsky
        self.bgd_r=bgd_r
        self.verbose=verbose
        self.bgd=bgd
        super().__init__()

    def calc_peaks(self,im):
        """
        Determine peak pixel value for each source in xy
        """
        x=self.sourcelist["xcentroid"]
        y=self.sourcelist["ycentroid"]
        apertures=CircularAperture(np.array((x,y)).T,2).to_mask()
        peaks=np.full(len(x),np.nan)
        for i,mask in enumerate(apertures):
            peaks[i]=np.nanmax( mask.multiply(im) )
        return peaks
    
    def log(self,msg):
        if self.verbose: printf(msg)


    """
    def calc_rlist(self,data):

        _NCALC_BGDR=5.0
        FUDGE=1/2.0
        N=int(len(self.sourcelist)/_NCALC_BGDR)
        ii=np.random.choice( len(self.sourcelist), size=N)
        sources=self.sourcelist[ii]
        shape=int(self.fwhm)*5
        if not shape%2: shape+=1

        mask=np.isnan(data)|np.isinf(data)
        _,median,_=sigma_clipped_stats(data, sigma=3)
        ro= 0.5*self.fwhm * np.sqrt( np.log(median/sources["flux"]))
        ro*=FUDGE

        with open("test.reg",'w') as fp:
            for i in range(len(sources)):
                if not np.isnan(ro[i]):
                    src=sources[i]
                    fp.write("fk5;circle %f %f %fi\n"%( src["RA"], src["DEC"], ro[i]))
        x=np.array([min(sources["flux"]), max(sources["flux"])])
    """

    def __call__(self, data, axis=None, masked=False):
        if self.sourcelist is None or data is None: return self.bgd
        _data=np.copy(data)
        X,Y=np.ogrid[:data.shape[1], :data.shape[0]]

        FUDGE=1
        DEFAULT_R=2*self.fwhm

        if self.bgd_r and self.bgd_r>0:
            self.log("-> using BGD_R=%g masking aperture radii\n"%self.bgd_r)
            rlist=self.bgd_r*np.ones(len(self.sourcelist)) 

        else:
            #peaks=self.calc_peaks(data)
            #rlist=np.sqrt(peaks**0.7)*self.fwhm/1.5 ## <-- that works but hmm
            if "flux" in self.sourcelist.colnames:
                self.log("-> calculating source aperture mask radii\n")
                _,median,_=sigma_clipped_stats(data,sigma=self.sigsky)
                rlist= FUDGE * self.fwhm/2 * np.sqrt( np.log( median/self.sourcelist["flux"]) )
                rlist[np.isnan(rlist)]=DEFAULT_R
            else:
                warn()
                perror("Unable to caluclate aperture mask sizes, add '-A' to starbug command.\n")
                rlist=DEFAULT_R*np.ones(len(self.sourcelist))

        D=50
        load=loading(len(self.sourcelist), msg="masking sources", res=10)#len(self.sourcelist)/1000)
        for r,src in zip(rlist,self.sourcelist):

            rin=1.5*r
            rout=rin+1

            x=round(src["xcentroid"])
            y=round(src["ycentroid"])
            _X=X[ max( x-D,0):min(x+D,data.shape[1])]
            _Y=Y[ :,max( y-D,0):min(y+D,data.shape[0])]

            R=np.sqrt((_X-src["xcentroid"])**2+(_Y-src["ycentroid"])**2)

            mask=(R<r)
            annuli_mask=( (R>rin) & (R<rout) )

            tmp=_data[_Y,_X]
            tmp[mask]=np.median(data[_Y,_X][annuli_mask])
            _data[_Y,_X]=tmp

            load()
            if self.verbose: load.show() ## This will slow the thing down quite a lot
        if self.verbose: printf("-> estimating bgd2d\n")
        self.bgd=Background2D(_data, self.boxsize).background
        return self.bgd

    def calc_background(self,data, axis=None, masked=None):
        if self.bgd is None: self.__call__(data)
        return self.bgd

class _grouping(DAOGroup):
    """
    Overwritten DAOGroup that just holds the number of groups
    for use in verbose loading of psfphot routine
    >>> This is now a bit redundant after photoutils added progress_bar=true

    Issue:6 recursion during fitting.
        >>> the fitting seems to include this recursive step within a given source of sources
            if the group contains more members than the system recursion limit then it will
            unceremonially crash on a recursion error.
        >>> for now at least, I will give a warning that this would have occurred, but then 
            ill override the recursion limit and avoid the crash. Hopefully
    """
    logfile=None
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ngroups=0
        if os.getenv("SBIIDEBUG"):
            self.logfile=open(os.getenv("SBIIDEBUG"),'a')
            self.logfile.write("# PSF Source Grouping\n")

    def __call__(self, *args):
        res=super().__call__(*args)
        self.ngroups=max(res["group_id"])

        #### hacking recursion error
        for gid in set(res["group_id"]):
            n_members=sum(res["group_id"]==gid)
            if self.logfile:
                self.logfile.write("GID:%d (%d)\n"%(gid,n_members))

            ## It seems to not quite hit the recursion limit
            ## Crashed on 980 with a limit of 1000, so im going to try 90%
            if n_members > (0.9*sys.getrecursionlimit()):
                warn()
                perror("This run will exceed the recursion depth of the system. "
                       "Starbug will intervene and override the recursion limit but "
                       "the parameter \"CRIT_SEP\" should be reduced to avoid this.\n"
                       "Setting recursion limit %d -> %d\n"%(sys.getrecursionlimit(), int(2.0*n_members)))
                sys.setrecursionlimit(int(2.0*n_members))
        if self.logfile: self.logfile.close()
        return res

class _fitmodel(LevMarLSQFitter):
    load=None
    def __init__(self, grouper=None, verbose=1):
        super().__init__()
        self.grouper=grouper
        if verbose:
            self.load=loading(1, msg="fitting psfs")

    def __call__(self, *args, **kwargs):
        if self.grouper and self.load:
            self.load.setlen(self.grouper.ngroups)
            if self.load is not None:
                self.load()
                self.load.show()
        return super().__call__(*args,**kwargs)


class PSFPhot_Routine(PSFPhotometry):
    """
    PSF Photometry routine called by starbug
    """
    def __init__(self, psf_model, fitshape, apphot_r=3,
            force_fit=False, background=None, verbose=1):

        self.verbose=verbose
        self.force_fit=force_fit        
        self.background=background

        #group_maker=_grouping(crit_separation=8)#crit_separation)
        #bkg_estimator=BackGround_Estimate_Routine(None, bgd=background)
        #fitter=_fitmodel(grouper=group_maker, verbose=verbose)
        
        if force_fit:
            psf_model.x_0.fixed=True
            psf_model.y_0.fixed=True
            #psf_model.fixed.update( {"x_0":True,"y_0":True} )
            #if fitter.load is not None: fitter.load.msg+=" (forced)"

        super().__init__(psf_model=psf_model, fit_shape=fitshape, finder=None,
                progress_bar=verbose, aperture_radius=apphot_r,
                grouper=SourceGrouper(8))
                #localbkg_estimator=bkg_estimator, fitter=fitter)
        """
        super().__init__(grouper=group_maker, localbkg_estimator=bkg_estimator,
                psf_model=psf_model, fit_shape=fitshape,
                finder=None, fitter=fitter)
        """

    def __call__(self,*args,**kwargs): return self.do_photometry(*args,**kwargs)
    def do_photometry(self, image, mask=None, init_params=None, progress_bar=False):
        """
        """ 

        if init_params is None or len(init_params)==0:
            perror("Must include source list\n")
            return None

        if self.background is not None: image=image-self.background
        if self.verbose: printf("-> fitting %d sources\n"%len(init_params))
        cat=super().__call__(image, mask=mask, init_params=init_params)

        d=np.sqrt((cat["x_init"]-cat["x_fit"])**2.0 + (cat["y_init"]-cat["y_fit"])**2.0)
        cat.add_column(Column(d,name="xydev"))

        if "flux_err" not in cat.colnames:
            cat.add_column(Column(np.full(len(cat),np.nan), name="eflux"))
            perror("NO ERRORS??\n")
        else: cat.rename_column("flux_err","eflux")
        
        keep=["x_fit","y_fit","flux_fit","eflux","xydev","qfit"]
        #cat=hstack(( init_params, cat[keep]))
        return hstack((init_params, cat[keep]))

class ArtificialStar_Routine(object):
    """
    docstring
    """
    def __init__(self, detector, psffitter, psf):
        """
        detector - Detection class that fits the StarFinder base class
        psffitter - PSF fitting class that fits the IterativelySubtractedPSFPhotometry base class
        psf - DiscretePRF
        """
        self.detector=detector
        self.psffitter=psffitter
        self.psf=psf

        print("WARNING: THIS IS UNDER DEVELOPMENT")


    def run(self, image, ntests=1000, subimage_size=500,  sources=None, fwhm=1,
            flux_range=(0,1e5), separation_thresh=2, save_progress=1):
        """
        run artificial star testing on an image
        parameters:
             ntests - number of tests to conduct
             subimage_size - size of the cropped subimage 
             sources - precalculated positions to test stars
                astropy table with x_0 y_0 flux columns
            fwhm - FWHM of the stars to be added,
                this is used to ensure the source isnt too close to the border
            flux_range - range of fluxes to test
            separation_thresh - number pixels above which the separation is too high and the artificial star failed
            save_progress - periodically save the catalogue 



            for each star, plaxe it into a copy of the base image.
            however, the base image should be cropped in and around that new star
        """
        #printf("starting artificial star testing..\n")
        shape=np.array(image.shape)
        psfsize=self.psf.shape
        if np.any( subimage_size>shape):
            warn()
            perror("subimage_size bigger than image dimensions\n")
            subimage_size=min(shape)
        subimage_size=int(subimage_size)

        if not sources:
            x_range=[ 2.0*fwhm, shape[0]-(2.0*fwhm)]
            y_range=[ 2.0*fwhm, shape[1]-(2.0*fwhm)]
            sources=make_random_models_table(int(ntests), {"x_0":x_range, "y_0":y_range, "flux":flux_range}, seed=int(time.time()))
        
        sources.add_column(Column(np.zeros(len(sources)), name="outflux"))
        sources.add_column(Column(np.zeros(len(sources)), name="x_det"))
        sources.add_column(Column(np.zeros(len(sources)), name="y_det"))
        sources.add_column(Column(np.zeros(len(sources)), name="status"))

        load=loading(len(sources), msg="artificial star tests")
        load.show()
        for n,src in enumerate(sources):

            subx=0
            suby=0
            if subimage_size>0:
                ## !! I might change this to be PSFSIZE not 2FWHM
                subx = np.random.randint( max(0, src['x_0']+(2*fwhm)-subimage_size), min(shape[0]-subimage_size, src['x_0']-(2*fwhm)))
                suby = np.random.randint( max(0, src['y_0']+(2*fwhm)-subimage_size), min(shape[1]-subimage_size, src['y_0']-(2*fwhm)))
                #subx = np.random.randint( max(0, src['x_0']+(psfsize[0]/2)-subimage_size), min(shape[0]-subimage_size, src['x_0']-(psfsize[0]/2)))
                #suby = np.random.randint( max(0, src['y_0']+(psfsize[1]/2)-subimage_size), min(shape[1]-subimage_size, src['y_0']-(psfsize[1]/2)))

            src_mod=Table(src)# src mod translates the position within the subimage
            src_mod["x_0"]-=subx
            src_mod["y_0"]-=suby
            sky=image[subx:subx+subimage_size,suby:suby+subimage_size]
            base=np.copy(sky)+ make_model_sources_image(2*[subimage_size], self.psf, src_mod)
            #base=np.copy(image)+make_model_sources_image(shape, self.psf, Table(src))

            detections=self.detector(base)
            detections.rename_column("xcentroid", "x_0")
            detections.rename_column("ycentroid", "y_0")

            separations=(src_mod["x_0"]-detections["x_0"])**2+(src_mod["y_0"]-detections["y_0"])**2
            best_match=np.argmin(separations)
            if np.sqrt(separations[best_match]) <= separation_thresh:
                psftab= self.psffitter(base, init_guesses=detections)
                index=np.where( psftab["id"]==detections[best_match]["id"])

                sources[n]["outflux"]=psftab[index]["flux_fit"]
                sources[n]["x_det"]=psftab[index]["x_0"]+subx
                sources[n]["y_det"]=psftab[index]["y_0"]+suby

                if abs(sources[n]["outflux"]-sources[n]["flux"]) < (sources[n]["flux"]/100.0):
                    sources[n]["status"]=1 # star matched
            load()
            load.show()

            if save_progress and not n%10:
                export_table(sources[0:n], fname="/tmp/artificial_stars.save")

        return sources


class SourceProperties:
    status=0
    def __init__(self, image, sourcelist, filter=None, verbose=1):
        self.image=image
        self.sourcelist=None
        self.verbose=verbose
    
        if sourcelist and type(sourcelist) in (Table,QTable):
            if len( set(("xcentroid","ycentroid")) & set(sourcelist.colnames))==2:
                self.sourcelist=Table(sourcelist[["xcentroid","ycentroid"]])
            elif len( set(("x_0","y_0")) & set(sourcelist.colnames))==2:
                self.sourcelist=Table(sourcelist[["x_0","y_0"]])
                self.sourcelist.rename_columns( ("x_0","y_0"), ("xcentroid","ycentroid"))
            else: perror("no posisional columns in sourcelist\n")
        else: perror("bad sourcelist type: %s\n"%type(sourcelist))


    def __call__(self, do_crowd=1, **kwargs):
        """

        """
        out=Table()#self.sourcelist.copy()
        if do_crowd: ## This can be slow
            out=hstack((out,Table([self.calculate_crowding(**kwargs)],names=["crowding"])))

        out=hstack((out,self.calculate_geometry(**kwargs)))
        return out

    def calculate_crowding(self,N=10, **kwargs):
        """
        Crowding Index: Sum of magnitude of separation of N closest sources
        """
        if self.sourcelist is None: 
            perror("no sourcelist\n")
            return None

        crowd=np.zeros(len(self.sourcelist))
        load=loading(len(self.sourcelist),msg="calculating crowding", res=10)

        for i,src in enumerate(self.sourcelist):
            dist=np.sqrt( (src["xcentroid"]-self.sourcelist["xcentroid"])**2 + (src["ycentroid"]-self.sourcelist["ycentroid"])**2 )
            dist.sort()
            crowd[i]= sum( dist[1:N])
            load()
            if self.verbose: load.show()
        return crowd

    def calculate_geometry(self, fwhm=2, **kwargs):
        """

        """
        if self.sourcelist is None:
            perror("no sourcelist\n")
            return None
        if self.verbose: printf("-> measuring source geometry\n")
        xycoords=np.array((self.sourcelist["xcentroid"], self.sourcelist["ycentroid"])).T

        daofind=DAOStarFinder(-np.inf, fwhm, sharplo=-np.inf, sharphi=np.inf, roundlo=-np.inf, roundhi=np.inf, xycoords=xycoords, peakmax=np.inf)
        return daofind._get_raw_catalog(self.image).to_table()






if __name__=="__main__":
    tab=Table.read("/dat/ngc346/jwst/stage2/F444W/jw01227002001_02105_00001_nrcalong_cal_destrip-ap.fits")
    crowding=SourceProperties(None, tab).calculate_crowding()
    ii=np.argsort(crowding)[len(crowding)//2:]


    import matplotlib.pyplot as plt
    plt.scatter( tab["RA"],tab["DEC"], c='k')
    plt.scatter( tab["RA"][ii],tab["DEC"][ii], c='r')
    plt.show()

    #sb=StarbugBase("/home/conor/dat/NGC346/MIRAGE/Pipeline_Level3/ngc346-f115w-mosaic_i2d-cropped.fits")
    #sb.detect()
    #sb.photometry()
    #sb.export()
    #sb.artificial_stars()
    #sb.detect()
    #sb.photometry()
    #sb.export()




       
