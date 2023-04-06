"""
Routines for Starbug
"""
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
from astropy.table import Column, Table, QTable, hstack
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.convolution import RickerWavelet2DKernel

from photutils.background import Background2D, BackgroundBase
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.detection import StarFinderBase, DAOStarFinder
from photutils.psf.groupstars import DAOGroup
from photutils.psf import BasicPSFPhotometry

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
                        sharplo=0.2, sharphi=1, roundlo=-1, roundhi=1,
                        verbose=0, cleansrc=1, bgd2d=1, boxsize=2):
        self.sig_src=sig_src
        self.sig_sky=sig_sky
        self.fwhm = fwhm
        self.sharphi=sharphi
        self.sharplo=sharplo
        self.roundhi=roundhi
        self.roundlo=roundlo
        self.cleansrc=cleansrc

        #self.match_threshold=u.Quantity(match_threshold)*u.dimensionless_unscaled
        self.catalogue=Table()
        self.verbose=verbose

        self.bgd2d=bgd2d
        self.boxsize=boxsize
        #self.filtersize=filtersize

    def detect(self, data, bkg_estimator=None, xycoords=None):
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

        _,_,std=sigma_clipped_stats(data,sigma=self.sig_sky)
        daofind=DAOStarFinder(std*self.sig_src, self.fwhm, sharplo=self.sharplo, sharphi=self.sharphi,
                roundlo=self.roundlo, roundhi=self.roundhi, peakmax=np.inf, xycoords=xycoords)
        return daofind(data - bkg)

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
        for src,sep in zip(cat,separation.to_value()):
            if sep>self.fwhm:#1 pixel?
                base.add_row(src)
                added+=1
        return added

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

        self.match(self.catalogue, self.detect(data, self._bkg2d))
        if self.verbose: printf("-> [BGD2D] pass: %d sources\n"%len(self.catalogue))

        ## 2nd order differential detection
        kernel=RickerWavelet2DKernel(1)
        conv=convolve(data, kernel)
        corr=match_template(conv/np.amax(conv), kernel.array)
        _detections=self.detect(corr)
        _detections["xcentroid"]+=kernel.shape[0]//2
        _detections["ycentroid"]+=kernel.shape[0]//2
        self.match(self.catalogue, _detections)
        if self.verbose: printf("-> [CONVL] pass: %d sources\n"%len(self.catalogue))

        ## Now with xycoords DAOStarfinder will refit the sharp and round values at the detected locations
        #self.catalogue=self.detect(data, xycoords=np.array([self.catalogue["xcentroid"],self.catalogue["ycentroid"]]).T)#, clean=0)
        tmp=SourceProperties(data,self.catalogue, verbose=self.verbose).calculate_geometry(self.fwhm)
        if tmp: self.catalogue=tmp

        mask=(~np.isnan(self.catalogue["xcentroid"]) & ~np.isnan(self.catalogue["ycentroid"]))
        self.catalogue.remove_rows(~mask)

        mask = ((self.catalogue["sharpness"]>self.sharplo)
               &(self.catalogue["sharpness"]<self.sharphi)
               &(self.catalogue["roundness1"]>self.roundlo)
               &(self.catalogue["roundness1"]<self.roundhi)
               &(self.catalogue["roundness2"]>self.roundlo)
               &(self.catalogue["roundness2"]<self.roundhi))
        if self.cleansrc: 
            if self.verbose: printf("-> cleaning %d+ unlikley point sources\n"%sum(~mask))
            self.catalogue.remove_rows(~mask)
        
        if self.verbose: printf("Total: %d sources\n"%len(self.catalogue))

        self.catalogue.replace_column("id", Column(range(1,1+len(self.catalogue))))

        return self.catalogue

class APPhot_Routine():
    """
    Aperture photometry called by starbug
    Given photometry radius, sky annuli radii rad_inner rad_outer
    """
    def __init__(self, radius, sky_in, sky_out, encircled_energy=50, fit_radius=1, verbose=0):
        self.radius=radius
        self.sky_in=sky_in
        self.sky_out=sky_out
        self.encircled_energy=encircled_energy
        self.fit_radius=fit_radius
        self.catalogue=Table(None)#, names=["ap_flux_r%d"%n for n in range(len(radii))]+["sky_median"])

        self.verbose=verbose


    def __call__(self, detections, image, **kwargs):
        return self.run(detections, image, **kwargs)

    #def run(self, detections, image, apcorr=1.0, sig_sky=3):
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
        annulus_aperture=CircularAnnulus(pos, r_in=self.sky_in, r_out=self.sky_out)
        phot=aperture_photometry(image, apertures, error=error, mask=mask)
        
        self.catalogue=Table(np.full((len(pos),3),np.nan),names=("flux","eflux","sky"))

        self.log("-> calculating photometric errors\n")
        for i,mask in enumerate(annulus_aperture.to_mask(method="center")):
            dat=np.array(mask.multiply(image),dtype=float)
            dat=sigma_clip(dat[dat>0 & np.isfinite(dat)], sigma=sig_sky)
            if len(dat): ##sometimes all the surrounding pixels are nan OR above SIGSKY value??
                #self.catalogue["sky"][i]=mode(dat)[0]
                self.catalogue["sky"][i]=np.ma.median(dat)


        self.catalogue["eflux"]=phot["aperture_sum_err"]
        self.catalogue["flux"]=apcorr*(phot["aperture_sum"] - (self.catalogue["sky"]*apertures.area))

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
        if not table_fname: return None
        tmp=Table.read(table_fname, format="fits")
        t_apcorr=tmp[tmp["filter"]==filter][:-1] ## SKIPPING FINAL VALUE

        fn=lambda x,a,b,c: a*np.exp(-b*x)+c
        popt,_=curve_fit(fn, t_apcorr["radius"], t_apcorr["apcorr"])
        apcorr= fn( radius,*popt)
        if verbose: printf("-> estimating aperture correction: %.3g..\n"%apcorr)

        popt,_=curve_fit(fn,t_apcorr["radius"],t_apcorr["eefraction"])
        eefrac= fn(radius,*popt)
        if verbose: printf("-> effective encircled energy fraction: %.3g\n"%eefrac)
        return apcorr


    @staticmethod
    def apcorr_from_encenergy(filter, encircled_energy, table_fname=None, verbose=0):
        """
        Rather than fitting radius to the APCORR CRDS, use the closes Encircled energy value
        """
        if not table_fname: return None
        tmp=Table.read(table_fname, format="fits")
        t_apcorr=tmp[tmp["filter"]==filter]

        line=t_apcorr[(np.abs(t_apcorr["eefraction"]-encircled_energy)).argmin()]
        if verbose:
            printf("-> best matching encircled energy %.1f, with radius %g pixels\n"%(line["eefraction"],line["radius"]))
            printf("-> using aperture correction: %f\n"%line["apcorr"])

        return line["apcorr"], line["radius"]

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
    def __init__(self, sourcelist, boxsize=2, fwhm=2, verbose=0, bgd=None):#mask_r0=7, mask_r1=9
        self.sourcelist=sourcelist
        self.boxsize=boxsize
        self.fwhm=fwhm
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

    def __call__(self, data, axis=None, masked=False):
        if self.sourcelist is None or data is None: return self.bgd
        _data=np.copy(data)
        X,Y=np.ogrid[:data.shape[1], :data.shape[0]]
        peaks=self.calc_peaks(data)

        rlist=np.sqrt(peaks**0.7)*self.fwhm/1.5 ## <-- that works but hmm
        D=50
        load=loading(len(self.sourcelist), msg="masking sources", res=10)
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
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ngroups=0
    def __call__(self, *args):
        res=super().__call__(*args)
        self.ngroups=max(res["group_id"])
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


class PSFPhot_Routine(BasicPSFPhotometry):
    """
    PSF Photometry routine called by starbug
    """
    def __init__(self, crit_separation, psf_model, fitshape,
            force_fit=False, dposition_threshold=2, background=None, verbose=1):

        self.verbose=verbose
        self.force_fit=force_fit        
        self.dpos_thresh=dposition_threshold

        group_maker=_grouping(crit_separation=crit_separation)
        bkg_estimator=BackGround_Estimate_Routine(None, bgd=background)
        fitter=_fitmodel(grouper=group_maker, verbose=verbose)
        
        if force_fit:
            psf_model.fixed.update( {"x_0":True,"y_0":True} )
            if fitter.load is not None: fitter.load.msg+=" (forced)"

        super().__init__(group_maker=group_maker, bkg_estimator=bkg_estimator,
                psf_model=psf_model, fitshape=fitshape,
                finder=None, fitter=fitter)

    #def _bkg(self, axis=None,masked=None):
        #return self.background

    def do_photometry(self, image, mask=None, init_guesses=None, progress_bar=False):
        """
        """ 
        _cat=None
        _fixcat=None

        if init_guesses is None:
            perror("Must include source list\n")
            return None

        if self.verbose: printf("-> fitting %d sources\n"%len(init_guesses))
        cat=super().do_photometry(image, mask=mask, init_guesses=init_guesses, progress_bar=False)

        if "flux_unc" not in cat.colnames:
            cat.add_column(Column(np.full(len(cat),np.nan), name="eflux"))
            perror("NO ERRORS??\n")
        else: cat.rename_column("flux_unc","eflux")

        return cat

class Cleaning_Routine(object):
    """
    docstring...
    """
    sharp_mu=0
    sharp_sig=0

    round1_mu=0
    round1_sig=0

    round2_mu=0
    round2_sig=0

    verbose=1

    def __init__(self, cat, verbose=1 ):
        self.cat=cat
        self.verbose=verbose
        self.fit()

    def log(self,msg):
        if self.verbose: printf(msg)

    def fit(self):
        if "sharpness" in self.cat.colnames:
            self.sharp_mu, self.sharp_sig   = norm.fit(self.cat["sharpness"])
            self.log("fit 'sharpness': mu=%.2g sig=%.2g\n"%(self.sharp_mu, self.sharp_sig))

        if "roundness1" in self.cat.colnames:
            self.round1_mu, self.round1_sig = norm.fit(self.cat["roundness1"])
            self.log("fit 'roundness1': mu=%.2g sig=%.2g\n"%(self.round1_mu, self.round1_sig))

        if "roundness2" in self.cat.colnames:
            self.round2_mu, self.round2_sig = norm.fit(self.cat["roundness2"])
            self.log("fit 'roundness2': mu=%.2g sig=%.2g\n"%(self.round2_mu, self.round2_sig))

    def _remove_rows(self,cat,indices):
        length=len(cat)
        cat.remove_rows(indices)
        if "id" in cat.colnames: cat.replace_column("id", Column(range(1,1+len(cat))))
        return length-len(cat)

    def run(self, mag_unc=1, sharp_sig_hi=2, sharp_sig_lo=2,
                                        round_sig_hi=2, round_sig_lo=2):
        """
        if "ap_flux" in self.cat.colnames:
            dr=self._remove_rows(self.cat, self.cat["ap_flux"]<0)
            self.log("-> cut %d sources with ap_phot flux<0\n"%dr)
        """

        if "mag_unc" in self.cat.colnames:
            dr=self._remove_rows(self.cat, self.cat["mag_unc"]>mag_unc)
            self.log("-> cut %d sources with mag_error>%g\n"%(dr,mag_unc))

        if "sharpness" in self.cat.colnames:
            dr =self._remove_rows(self.cat, self.cat["sharpness"]>(self.sharp_mu + (sharp_sig_hi*self.sharp_sig)))
            dr+=self._remove_rows(self.cat, self.cat["sharpness"]<(self.sharp_mu - (sharp_sig_lo*self.sharp_sig)))
            self.log("-> cut %d sources on sharpness\n"%dr)

        if "roundness1" in self.cat.colnames:
            dr =self._remove_rows(self.cat, self.cat["roundness1"]>(self.round1_mu + (round_sig_hi*self.round1_sig)))
            dr+=self._remove_rows(self.cat, self.cat["roundness1"]<(self.round1_mu - (round_sig_lo*self.round1_sig)))
            self.log("-> cut %d sources on roundness1\n"%dr)

        if "roundness2" in self.cat.colnames:
            dr =self._remove_rows(self.cat, self.cat["roundness2"]>(self.round2_mu + (round_sig_hi*self.round2_sig)))
            dr+=self._remove_rows(self.cat, self.cat["roundness2"]<(self.round2_mu - (round_sig_lo*self.round2_sig)))
            self.log("-> cut %d sources on roundness2\n"%dr)

        self.log("final catalogue size: %d\n"%len(self.cat))
        return self.cat




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




       
