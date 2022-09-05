"""
Routines for Starbug
"""
import matplotlib.pyplot as plt
import os
import sys
import time
import numpy as np
import pkg_resources
from scipy.stats import norm, mode
from scipy.optimize import curve_fit

from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats, sigma_clip
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Column, Table
import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter

from photutils.background import Background2D, SExtractorBackground, BackgroundBase
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.detection import StarFinderBase, DAOStarFinder
from photutils.psf.groupstars import DAOGroup
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.psf.sandbox import DiscretePRF
from photutils.datasets import make_model_sources_image, make_random_models_table
from starbug2.utils import loading, split_fname, tabppend, printf, perror, extnames
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
    progress=0
    def __init__(self,  sig_src=5, sig_sky=3, fwhm=1,
                        sharplo=0.2, sharphi=1, roundlo=-1, roundhi=1,
                        wcs=None, match_threshold=1.0, verbose=0,
                        boxsize=2, filtersize=3):
        self.sig_src=sig_src
        self.sig_sky=sig_sky
        self.fwhm = fwhm
        self.sharphi=sharphi
        self.sharplo=sharplo
        self.roundhi=roundhi
        self.roundlo=roundlo

        self.wcs=wcs
        self.match_threshold=match_threshold*u.arcsec
        self.catalogue=None
        self.verbose=verbose

        self.boxsize=boxsize
        self.filtersize=filtersize

    def detect(self, data, bkg_estimator=None, xycoords=None, clean=1):
        """
        docstring
        """
        bkg=np.zeros(data.shape)
        if bkg_estimator:
            bkg=bkg_estimator(data)
            """
            try:
                fits.PrimaryHDU(data=data-bkg).writeto("/tmp/bgd.fits",overwrite=True)
            except: pass
            """

        mean,median,std=sigma_clipped_stats(data,sigma=self.sig_sky)
        if clean:
            daofind=DAOStarFinder(std*self.sig_src, self.fwhm, sharplo=self.sharplo, sharphi=self.sharphi,
                    roundlo=self.roundlo, roundhi=self.roundhi, xycoords=xycoords)
        else: ## i dont want it to lose sharp/round values on the final run
            daofind=DAOStarFinder(std*self.sig_src, self.fwhm, sharplo=-999, sharphi=999,
                    roundlo=-999, roundhi=999, xycoords=xycoords)

        return daofind(data - bkg)

    def _bkg2d(self, data):
        return Background2D(data, self.boxsize, filter_size=self.filtersize).background

    def match(self, base, cat):
        """
        Internal function to class
        Used to match detenctoins from separate background subtracted images
        into the main catalogue. This will append a source if its matched separation
        is above the threshold self.match_threshold
        """
        added=0
        b_ra, b_dec = self.wcs.all_pix2world( base["xcentroid"], base["ycentroid"], 0)
        c_ra, c_dec = self.wcs.all_pix2world( cat["xcentroid"], cat["ycentroid"], 0)

        base_sky=SkyCoord(b_ra*u.degree, b_dec*u.degree)
        cat_sky=SkyCoord(c_ra*u.degree, c_dec*u.degree)

        _,separation,_=cat_sky.match_to_catalog_3d(base_sky)
        for src,sep in zip(cat,separation):
            if sep>self.match_threshold:
                base.add_row(src)
                added+=1
        return added

    def find_stars(self, data, mask=None):
        self.catalogue=self.detect(data)
        if self.verbose: printf("-> [PLAIN] pass: %d sources\n"%len(self.catalogue))
        sigma_clip = SigmaClip(sigma=self.sig_sky, maxiters=10)

        #self.match(self.catalogue, self.detect(data, MedianBackground(sigma_clip=sigma_clip)))

        self.match(self.catalogue, self.detect(data, SExtractorBackground(sigma_clip=sigma_clip)))
        if self.verbose: printf("-> [SExTr] pass: %d sources\n"%len(self.catalogue))

        self.match(self.catalogue, self.detect(data, self._bkg2d))
        if self.verbose: printf("-> [BGD2D] pass: %d sources\n"%len(self.catalogue))

        ## Now with xycoords DAOStarfinder will refit the sharp and round values at the detected locations
        self.catalogue=self.detect(data, xycoords=np.array([self.catalogue["xcentroid"],self.catalogue["ycentroid"]]).T)#, clean=0)

        if self.verbose: printf("Total: %d sources\n"%len(self.catalogue))

        self.catalogue.replace_column("id", Column(range(1,1+len(self.catalogue))))
        ra,dec=self.wcs.all_pix2world(self.catalogue["xcentroid"], self.catalogue["ycentroid"],0)
        self.catalogue.add_column( Column(ra, name="RA"), index=1)
        self.catalogue.add_column( Column(dec, name="DEC"), index=2)

        return self.catalogue

class APPhot_Routine():
    """
    Aperture photometry called by starbug
    Given photometry radius, sky annuli radii rad_inner rad_outer
    """
    def __init__(self, radius, sky_in, sky_out, verbose=0):
        self.radius=radius
        self.sky_in=sky_in
        self.sky_out=sky_out
        self.catalogue=Table(None)#, names=["ap_flux_r%d"%n for n in range(len(radii))]+["sky_median"])

        self.verbose=verbose


    def __call__(self, detections, image, **kwargs):
        return self.run(detections, image, **kwargs)

    def run(self, detections, image, apcorr=1.0, sig_sky=3):
        """
        Forced aperture photometry on a list of detections
        detections are a astropy.table.Table with columns xcentroid ycentroid or x_0 y_0
        This will add extra columns into this table ap_flux ap_sky
        """
        try:
            pos=[(line["xcentroid"],line["ycentroid"]) for line in detections]
        except:
            pos=[(line["x_0"],line["y_0"]) for line in detections]

        area=1.0 if image[4].name!="AREA" else image["AREA"].data
        mask= image["DQ"].data &( DQ_DO_NOT_USE| DQ_SATURATED | DQ_JUMP_DET)

        data= image["SCI"].data * area
        data[mask]=np.nan
        error=np.sqrt(data)

        apertures=CircularAperture(pos,self.radius)
        annulus_aperture=CircularAnnulus(pos, r_in=self.sky_in, r_out=self.sky_out)
        phot=aperture_photometry(data, apertures, error=error, mask=mask)
        
        self.catalogue=Table(np.full((len(pos),3),np.nan),names=("flux","eflux","sky"))
        bkg_stats=np.full((len(annulus_aperture),2), np.nan) # (mode, stdev) 

        self.log("-> calculating photometric errors\n")
        load=loading(len(apertures),"error?",res=10)
        for i,mask in enumerate(annulus_aperture.to_mask(method="center")):
            dat=mask.multiply(data)
            dat=sigma_clip(dat[dat>0 & np.isfinite(dat)], sigma=sig_sky)
            bkg_stats[i]=sigma_clipped_stats(dat)[1:]
            #bkg_stats[i]=( mode(dat)[0], np.nanstd(dat)) ##Take the mode to estimate background level


        err=np.full( (len(phot),3), np.nan)
        err[:,0]=phot["aperture_sum_err"].value      ##POISSON ERR
        err[:,1]=apertures.area * bkg_stats[:,1]**2  ##SKY_ERR
        err[:,2]=(apertures.area**2) * (bkg_stats[:,1]**2) / annulus_aperture.area ##MEANSKY_ERR
        self.catalogue["eflux"]=np.sqrt(np.nansum(err,axis=1))

        self.catalogue["sky"]=bkg_stats[:,0]
        self.catalogue["flux"]=apcorr*(phot["aperture_sum"] - (self.catalogue["sky"]*apertures.area))

        if "DQ" in extnames(image):
            self.log("-> flagging unlikely sources\n")
            col=Column(np.full(len(apertures),SRC_GOOD), dtype=np.uint32, name="flag")
            for i, mask in enumerate(apertures.to_mask(method="center")):
                dat=np.array(mask.multiply(image["DQ"].data),np.uint32)
                if np.sum( dat & (DQ_DO_NOT_USE|DQ_SATURATED)): col[i]|=SRC_BAD
                if np.sum( dat & DQ_JUMP_DET): col[i]|=SRC_JMP
            self.catalogue.add_column(col)

        return self.catalogue

    def calc_apcorr(self, filter, radius, table_fname=None):
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
        self.log("-> estimating aperture correction: %.3g..\n"%apcorr)

        popt,_=curve_fit(fn,t_apcorr["radius"],t_apcorr["eefraction"])
        eefrac= fn(radius,*popt)
        self.log("-> effective encircled energy fraction: %.3g\n"%eefrac)
        return apcorr

    def log(self,msg):
        if self.verbose: 
            printf(msg)
            sys.stdout.flush()


class _grouping(DAOGroup):
    """
    Overwritten DAOGroup that just holds the number of groups
    for use in verbose loading of psfphot routine
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
            self.load()
            self.load.show()
        return super().__call__(*args,**kwargs)

class BackGround_Estimate_Routine(BackgroundBase):
    bgd=None
    sourcelist=[]
    def __init__(self, sourcelist, boxsize=2, fwhm=2, verbose=0):#mask_r0=7, mask_r1=9
        self.sourcelist=sourcelist
        self.boxsize=boxsize
        self.fwhm=fwhm
        self.verbose=verbose

    def __call__(self, data):
        _data=np.copy(data)
        X,Y=np.ogrid[:data.shape[1], :data.shape[0]]
        rlist=np.sqrt(self.sourcelist["peak"]**0.7)*self.fwhm
        D=50
        load=loading(len(self.sourcelist), msg="masking sources")
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

            #fp.write("circle %f %f %f # color=red  \n"%( src["xcentroid"]+1, src["ycentroid"]+1, r))
            #fp.write("circle %f %f %f # color=green\n"%( src["xcentroid"]+1, src["ycentroid"]+1, rin))
            #fp.write("circle %f %f %f # color=green\n"%( src["xcentroid"]+1, src["ycentroid"]+1, rout))

            load()
            load.show() ## This will slow the thing down quite a lot
        self.bgd=Background2D(_data, self.boxsize).background
        return self.bgd

    def calc_background(self,data, axis=None, masked=None):
        if not self.bgd: self.__call__(data)
        return self.bgd


class PSFPhot_Routine(IterativelySubtractedPSFPhotometry):
    """
    PSF Photometry routine called by starbug
    It almost exactly matches DAOPhotRoutine however has some
    verbose progress bar outputs so you can check if the program is
    still running
    """
    def __init__(self, crit_separation, fwhm, psf_model, fitshape,
                 sig_sky=3, sig_src=5,
                 sharplo=0.2, sharphi=1.0, roundlo=-1.0, roundhi=1.0,
                 background=None, wcs=None, verbose=1):

        group_maker=_grouping(crit_separation=crit_separation)

        #self.background=background
        bkg_estimator=BackGround_Estimate_Routine(None)
        bkg_estimator.bgd=background

        finder=Detection_Routine(sig_src=sig_src, sig_sky=sig_sky, fwhm=fwhm,
                sharplo=sharplo, sharphi=sharphi, roundlo=roundlo, roundhi=roundhi,
                wcs=wcs)
        fitter=_fitmodel(grouper=group_maker, verbose=verbose)

        super().__init__(group_maker=group_maker, bkg_estimator=bkg_estimator,
                psf_model=psf_model, fitshape=fitshape,
                finder=finder, fitter=fitter, niters=1)

    def _bkg(self, axis=None,masked=None):
        return self.background

    def do_photometry(self, image, init_guesses=None):
        """
        if(init_guesses):
            init_guesses.rename_column("xcentroid", "x_0")
            init_guesses.rename_column("ycentroid", "y_0")
        """ 
        #if init_guesses: self.bkg_estimator.load_sources(init_guesses)
        cat=super().do_photometry(image, init_guesses) 

        cat.remove_columns(("flux_0", "group_id", "x_fit", "y_fit"))
        cat.remove_rows( cat["flux_fit"]<=0)
        cat.add_column(-2.5*np.log10(cat["flux_fit"]), name="mag_fit")
        cat.add_column((2.5/np.log(10))*(cat["flux_unc"]/cat["flux_fit"]), name="mag_unc")
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
            """
            import matplotlib.pyplot as plt
            plt.hist(self.cat["roundness1"],bins=40, density=True)
            x=np.linspace(min(self.cat["roundness1"]), max(self.cat["roundness1"]), 100)
            plt.plot(x, norm.pdf(x,self.round1_mu, self.round1_sig))
            plt.show()
            """


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
        psfsize=self.psf.prf_shape
        if np.any( subimage_size>shape):
            perror("warning: subimage_size bigger than image dimensions\n")
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

if __name__=="__main__":
    pass
    #sb=StarbugBase("/home/conor/dat/NGC346/MIRAGE/Pipeline_Level3/ngc346-f115w-mosaic_i2d-cropped.fits")
    #sb.detect()
    #sb.photometry()
    #sb.export()
    #sb.artificial_stars()
    #sb.detect()
    #sb.photometry()
    #sb.export()
    
