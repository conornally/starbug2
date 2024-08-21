import numpy as np
from photutils.datasets import make_model_sources_image, make_random_models_table
from photutils.psf import FittableImageModel
from astropy.table import Table,hstack,vstack
from scipy.optimize import curve_fit

from starbug2.utils import perror, cropHDU, get_MJysr2Jy_scalefactor, warn
from starbug2.matching import GenericMatch

class Artificial_StarsIII():
    """
    ast
    """
    def __init__(self, starbug, index=-1):
        ## Initialis the starbug instance
        self.starbug=starbug
        _=self.starbug.image
        _=self.starbug.load_psf()

        self.psf=FittableImageModel(self.starbug.psf)
        self.index=index

    def __call__(self,*args,**kwargs): return self.auto_run(*args,**kwargs)

    def auto_run(self, ntests, stars_per_test=1, subimage_size=-1, mag_range=(18,27),
            loading_buffer=None, autosave=-1):
        """
        The main entry point into the artificial star test
        This handles everything except the results compilation at the end,

        Parameters
        ----------
        ntests : int
            Number of tests to run

        stars_per_test : int
            Number of stars to inject per test

        subimage_size : int
            in prep.

        mag_range : tuple,list
            Length two list or tuple containing the magnitude range of
            injected stars. These will be uniformly sampled from within 
            this range.

        loading_buffer : numpy.ndarray
            Length 3 array of shared memory to increment a loading bar
            between multiple subprocesses..

        autosave : int
            Auto quick saving output frequency
            
        Returns
        -------
        test_result : astropy.Table
            Full raw test results. Injected initial properties with measured values
        """

        test_result=Table(np.full((ntests*stars_per_test,8),np.nan), names=["x_0","y_0","mag","flux","x_det","y_det","flux_det", "status"])
        scalefactor= get_MJysr2Jy_scalefactor(self.starbug.image)
        base_image=self.starbug._image.copy()
        base_shape=np.copy(self.starbug.image.shape)
        stars_per_test=int(stars_per_test)
        passed=0

        ZP = self.starbug.options.get("ZP_MAG") if self.starbug.options.get("ZP_MAG") else 0
        buffer=0

        if mag_range[0]-mag_range[1] >=0:
            warn("Detected magnitude range in wrong order, put bright limit first\n")
            return None

        if any(base_shape < subimage_size):
            subimage_size=min(base_shape)
            perror("subimage size greater than image size, setting to 'safe' value %d.\n"%subimage_size)
            
        for test in range(1,int(ntests)+1):
            centre=0
            centre= (base_shape[0]*np.random.random(), base_shape[1]*np.random.random())

            #image=cropHDU( base_image.__deepcopy__(), (0,-1), (0,-1) )
            image=base_image.__deepcopy__()
            #image=self.create_subimage( base_image.__deepcopy__(), subimage_size, position=centre, hdu=self.st

            shape=image[self.starbug._nHDU].shape

            sourcelist= make_random_models_table( stars_per_test, { "x_0":[buffer,shape[0]-buffer],
                                                                    "y_0":[buffer,shape[1]-buffer],
                                                                    "mag":mag_range}) 
            sourcelist.add_column( 10.0 ** ( (ZP-sourcelist["mag"])/2.5 ) , name="flux")
            sourcelist.remove_column("id")

            #image[self.starbug._nHDU].data*=0
            star_overlay=make_model_sources_image( shape, self.psf, sourcelist)/scalefactor
            image[self.starbug._nHDU].data+=star_overlay
            self.starbug._image=image
            
            n=len(sourcelist)
            result=self.single_test(image, sourcelist)
            passed+=sum(result["status"])
            test_result[(test-1)*stars_per_test: test*stars_per_test]=result

            if loading_buffer is not None:
                loading_buffer[0]+=1
                loading_buffer[2]=int(100*passed/(test*stars_per_test))

            if autosave>0 and not test%autosave:
                test_result.write("sbast-autosave%d.tmp"%self.index, overwrite=True, format="fits")
            del image # is this neccessary?
        return test_result

    def single_test(self, image, contains):
        """
        Conduct a single test on an image with a set of initial source properties

        Parameters
        ----------
        image : numpy.ndarray
            2D image array to conduct test on

        contains : table
            Table of initial source properties to be injected into the image.
            This table must contain the columns ("x_0","y_0","flux")

        Returns
        -------
        result : Table
            Table hoizontally stacked with the initial inputs and the detection and 
            photometric results. Plus column named "status", an integer flag as to 
            whether the source was detected or not.
        """
        NULL=0
        DETECT=1
        test_result=Table(np.full((len(contains),4),np.nan), names=["x_det","y_det","flux_det","status"])

        threshold=2
        if not self.starbug.detect(): #Run detection on the image
            det=self.starbug.detections
            for i, src in enumerate(contains): #Check for detection in output
                separations=np.sqrt( (src["x_0"]-det["xcentroid"])**2 + (src["y_0"]-det["ycentroid"])**2)
                best_match=np.argmin(separations)
                if separations[best_match]<threshold:
                    test_result["x_det"][i]=det["xcentroid"][best_match]
                    test_result["y_det"][i]=det["ycentroid"][best_match]
                    test_result["flux_det"][i]=det["flux"][best_match]
                    test_result["status"][i]=DETECT
                else: test_result["status"][i]=NULL

            if sum(test_result["status"]) and not self.starbug.bgd_estimate():  # Run background estim if there were detections
                self.starbug.detections = test_result

                if not self.starbug.photometry(): # Run PSF photometry on detected sources
                    self.starbug.psfcatalogue.rename_columns(("x_init","y_init","xydev"),("_x_init","_y_init","_xydev"))
                    matched=GenericMatch(threshold=threshold)([contains, self.starbug.psfcatalogue], cartesian=True)
                    test_result["flux_det"] = matched[:len(test_result)]["flux_2"]

        return hstack((contains,test_result))


    def create_subimage(self, image, size, position=(0,0), hdu=1, buffer=0):
        """
        probably to be deprecated
        """
        subimage=None
        imshape=00
        if size<=0: return image,0,0
        x_edge=0
        y_edge=0
        if any(imshape < size):
            size=min(imshape)
            perror("subimage size greater than image size, setting to 'safe' value %d.\n"%size)

        x_edge = int(max( position[0]-(size/2), buffer ))
        y_edge = int(max( position[1]-(size/2), buffer ))
        x_end =  int(min( position[0]+(size/2), imshape[0]-buffer))
        y_end =  int(min( position[1]+(size/2), imshape[1]-buffer))

        return cropHDU(image,xlim=(x_edge,x_end),ylim=(y_edge,y_end)), x_edge, y_edge

    @staticmethod
    def get_completeness(test_result):
        """
        Compile the results into magnitude binned values of recovery fraction
        and flux error

        Parameters
        ----------
        test_result : table
            The output from auto_run

        Returns
        -------
        result : astropy Table
            Table containing percent completeness as a function of magnitude
        """

        bins = np.arange( np.floor(min(test_result["mag"])), np.ceil(max(test_result["mag"])), 0.1)
        percs= np.zeros(len(bins))
        errors=np.zeros(len(bins))
        offsets=np.zeros(len(bins))
        means =np.zeros(len(bins))
        
        ibins = np.digitize( test_result["mag"], bins=bins)
        for i in range(max(ibins)):
            binned=test_result[ (ibins==i) ]
            if binned: percs[i]=float(sum(binned["status"]))/len(binned)

            mag_inj= -2.5*np.log10( binned["flux"])
            mag_det= -2.5*np.log10( binned["flux_det"])
            errors[i]=np.nanstd( mag_inj-mag_det )
            means[i]=np.nanmean( mag_inj-mag_det )
            offsets[i]=np.nanmedian(binned["flux"]/binned["flux_det"])


        out=Table( [bins,percs,errors,offsets], names=("mag","rec","err","off"), dtype=(float,float,float,float))
        return out

    @staticmethod
    def get_spatialcompleteness(test_result,image,res=10):
        """
        Produce an image array showing the spatially dependant recovery fraction 

        Parameters
        ----------
        test_result : table
            The output from auto_run

        image : numpy.ndarry
            2D image array to take the shape from 

        res : int
            The resolution of the spatial bins

        Returns
        -------
        percs : numpy.ndarray
            A 2D array the same shape as the image input, pixel values
            show the fraction of injected sources recovered in this bin
        """
        xbins=np.arange(min(test_result["x_0"]),max(test_result["x_0"]), int(res))
        ybins=np.arange(min(test_result["y_0"]),max(test_result["y_0"]), int(res))
        percs=np.zeros(image.shape)
        for xi in xbins[:-1]:
            xo=xi+res
            for yi in ybins[:-1]:
                yo=yi+res
                mask=(test_result["x_0"]>=xi) & (test_result["x_0"]<xo) &(test_result["y_0"]>=yi) & (test_result["y_0"]<yo) 
                binned=test_result[mask] 
                if len(binned): percs[int(xi):int(xo),int(yi):int(yo)]=float(sum(binned["status"])/len(binned))
        return percs

    @staticmethod
    def estim_completeness(raw):
        """
        Estimate the completenss level of the artificial star test
        
        Parameters
        ----------
        raw : astropy Table
            Output of Artificial_Stars.get_completeness, table must contain columns (mag, rec)

        Returns
        -------
        fit : list
            The fitting parameters to the logistic curve f(x)=l/(1+exp(-k(x-xo)))
            fit=[l,xo,k]

        complete : list
            Magnitude of 70% and 50% completeness 
        """
        fit=[None,None,None]
        compl=[None,None,None]
        fn_i=lambda y,l,k,xo: xo-(np.log((l/y)-1)/k)

        if len(set(raw.colnames) & set(("mag","rec")))==2:
            fit,_=curve_fit(Artificial_StarsIII.scurve, raw["mag"], raw["rec"], [1, -1,np.median(raw["mag"])])
            compl=(fn_i(0.9,*fit),fn_i(0.7,*fit),fn_i(0.5,*fit))
        else: perror("Input table must have columns 'mag' and 'rec'\n")
        return fit,compl

    @staticmethod
    def scurve(x,l,k,xo):
        """
        S-curve function to fit completeness results to

        f(x)=l/(1+exp(-k(x-xo)))

        Parameters
        ----------
        x : list
            Magnitude range to input into function

        l,xo,k : float
            Function parameters

        Returns
        -------
        f(x) : float
        """
        return l/(1+np.exp(-k*(x-xo)))


#if __name__=="__main__":
#    from astropy.table import Table
#    import matplotlib.pyplot as plt
#    a=Table(fits.open("/home/conor/sci/docs/starbug/ws/F444W-ast.fits")[1].data)
#
#    fit,compl=Artificial_StarsIII.estim_completeness(a)
#    print(fit,compl)









#class Artificial_Stars(object):
#
#    def __init__(self, psf=None, detector=None, photometry=None):
#        self.psf=psf
#        self.detector=detector
#        self.photometry=photometry
#
#
#    def __call__(self, *args, **kwargs): return self.run_auto(*args,**kwargs)
#
#    def run_auto(self, data, background=None, ntests=100, stars_per_test=1,
#            subimage_size=100, buffer=0, flux_range=(1,1e6)):
#        """
#        """
#        load=loading(ntests, msg="artificial star testing")
#        test_result=None
#
#        ntests=int(ntests)
#        passed=0
#        stars_per_test=int(stars_per_test)
#
#        test_result=Table(None, names=["x_0","y_0","flux","x_det","y_det","flux_det", "status"])
#        for test in range(1,ntests+1):
#            x_edge=0
#            y_edge=0
#                
#            centre= (data.shape[0]*np.random.random(), data.shape[1]*np.random.random())
#            image,x_edge,y_edge=self.create_subimage(data,subimage_size, position=centre)
#            sourcelist= make_random_models_table( stars_per_test, { "x_0":[buffer,image.shape[0]-buffer],
#                                                                    "y_0":[buffer,image.shape[1]-buffer],
#                                                                    "flux":flux_range})
#            star_overlay=make_model_sources_image( image.shape, self.psf, sourcelist)    
#            
#            if background is not None:
#                # This feels like a bodge
#                self.photometry.background,_,_=self.create_subimage(background, subimage_size, position=centre)
#
#            if image.shape==star_overlay.shape:
#                result=self.single_test( image+star_overlay, sourcelist, threshold=2)
#                result["x_0"]+=x_edge
#                result["y_0"]+=y_edge
#                result["x_det"]+=x_edge
#                result["y_det"]+=y_edge
#                passed+=result["status"]
#                test_result=vstack((test_result,result))
#
#            load.msg="recovering %d%%"%(100*passed/test)
#            load()
#            load.show()
#
#        return test_result
#
#    def single_test(self, image, contains, background=None, threshold=2):
#        """
#        One single artifical star test. 
#
#        This will detect on the supplied image and check if the
#        input stars are recovered. Photometry will be conducted on 
#        recovered sources, if `self.photomotry` is not None
#
#        Parameters
#        ----------
#        image : 2d numpy array
#            The image (or subimage) onto which to place the source
#
#        contains : `astropy.table.Table`
#            A list of sources that have been added to the image.
#        
#        background : 
#
#        threshold : 
#
#        Results
#        -------
#        test_result : `astropy.table.Table`
#            A list of the results.
#        """
#        DETECT=1
#        NULL=0
#        test_result=Table(np.full((len(contains),4),np.nan), names=["x_det","y_det","flux_det", "status"])
#        detections=self.detector(image)
#        detections.rename_columns(("xcentroid","ycentroid"),("x_0","y_0"))
#
#        if detections:
#            for i,src in enumerate(contains):
#                separations=np.sqrt( (src["x_0"]-detections["x_0"])**2 + (src["y_0"]-detections["y_0"])**2)
#                best_match=np.argmin(separations)
#                if separations[best_match]<threshold:
#                    test_result["x_det"][i]=detections["x_0"][best_match]
#                    test_result["y_det"][i]=detections["y_0"][best_match]
#                    test_result["status"][i]=DETECT
#
#                    if self.photometry:
#                        # I think this could all be done simultaneously
#                        test_result["flux_det"][i]=self.photometry(image, Table(detections[["x_0","y_0"]][best_match]))["flux"]
#                else:
#                    test_result["status"][i]=NULL
#
#        return hstack((contains,test_result))
#
#    def create_subimage(self, image, size, position=(0,0), method="centre", buffer=0):
#        """
#        Create a subimage from a larger image
#
#        Parameters
#        ----------
#
#        image : 2d np.array
#            base image to cut out of
#
#        size : int, tuple
#            side length of subimage. The result will be a square unless size
#            is a tuple (size x, size y)
#
#        poisition : tuple
#            point which must be contained within the subimage
#
#        method : string
#            method to create the image around the point:
#            -   centre : point will be centred in the image (if possible)
#            -   random : point will exist somewhere within the image
#
#        buffer : int
#            Buffer around the edge which the including point cannot enter
#        """
#        subimage=None
#        imshape=np.array(image.shape)
#        x_edge=0
#        y_edge=0
#
#
#        if size<=0: return image,0,0
#
#        if any(imshape < size):
#            size=min(imshape)
#            perror("subimage size greater than image size, setting to 'safe' value %d.\n"%size)
#
#        if buffer <0 or buffer>size/2:
#            buffer=0
#            perror("buffer must be >=0 and < size/2, setting to 'safe' value zero.\n")
#
#        if False: pass ## position check
#
#        if method=="centre":
#            x_edge = int(max( position[0]-(size/2), buffer ))
#            y_edge = int(max( position[1]-(size/2), buffer ))
#            x_end =  int(min( position[0]+(size/2), imshape[0]-buffer))
#            y_end =  int(min( position[1]+(size/2), imshape[1]-buffer))
#
#
#        elif method=="random":
#            #|----------------|
#            #|      |=======| |
#            #|      |       | |
#            #|      |  x    | |
#            #|      |       | |
#            #|      |=======| |
#            #|                |
#            #|                |
#            #|----------------|
#            perror("not impleneted\n")
#            raise NotImplementedError 
#
#
#        subimage=image[ x_edge:x_end,y_edge:y_end]
#        return subimage, x_edge, y_edge
#        
#    def periodic_export(self,fname=""):
#        return 0
#
#
#def fnxep(x,a,b,c):
#    return a*np.exp(b*x+c)
