import numpy as np
from photutils.datasets import make_model_sources_image, make_random_models_table
from astropy.coordinates import SkyCoord
from astropy.table import Table,hstack,vstack
from starbug2.routines import Detection_Routine
from starbug2.utils import perror, loading, cropHDU, get_MJysr2Jy_scalefactor, warn
from starbug2.matching import GenericMatch

from astropy.io import fits

from photutils.psf import FittableImageModel
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit



class Artificial_StarsIII(object):
    def __init__(self, starbug):
        ## Initialis the starbug instance
        self.starbug=starbug
        _=self.starbug.image
        _=self.starbug.load_psf()

        self.psf=FittableImageModel(self.starbug.psf)

    def __call__(self,*args,**kwargs): return self.run_auto(*args,**kwargs)

    def auto_run(self, ntests, stars_per_test=1, subimage_size=-1, mag_range=(18,27),
            loading_buffer=None):
        """
        """

        test_result=Table(None, names=["x_0","y_0","flux","mag","x_det","y_det","flux_det", "status"])
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

            star_overlay=make_model_sources_image( shape, self.psf, sourcelist)/scalefactor
            image[self.starbug._nHDU].data+=star_overlay
            self.starbug._image=image
            
            n=len(sourcelist)
            result=self.single_test(image, sourcelist)
            passed+=sum(result["status"])
            test_result=vstack((test_result,result))

            if loading_buffer is not None:
                loading_buffer[0]+=1
                loading_buffer[2]=int(100*passed/(test*stars_per_test))
        return test_result

    def single_test(self, image, contains):
        """
        """
        NULL=0
        DETECT=1
        test_result=Table(np.full((len(contains),4),np.nan), names=["x_det","y_det","flux_det","status"])

        threshold=2
        if not self.starbug.detect():
            det=self.starbug.detections
            """
            match=GenericMatch(threshold=2, colnames=None )
            full= match([contains, self.starbug.detections], cartesian=True)[:len(contains)]
            status = ~np.isnan(full["flux_2"])

            test_result["x_det"] = full["xcentroid_2"]
            test_result["y_det"] = full["ycentroid_2"]
            test_result["flux_det"]= full["flux_2"]
            test_result["status"]=status.astype(int)

            self.starbug.detections=contains
            self.starbug.psfcatalogue=None
            self.starbug.photometry()
            """

            for i, src in enumerate(contains):
                separations=np.sqrt( (src["x_0"]-det["xcentroid"])**2 + (src["y_0"]-det["ycentroid"])**2)
                best_match=np.argmin(separations)
                if separations[best_match]<threshold:
                    test_result["x_det"][i]=det["xcentroid"][best_match]
                    test_result["y_det"][i]=det["ycentroid"][best_match]
                    test_result["flux_det"][i]=det["flux"][best_match]
                    test_result["status"][i]=DETECT

                    """
                    print(best_match, len(det))
                    self.starbug.detections = Table(det[best_match])
                    self.starbug.photometry()
                    test_result["flux_det"][i]=self.starbug.psfcatalogue["flux"]
                    """
                else: test_result["status"][i]=NULL

            #test_result.remove_rows( (np.isnan(test_result["x_det"])|np.isnan(test_result["y_det"])) )
            self.starbug.detections = test_result

            if sum(test_result["status"]):

                #if not self.starbug.bgd_estimate() and not self.starbug.photometry():
                if not self.starbug.photometry():
                    #self.starbug.psfcatalogue
                    #test_result["flux_det"][test_result["status"].value.astype(bool)]=self.starbug.psfcatalogue["flux"]

                    self.starbug.psfcatalogue.rename_columns(("x_init","y_init","xydev"),("_x_init","_y_init","_xydev"))
                    matched=GenericMatch(threshold=threshold)([contains, self.starbug.psfcatalogue], cartesian=True)
                    #print(self.starbug.psfcatalogue[["x_init","y_init","x_fit","y_fit"]])
                    #print(matched[["Catalogue_Number_2","x_0_1","y_0_1","x_fit_2","y_fit_2"]])
                    #print(len(test_result),len(matched))
                    test_result["flux_det"] = matched[:len(test_result)]["flux_2"]

        return hstack((contains,test_result))


    def create_subimage(self, image, size, position=(0,0), hdu=1, buffer=0):
        """
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

    def get_magerr(self, test_result):
        return None

    @staticmethod
    def get_completeness(test_result):
        """

        Returns:
        --------
        result : astropy Table
            Table containing percent completeness as a function of magnitude
        """

        bins = np.arange( np.floor(min(test_result["mag"])), np.ceil(max(test_result["mag"])), 0.1)
        percs= np.zeros(len(bins))
        errors=np.zeros(len(bins))
        means =np.zeros(len(bins))
        
        ibins = np.digitize( test_result["mag"], bins=bins)
        for i in range(max(ibins)):
            binned=test_result[ (ibins==i) ]
            if binned: percs[i]=float(sum(binned["status"]))/len(binned)

            mag_inj= -2.5*np.log10( binned["flux"])
            mag_det= -2.5*np.log10( binned["flux_det"])
            errors[i]=np.nanstd( mag_inj-mag_det )
            means[i]=np.nanmean( mag_inj-mag_det )


        #_s=curve_fit(fnxep,bins,percs)
        #a,b,c=_s[0]

        out=Table( [bins,percs,errors], names=("mag","rec","err"), dtype=(float,float,float))
        return out


        """
        mask= np.isnan(test_result["flux_det"])
        test_result=test_result[~mask]
        m,c=np.polyfit(test_result["flux_inj"], test_result["flux_det"],1)

        import matplotlib.pyplot as plt
        fig,(ax1,ax2)=plt.subplots(1,2)
        ax1.plot(bins,percs,c='k')
        x=np.linspace(min(bins),max(bins),100)
        ax1.plot(x,fnxep(x,a,b,c))
        #y=savgol_filter( percs, window_length=10, polyorder=2)
        #ax1.plot(bins,y, c='cyan')

        ax2.plot(bins,errors)

        plt.tight_layout()
        plt.show()
        """



class Artificial_Stars(object):

    def __init__(self, psf=None, detector=None, photometry=None):
        self.psf=psf
        self.detector=detector
        self.photometry=photometry


    def __call__(self, *args, **kwargs): return self.run_auto(*args,**kwargs)

    def run_auto(self, data, background=None, ntests=100, stars_per_test=1,
            subimage_size=100, buffer=0, flux_range=(1,1e6)):
        """
        """
        load=loading(ntests, msg="artificial star testing")
        test_result=None

        ntests=int(ntests)
        passed=0
        stars_per_test=int(stars_per_test)

        test_result=Table(None, names=["x_0","y_0","flux","x_det","y_det","flux_det", "status"])
        for test in range(1,ntests+1):
            x_edge=0
            y_edge=0
                
            centre= (data.shape[0]*np.random.random(), data.shape[1]*np.random.random())
            image,x_edge,y_edge=self.create_subimage(data,subimage_size, position=centre)
            sourcelist= make_random_models_table( stars_per_test, { "x_0":[buffer,image.shape[0]-buffer],
                                                                    "y_0":[buffer,image.shape[1]-buffer],
                                                                    "flux":flux_range})
            star_overlay=make_model_sources_image( image.shape, self.psf, sourcelist)    
            
            if background is not None:
                # This feels like a bodge
                self.photometry.background,_,_=self.create_subimage(background, subimage_size, position=centre)

            if image.shape==star_overlay.shape:
                result=self.single_test( image+star_overlay, sourcelist, threshold=2)
                result["x_0"]+=x_edge
                result["y_0"]+=y_edge
                result["x_det"]+=x_edge
                result["y_det"]+=y_edge
                passed+=result["status"]
                test_result=vstack((test_result,result))

            load.msg="recovering %d%%"%(100*passed/test)
            load()
            load.show()

        return test_result

    def single_test(self, image, contains, background=None, threshold=2):
        """
        One single artifical star test. 

        This will detect on the supplied image and check if the
        input stars are recovered. Photometry will be conducted on 
        recovered sources, if `self.photomotry` is not None

        Parameters
        ----------
        image : 2d numpy array
            The image (or subimage) onto which to place the source

        contains : `astropy.table.Table`
            A list of sources that have been added to the image.
        
        background : 

        threshold : 

        Results
        -------
        test_result : `astropy.table.Table`
            A list of the results.
        """
        DETECT=1
        NULL=0
        test_result=Table(np.full((len(contains),4),np.nan), names=["x_det","y_det","flux_det", "status"])
        detections=self.detector(image)
        detections.rename_columns(("xcentroid","ycentroid"),("x_0","y_0"))

        if detections:
            for i,src in enumerate(contains):
                separations=np.sqrt( (src["x_0"]-detections["x_0"])**2 + (src["y_0"]-detections["y_0"])**2)
                best_match=np.argmin(separations)
                if separations[best_match]<threshold:
                    test_result["x_det"][i]=detections["x_0"][best_match]
                    test_result["y_det"][i]=detections["y_0"][best_match]
                    test_result["status"][i]=DETECT

                    if self.photometry:
                        # I think this could all be done simultaneously
                        test_result["flux_det"][i]=self.photometry(image, Table(detections[["x_0","y_0"]][best_match]))["flux"]
                else:
                    test_result["status"][i]=NULL

        return hstack((contains,test_result))

    def create_subimage(self, image, size, position=(0,0), method="centre", buffer=0):
        """
        Create a subimage from a larger image

        Parameters
        ----------

        image : 2d np.array
            base image to cut out of

        size : int, tuple
            side length of subimage. The result will be a square unless size
            is a tuple (size x, size y)

        poisition : tuple
            point which must be contained within the subimage

        method : string
            method to create the image around the point:
            -   centre : point will be centred in the image (if possible)
            -   random : point will exist somewhere within the image

        buffer : int
            Buffer around the edge which the including point cannot enter
        """
        subimage=None
        imshape=np.array(image.shape)
        x_edge=0
        y_edge=0


        if size<=0: return image,0,0

        if any(imshape < size):
            size=min(imshape)
            perror("subimage size greater than image size, setting to 'safe' value %d.\n"%size)

        if buffer <0 or buffer>size/2:
            buffer=0
            perror("buffer must be >=0 and < size/2, setting to 'safe' value zero.\n")

        if False: pass ## position check

        if method=="centre":
            x_edge = int(max( position[0]-(size/2), buffer ))
            y_edge = int(max( position[1]-(size/2), buffer ))
            x_end =  int(min( position[0]+(size/2), imshape[0]-buffer))
            y_end =  int(min( position[1]+(size/2), imshape[1]-buffer))


        elif method=="random":
            #|----------------|
            #|      |=======| |
            #|      |       | |
            #|      |  x    | |
            #|      |       | |
            #|      |=======| |
            #|                |
            #|                |
            #|----------------|
            perror("not impleneted\n")
            raise NotImplementedError 


        subimage=image[ x_edge:x_end,y_edge:y_end]
        return subimage, x_edge, y_edge
        
    def periodic_export(self,fname=""):
        return 0


def fnxep(x,a,b,c):
    return a*np.exp(b*x+c)
