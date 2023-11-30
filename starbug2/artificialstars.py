import numpy as np
from photutils.datasets import make_model_sources_image, make_random_models_table
from astropy.coordinates import SkyCoord
from astropy.table import Table,hstack,vstack
from starbug2.routines import Detection_Routine
from starbug2.utils import perror, loading


class Artificial_Stars(object):

    def __init__(self, psf=None, detector=None, photometry=None):
        self.psf=psf
        self.detector=detector
        self.photometry=photometry


    def __call__(self, *args, **kwargs): return self.run_auto(*args,**kwargs)

    def run_auto(self, data, ntests=100, stars_per_test=1,
            subimage_size=100, buffer=0, 
            flux_range=(1,1e6)):
        """
        """
        load=loading(ntests, msg="artificial star testing")
        test_result=None

        if subimage_size>0:

            test_result=Table(None, names=["x_0","y_0","flux","x_det","y_det","flux_det", "status"])
            for test in range(ntests):
                x_edge=0
                y_edge=0
                    
                centre= (data.shape[0]*np.random.random(), data.shape[1]*np.random.random())
                image,x_edge,y_edge=self.create_subimage(data,subimage_size, position=centre)
                sourcelist= make_random_models_table( stars_per_test, { "x_0":[buffer,image.shape[0]-buffer],
                                                                        "y_0":[buffer,image.shape[1]-buffer],
                                                                        "flux":flux_range})
                star_overlay=make_model_sources_image( image.shape, self.psf, sourcelist)    

                if image.shape==star_overlay.shape:
                    result=self.single_test( image+star_overlay, sourcelist, threshold=2)
                    result["x_0"]+=x_edge
                    result["y_0"]+=y_edge
                    result["x_det"]+=x_edge
                    result["y_det"]+=y_edge
                    test_result=vstack((test_result,result))
                load();load.show()

        
        return test_result

    def single_test(self, image, contains, threshold=2):
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
                        test_result["flux_det"][i]=self.photometry(Table(detections[["x_0","y_0"]][best_match]), image)["flux"]
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

