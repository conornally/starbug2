import numpy as np
from astropy.table import Table
from astropy.io import fits
from starbug2 import routines
from starbug2.artificialstars import Artificial_Stars
from starbug2.utils import export_table
from photutils.psf import FittableImageModel


class Test_Artificial_Stars():
    im=fits.open("tests/dat/image.fits")["SCI"].data
    contains=Table([[50.8, 57.3, 90],
                    [27.0, 73.4, 50],
                    [1.6E-4, 1.8E-5,0]], names=["x_0","y_0","flux_0"])

    def test_single_test(self):
        det=routines.Detection_Routine()
        art=Artificial_Stars(detector=det)

        res=art.single_test(self.im.copy(), self.contains)
        assert type(res)==Table
        assert len(res)==3

    def test_single_testphot(self):
        det=routines.Detection_Routine()
        phot=routines.APPhot_Routine(3,4,5)
        art=Artificial_Stars(detector=det, photometry=phot)
        res=art.single_test(self.im.copy(),self.contains)
        assert type(res)==Table
        assert len(res)==3
        assert np.nansum(res["flux_det"])>0

    def test_create_subimage(self):
        art=Artificial_Stars()
        return 

    def test_run(self):
        psf_model=FittableImageModel(fits.open("tests/dat/psf.fits")[0].data)
        art=Artificial_Stars(detector=routines.Detection_Routine(), photometry=routines.APPhot_Routine(2,3,4), psf=psf_model)
        pos= art.detector(self.im.copy())
        _flux_range=art.photometry(pos, self.im.copy())["flux"]
        flux_range = ( np.nanmin(_flux_range), np.nanmax(_flux_range) )
        print(flux_range)


        assert art.run_auto(self.im.copy(), 10, stars_per_test=1, flux_range=flux_range) is not None
        assert art.run_auto(self.im.copy(), 10, stars_per_test=5, flux_range=flux_range) is not None
        assert art.run_auto(self.im.copy(), 1, stars_per_test=10, flux_range=flux_range) is not None

        assert art.run_auto(self.im.copy(), 10, stars_per_test=5, subimage_size=50, flux_range=flux_range) is not None
        assert art.run_auto(self.im.copy(), 1, stars_per_test=100, subimage_size=50, flux_range=flux_range) is not None

        assert art.run_auto(self.im.copy(), 10, stars_per_test=1, buffer=5, flux_range=flux_range) is not None
        assert art.run_auto(self.im.copy(), 10, stars_per_test=1, buffer=5, subimage_size=20, flux_range=flux_range) is not None

        #assert art.run_auto(self.im.copy(), 10, stars_per_test=1) is not None

        result= art.run_auto(self.im.copy(), 1000, stars_per_test=5, buffer=5, subimage_size=50, flux_range=flux_range)
        assert result is not None
        export_table(result, "/tmp/art.fits")


