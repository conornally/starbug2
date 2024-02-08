from starbug2.starbug import StarbugBase

class Test_Implementation:
    image="tests/dat/image.fits"
    
    def test_artificial_stars(self):

        sb=StarbugBase(self.image, options={"SUBIMAGE":50, "NTESTS":100, "NSTARS":10, "VERBOSE":1})
        sb._image["SCI"].data*=0
        status= sb.artificial_stars()
        assert status==0
