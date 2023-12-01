import os,numpy as np
import pytest
from starbug2.matching import Matcher, CascadeMatch
from starbug2.utils import import_table
from starbug2.param import load_default_params
from astropy.table import Table

@pytest.fixture(autouse=True)
def init():

    if not os.path.exists("tests/dat/image-ap.fits"):
        os.system("starbug2 -Ds SIGSRC=10 tests/dat/image.fits")
        os.system("starbug2 -Ds SIGSRC=3 -otests/dat/image2.fits  tests/dat/image.fits")



def cats():
    t1=[[ 0.0, 0.0, 1.0, 0.1],
        [ 0.1, 0.1, 1.0, 0.1],
        [ 0.2, 0.1, 2.0, 0.2],
        [ 0.1, 0.2, 2.0, 0.2],
        [ 1.0, 1.0, 100, 0.1],
        [ 1.1, 1.0, 100, 0.1],
        #[ 2.0, 2.0, 201, 0.1],
        ]

    t2=[[ 0.0, 0.0, 1.1, 0.1],
        #[ 0.1, 0.1, 1.0, 0.1],
        [ 0.2, 0.1, 2.1, 0.2],
        [ 0.1, 0.2, 2.1, 0.2],
        [ 1.0, 1.0, 101, 0.1],
        [ 1.1, 1.0, 101, 0.1],
        [ 2.0, 2.0, 201, 0.1],
        ]

    cat1 = Table( np.array(t1), names=["RA","DEC","flux","eflux"])
    cat2 = Table( np.array(t2), names=["RA","DEC","flux","eflux"])
    return [cat1,cat2]

            

class Test_Matcher():

    def test_initialsing(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        options = load_default_params()

        m=Matcher( )
        assert m.colnames is None
        assert not m.filter 
        assert m.threshold.value == float(options.get("MATCH_THRESH"))
        assert m.verbose == options.get("VERBOSE")

        m=Matcher(fltr="MAG", colnames=["RA"], threshold=0.5, verbose=True)
        assert m.colnames == ["RA"]
        assert m.filter == "MAG"
        assert m.threshold.value == 0.5
        assert m.verbose == True

        assert isinstance(m.__str__(), str)

    def test_generic_match1(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        m=Matcher()

        out=m(cats)
        assert isinstance(out, Table)
        for name in cats[0].colnames:
            if name != "Catalogue_Number":
                assert "%s_1"%name in out.colnames
                assert "%s_2"%name in out.colnames
        assert len(out)>=len(cats[0])
        assert len(out)>=len(cats[1])
        assert m.filter=="F444W"

        out=m(cats, join_type="and")
        assert len(out)<=len(cats[0])
        assert len(out)<=len(cats[1])

    def test_generic_match2(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        m=Matcher( colnames=["RA"])
        out=m(cats)

        assert out.colnames == ["RA_1","RA_2"]

    def test_finishmatching(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        m=Matcher()
        out=m(cats)
        av=m.finish_matching(out)
        

        m=Matcher(colnames=["RA","DEC","flux"], fltr="F444W")
        av=m.finish_matching(m.match(cats))
        assert av.colnames==["RA","DEC","flux","stdflux","flag","F444W","eF444W","NUM"]

        m=Matcher(colnames=["RA","DEC","flux"])
        for c in cats: del c.meta["FILTER"]
        av=m.finish_matching(m.match(cats))
        assert av.colnames==["RA","DEC","flux","stdflux","flag","MAG","eMAG","NUM"]




    def test_vals(self):
        m=Matcher()
        out=m(cats())
        t=[[ 0.0, 0.0, 1.0, 0.1,   0.0, 0.0, 1.1, 0.1],
           [ 0.1, 0.1, 1.0, 0.1,   np.nan, np.nan, np.nan, np.nan],
           [ 0.2, 0.1, 2.0, 0.2,   0.2, 0.1, 2.1, 0.2],
           [ 0.1, 0.2, 2.0, 0.2,   0.1, 0.2, 2.1, 0.2],
           [ 1.0, 1.0, 100, 0.1,   1.0, 1.0, 101, 0.1],
           [ 1.1, 1.0, 100, 0.1,   1.1, 1.0, 101, 0.1],
           [np.nan,np.nan,np.nan,np.nan, 2.0, 2.0, 201, 0.1] ]
        c=Table(np.array(t), names=[ "RA_1","DEC_1","flux_1","eflux_1","RA_2","DEC_2","flux_2","eflux_2"])
        
        assert np.shape(c)==np.shape(out)
        for m in range(len(c)):
            for n in range(len(c[m])):
                a=c[m][n]
                b=out[m][n]
                assert np.isnan(a)==np.isnan(b)
                if not np.isnan(a) or not np.isnan(b):
                    assert a==b


class Test_Cascade:
    def test_cascadematch(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        m=CascadeMatch()
        out=m.match(cats)
    def test_vals(self):
        t=[[ 0.0, 0.0, 1.0, 0.1,   0.0, 0.0, 1.1, 0.1],
           [ 0.1, 0.1, 1.0, 0.1,   np.nan, np.nan, np.nan, np.nan],
           [ 0.2, 0.1, 2.0, 0.2,   0.2, 0.1, 2.1, 0.2],
           [ 0.1, 0.2, 2.0, 0.2,   0.1, 0.2, 2.1, 0.2],
           [ 1.0, 1.0, 100, 0.1,   1.0, 1.0, 101, 0.1],
           [ 1.1, 1.0, 100, 0.1,   1.1, 1.0, 101, 0.1],
           [2.0, 2.0, 201, 0.1,    np.nan,np.nan,np.nan,np.nan] ]
        c=Table(np.array(t), names=[ "RA_1","DEC_1","flux_1","eflux_1","RA_2","DEC_2","flux_2","eflux_2"])
        m=CascadeMatch()
        out=m.match(cats())

        assert np.shape(c)==np.shape(out)
        for m in range(len(c)):
            for n in range(len(c[m])):
                a=c[m][n]
                b=out[m][n]
                assert np.isnan(a)==np.isnan(b)
                if not np.isnan(a) or not np.isnan(b):
                    assert a==b


