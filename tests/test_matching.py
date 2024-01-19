import os,numpy as np
import pytest
from starbug2.matching import GenericMatch, CascadeMatch, BandMatch, parse_mask
from starbug2.utils import import_table
from starbug2.param import load_default_params
from starbug2.bin.main import starbug_main
from astropy.table import Table

@pytest.fixture(autouse=True)
def init():

    if not os.path.exists("tests/dat/image-ap.fits"):
        starbug_main("starbug2 -Ds SIGSRC=10 tests/dat/image.fits".split())
        starbug_main("starbug2 -Ds SIGSRC=3 -otests/dat/image2.fits  tests/dat/image.fits".split())



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

    cat1 = Table( np.array(t1), names=["RA","DEC","flux","eflux"], meta={"FILTER":'a'})
    cat2 = Table( np.array(t2), names=["RA","DEC","flux","eflux"], meta={"FILTER":'b'})
    return [cat1,cat2]

            

class Test_GenericMatch():

    def test_initialsing(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        options = load_default_params()

        m=GenericMatch( )
        assert m.colnames is None
        assert not m.filter 
        assert m.threshold.value == float(options.get("MATCH_THRESH"))
        assert m.verbose == options.get("VERBOSE")

        m=GenericMatch(fltr="MAG", colnames=["RA"], threshold=0.5, verbose=True)
        assert m.colnames == ["RA"]
        assert m.filter == "MAG"
        assert m.threshold.value == 0.5
        assert m.verbose == True

        assert isinstance(m.__str__(), str)

    def test_generic_match1(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        m=GenericMatch()

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
        print(out)
        assert len(out)<=len(cats[0])
        assert len(out)<=len(cats[1])

    def test_generic_match2(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        m=GenericMatch( colnames=["RA"])
        out=m(cats)

        assert out.colnames == ["RA_1","RA_2"]

    def test_finishmatching(self):
        cats = [ import_table(f) for f in ("tests/dat/image-ap.fits", "tests/dat/image2-ap.fits")]
        m=GenericMatch()
        out=m(cats)
        av=m.finish_matching(out)
        

        m=GenericMatch(colnames=["RA","DEC","flux"], fltr="F444W")
        av=m.finish_matching(m.match(cats))
        assert av.colnames==["RA","DEC","flux","stdflux","flag","F444W","eF444W","NUM"]

        m=GenericMatch(colnames=["RA","DEC","flux"])
        for c in cats: del c.meta["FILTER"]
        av=m.finish_matching(m.match(cats))
        assert av.colnames==["RA","DEC","flux","stdflux","flag","MAG","eMAG","NUM"]




    def test_vals(self):
        m=GenericMatch()
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
        print(out)

        assert np.shape(c)==np.shape(out)
        for m in range(len(c)):
            for n in range(len(c[m])):
                a=c[m][n]
                b=out[m][n]
                assert np.isnan(a)==np.isnan(b)
                if not np.isnan(a) or not np.isnan(b):
                    assert a==b

class Test_BandMatch:
    def test_init(self):
        filters = ["a","b","c"]
        m=BandMatch(fltr=filters)
        assert m.filter==["a","b","c"]

    def test_order_catalogue_JWSTMETA(self):
        a = Table(None, meta={"FILTER":'F115W'})
        b = Table(None, meta={"FILTER":'F187N'})
        c = Table(None, meta={"FILTER":'F770W'})

        m=BandMatch()
        assert m.filter==""
        assert m.order_catalogues( [a,c,b] ) == [a,b,c]
        assert m.filter==["F115W","F187N","F770W"]

    def test_order_catalogue_JWSTcolnames(self):
        a = Table(None, names=['F115W'])
        b = Table(None, names=['F187N'])
        c = Table(None, names=['F770W'])

        m=BandMatch()
        assert m.filter==""
        assert m.order_catalogues( [a,c,b] ) == [a,b,c]
        assert m.filter==["F115W","F187N","F770W"]

    def test_order_catalogue_filterMETA(self):
        a = Table(None, meta={"FILTER":'a'})
        b = Table(None, meta={"FILTER":'b'})
        c = Table(None, meta={"FILTER":'c'})

        m=BandMatch(fltr=["a","b","c"])
        assert m.order_catalogues( [a,c,b] ) == [a,b,c]

    def test_order_catalogue_filtercolnames(self):
        a = Table(None, names=['a'])
        b = Table(None, names=['b'])
        c = Table(None, names=['c'])

        m=BandMatch(fltr=["a","b","c"])
        assert m.order_catalogues( [a,c,b] ) == [a,b,c]

    def test_match(self):
        t1=[[1.,1.,1,1,0],
            [2.,2.,2,2,0],
            [3.,3.,3,3,0],
            #[4.,4.,4,4,0],
            ]
        t2=[[1.,1.,1,1,0],
            [2.,2.,2,2,0],
            #[3.,3.,3,3,0],
            [4.,4.,4,4,1],
            ]
        t3=[[1.,1.,1,1,0],
            #[2.,2.,2,2,0],
            #[3.,3.,3,3,0],
            [4.,4.,4,4,2],
            ]

        f=float
        cats = [Table(np.array(t1), names=["RA","DEC","A","NUM","flag"], dtype=[f,f,f,f,np.uint16], meta={"FILTER":"A"}),
                Table(np.array(t2), names=["RA","DEC","B","NUM","flag"], dtype=[f,f,f,f,np.uint16], meta={"FILTER":"B"}),
                Table(np.array(t3), names=["RA","DEC","C","NUM","flag"], dtype=[f,f,f,f,np.uint16], meta={"FILTER":"C"})]


        bm=BandMatch(fltr=["A","B","C"], threshold=[0.1,0.2])
        res=bm(cats)
        print(res)
        assert res.colnames==[ "RA","DEC","NUM","flag", "A","B","C"]
        #res=bm(cats, method="last")
        #print(res)
        res=bm(cats, method="bootstrap")

def test_parsemask():
    table=import_table("tests/dat/image-ap.fits")
    """
    tests=[ "F444W!=np.nan",
            "F444W==np.nan",
            "F444W>0",
            "(F444W>0)&(F444W<20)", ## Dont like this "syntax"
            "F444W+0"
            ]

    for test in tests:
        assert parse_mask(test,table) is not None
    """

def test_matchwithmasks():
    t1=[[0,0,1],
        [1,1,1],
        [2,2,1],
        [3,3,1]]
    t2=[[0,0,1],
        [1,1,1],
        [2,2,0],
        [3,3,1]]
    t3=[[0,0,1],
        [1,1,1],
        [2,2,0],
        [3,3,1]]
    cat1=Table(np.array(t1,float),names=["RA","DEC","a"])
    cat2=Table(np.array(t2,float),names=["RA","DEC","a"])
    cat3=Table(np.array(t3,float),names=["RA","DEC","a"])
    mask=[ np.array([True,True,False,True]), None, np.array([True,True,True,False])]

    res=GenericMatch().match([cat1,cat2,cat3], mask=mask)
    print(res)

        

