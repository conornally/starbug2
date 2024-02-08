import starbug2
from starbug2 import utils
import numpy as np
from astropy.table import Table, MaskedColumn
from astropy.io import fits

def test_strnktn():
    assert utils.strnktn( "", 3, 'a') == "aaa"
    assert utils.strnktn( "a", 3, 'a') == "aaaa"
    assert utils.strnktn( "a", 0, 'a') == "a"

def test_split_fname():
    fname="/path/to/file.fits"
    d,f,e = utils.split_fname(fname)
    assert d=="/path/to"
    assert f=="file"
    assert e==".fits"

    fname="file.fits"
    d,f,e = utils.split_fname(fname)
    assert d=="."
    assert f=="file"
    assert e==".fits"

    fname="file"
    d,f,e = utils.split_fname(fname)
    assert d=="."
    assert f=="file"
    assert e==""


def test_flux2mag():
    ## input shape
    assert len( utils.flux2mag(1, None, zp=1)[0]) ==1
    assert len( utils.flux2mag(np.ones(10), None, zp=1)[0]) ==10
    assert len( utils.flux2mag(np.full(10,np.nan), None, zp=1)[0]) ==10
    a,b=utils.flux2mag( np.empty(10), np.empty(10), zp=1)
    assert len(a)==len(b)
    a,b=utils.flux2mag( 1, 1, zp=1)
    assert len(a)==len(b)
    a,b=utils.flux2mag( 0, 0, zp=1)
    assert len(a)==len(b)

    ## normal fluxed
    flux=np.array(    [1, 100, 999, 123, 3.4, 87654, np.pi] )
    fluxerr=None
    mag,magerr=utils.flux2mag(flux,fluxerr, zp=1)
    assert np.all(np.equal(mag , -2.5*np.log10(flux)))

    ## boundary fluxes
    flux=np.array( [0, 0.0, -1, np.nan] )
    fluxerr=None
    mag,magerr=utils.flux2mag(flux,fluxerr,zp=1)
    assert np.isnan(mag).all()
    assert utils.flux2mag( np.inf )[0] == -np.inf ##should be -inf
    assert np.isnan(utils.flux2mag( -np.inf )[0]) ##Should be nan

    ##fluxerr
    flux=np.array( [1234, 1, 0.00001, 10])
    fluxerr=np.array( [1,100,123456,1.234567] )
    mag,magerr=utils.flux2mag( np.ones(flux.shape), fluxerr,zp=1)
    assert np.all(np.equal(magerr, 2.5*np.log10( 1.0+( fluxerr/np.ones(flux.shape)) ))) ##flux all 1
    mag,magerr=utils.flux2mag( flux,fluxerr,zp=1)
    assert np.all(np.equal(magerr, 2.5*np.log10( 1.0+( fluxerr/flux) ))) ## random fluxes

    ## boundary fluxerrs
    assert utils.flux2mag(1, None, zp=1)[1] ==0
    assert np.isnan(utils.flux2mag(1, np.nan, zp=1)[1])
    assert np.isnan(utils.flux2mag(1, -1, zp=1)[1])

    ##ZPs


def test_find_colnames():
    tab=Table(None, names=["A", "word", "word1", "word2", "notword", "_word"])
    res=utils.find_colnames(tab, "word")

    assert res is not None
    assert res == ["word", "word1", "word2"]
    assert utils.find_colnames(tab, "badmatch")==[]


def test_tabppend():
    base=Table( [[0,0], [0,0]], names=('a','b'))
    tab =Table( [[1,1], [1,1]], names=('a', 'b'))
    exp =Table( [[0,0,1,1],[0,0,1,1]], names=('a','b'))
    out=utils.tabppend(base,tab)
    assert np.all(out==exp)

    tab1=tab.copy()
    out=utils.tabppend(None, tab1)
    assert np.all( out==tab) ## tab is not a typo


def test_parse_unit():
    assert utils.parse_unit("10p") == (10, starbug2.PIX)
    assert utils.parse_unit("10s") == (10, starbug2.ARCSEC)
    assert utils.parse_unit("10m") == (10, starbug2.ARCMIN)
    assert utils.parse_unit("10d") == (10, starbug2.DEG)

    assert utils.parse_unit("10.1s") == (10.1, starbug2.ARCSEC)
    assert utils.parse_unit("-10.1s") == (-10.1, starbug2.ARCSEC)
    assert utils.parse_unit("0s") == (0, starbug2.ARCSEC)
    assert utils.parse_unit("0") == (0, None)

    assert utils.parse_unit("") == (None, None)
    assert utils.parse_unit("p") == (None, None)

def test_rmduplicates():
    lst=["a","b","b","c","b","c"]
    lst2=utils.rmduplicates(lst)
    assert lst2==["a","b","c"]

    assert utils.rmduplicates([]) == []
    assert utils.rmduplicates(["a"]) == ["a"]

def test_hcascade():
    t1=[[1,1,0],
        [2,2,0],
        [3,3,0],
        #[4,4,0]
        ]
    t2=[[1,1,0],
        [2,2,0],
        [3,3,1],
        [4,4,0]
        ]

    tables=[Table(np.array(t1), names=["A","B","flag"], dtype=[float,float,np.uint16]),
            Table(np.array(t2), names=["A","B","flag"], dtype=[float,float,np.uint16])]
    nan=MaskedColumn(None,dtype=float).info.mask_val
    nan=np.ma.masked
    nan=np.nan
    res=utils.hcascade(tables)
    test=Table( np.ma.array([ [1,1,0,1,1,0],
                            [2,2,0,2,2,0],
                            [3,3,0,3,3,1],
                            [4,4,0,nan,nan,0]]), 

                            dtype=[float,float,np.uint16,float,float,np.uint16], 
                            names=["A_1","B_1","flag_1","A_2","B_2","flag_2"])

    res=utils.fill_nan(res)
    assert np.shape(res)==np.shape(test)
    for m in range(len(res)):
        for n in range(len(res[m])):
            a=res[m][n]
            b=test[m][n]
            assert np.isnan(a)==np.isnan(b)
            if not np.isnan(a) or not np.isnan(b):
                assert a==b

def test_collapseheader():
    header=fits.Header( {"OK":0,
                         "PARAMFILE":"/PATH/TO/FILE/THAT/IS/TOO/LONG/FOR/A/HIERARCH/CARD",
                        "PARAMFILE2":"/PATH/TO/FILE/THAT/IS/TOO/LONG/FOR/A/HIERARCH/CARD"})

    h=utils.collapse_header(header)
    assert h["COMMENT"] is not None
    assert type(utils.collapse_header( {"a":"b"} ))==fits.Header
