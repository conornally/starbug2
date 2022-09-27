import numpy as np
from astropy.table import Table
from starbug2 import utils

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

def test_strnktn():
    assert utils.strnktn( "", 3, 'a') == "aaa"

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

def test_tabppend():
    base=Table( [[0,0], [0,0]], names=('a','b'))
    tab =Table( [[1,1], [1,1]], names=('a', 'b'))
    exp =Table( [[0,0,1,1],[0,0,1,1]], names=('a','b'))
    out=utils.tabppend(base,tab)
    assert np.all(out==exp)

    tab1=tab.copy()
    out=utils.tabppend(None, tab1)
    assert np.all( out==tab) ## tab is not a typo



