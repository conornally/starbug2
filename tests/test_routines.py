from starbug2 import routines
from astropy.table import Table, hstack
from astropy.io import fits
from astropy.wcs import WCS

#image=fits.open("tests/dat/image.fits")

#def test_Detection_Routine_none():
#    dt=routines.Detection_Routine()
#    assert dt.find_stars(None) is None
#
#def test_Detection_Routine_crashes():
#    dt=routines.Detection_Routine()
#    out=dt.find_stars(image["SCI"].data)
#    assert out is not None
#
#def test_BackGround_Estimate_Routine_none():
#    bg=routines.BackGround_Estimate_Routine(None)
#    assert bg(None) is None
#
#def test_BackGround_Estimate_Routine_crashes():
#    pass
#
#
#
#
#def test_SourceProperties_none():
#    sp=routines.SourceProperties(None,None)
#    assert sp.calculate_crowding() is None
#    assert sp.calculate_geometry(1)==None
#
#def test_SourceProperties_crashes():
#    ## yeh i know this isnt how youre supposed to test
#    ## im just verifying it doesnt crash
#    print("testing crashes")
#    slist=Table.read("tests/dat/image-ap.fits", format="fits")
#    sp=routines.SourceProperties(image["SCI"].data, slist[["xcentroid","ycentroid"]])
#    a=sp(fwhm=2.3)
#
#    assert "crowding" in a.colnames
#    assert len(a)==len(slist)
#
#    b=sp(fwhm=2.3, do_crowd=0)
#    assert "crowding" not in b.colnames
#    assert len(b)==len(slist)
#
#
#    
#
#
#
#
#    
