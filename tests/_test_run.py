from os import system as run
from starbug2.utils import wget


def test_start():
    assert run("starbug2 -h")==0
    assert run("starbug2 -vh")==0
    assert run("starbug2 --version")==0
    assert run("starbug2")==0

def test_detect():
    #wget("https://app.box.com/s/f2trqcln5mjug3rigs9246202ztuu5fh", "image.fits")
    assert run("starbug2 -v dat/image.fits")==0
    assert run("starbug2 -D dat/image.fits")==0
    assert run("starbug2 --detect dat/image.fits")==0
    assert run("starbug2 -D -sSIGSKY=3 -sSIGSRC=15 dat/image.fits")==0
    #run("rm image*.fits")

def test_bgd():
    assert run("starbug2 -d dat/image-ap.fits -B dat/image.fits")==0
    assert run("starbug2 -d dat/image-ap.fits --background dat/image.fits")==0
    assert run("starbug2 -vf -B dat/image.fits")==0

def test_psf():
    assert run("starbug2 -d dat/image-ap.fits -b dat/image-bgd.fits -P dat/image.fits")==0
    assert run("starbug2 -fP dat/image.fits")==0
    assert run("starbug2 -fBP dat/image.fits")==0
    assert run("starbug2 -fPs GEN_RESIDUAL=1 dat/image.fits")==0

def test_2ndrun():
    assert run("starbug2 -v dat/image-res.fits")==0
    assert run("starbug2 -D dat/image-res.fits")==0
    assert run("starbug2 -DB dat/image-res.fits")==0
    assert run("starbug2 -DBP dat/image-res.fits")==0
    assert run("starbug2 -vfP dat/image-res.fits")==0




