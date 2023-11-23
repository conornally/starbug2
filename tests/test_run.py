import os,glob
from starbug2.utils import wget
from starbug2.bin import EXIT_SUCCESS, EXIT_EARLY, EXIT_FAIL
from starbug2.bin.main import starbug_main

run = lambda s:starbug_main(s.split())


def test_start():
    clean()
    assert run("starbug2 -h")==EXIT_EARLY
    assert run("starbug2 -vh")==EXIT_EARLY
    assert run("starbug2 --version")==EXIT_EARLY
    assert run("starbug2")==EXIT_FAIL

def test_param():
    clean()
    assert run("starbug2 --local-param")==EXIT_EARLY
    assert run("starbug2 --update-param")==EXIT_EARLY
    assert run("starbug2 -p starbug.param dat/image.fits")==EXIT_SUCCESS
    #run("rm starbug.param")

def test_detect():
    clean()
    #wget("https://app.box.com/s/f2trqcln5mjug3rigs9246202ztuu5fh", "image.fits")
    assert run("starbug2 -v dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -D dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 --detect dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -D -sSIGSKY=3 -sSIGSRC=15 dat/image.fits")==EXIT_SUCCESS
    #run("rm image*.fits")

def test_bgd():
    clean()
    assert run("starbug2 -D dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -d dat/image-ap.fits -B dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -d dat/image-ap.fits --background dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -vf -B dat/image.fits")==EXIT_SUCCESS

def test_psf():
    clean()
    assert run("starbug2 -DB dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -d dat/image-ap.fits -b dat/image-bgd.fits -P dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -fP dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -fBP dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -fPs GEN_RESIDUAL=1 dat/image.fits")==EXIT_SUCCESS

def clean():
    files=glob.glob("dat/*")
    print(files)
    files.remove("dat/image.fits")
    for fname in files: os.remove(fname)



