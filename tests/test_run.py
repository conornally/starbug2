import os,glob
import pytest
from starbug2.utils import wget
from starbug2.bin import EXIT_SUCCESS, EXIT_EARLY, EXIT_FAIL, EXIT_MIXED
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
    assert run("starbug2 -p starbug.param tests/dat/image.fits")==EXIT_SUCCESS
    clean()
    #run("rm starbug.param")

def test_detect():
    clean()
    #wget("https://app.box.com/s/f2trqcln5mjug3rigs9246202ztuu5fh", "image.fits")
    assert run("starbug2 -v tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -D tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 --detect tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -D -sSIGSKY=3 -sSIGSRC=15 tests/dat/image.fits")==EXIT_SUCCESS
    clean()
    #run("rm image*.fits")

def test_bgd():
    clean()
    assert run("starbug2 -D tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -d tests/dat/image-ap.fits -B tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -d tests/dat/image-ap.fits --background tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -vf -B tests/dat/image.fits")==EXIT_SUCCESS
    clean()

def test_psf():
    clean()
    assert run("starbug2 -DB tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -d tests/dat/image-ap.fits -b tests/dat/image-bgd.fits -P tests/dat/image.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    assert run("starbug2 -fP tests/dat/image.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    assert run("starbug2 -d tests/dat/image-ap.fits -P tests/dat/image.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    assert run("starbug2 -fBP tests/dat/image.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    assert run("starbug2 -fPs GEN_RESIDUAL=1 tests/dat/image.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    clean()

def test_residual():
    clean()
    assert run("starbug2 -DB tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 -fSs GEN_RESIDUAL=1 tests/dat/image.fits")==EXIT_SUCCESS
    assert run("starbug2 tests/dat/image-res.fits")==EXIT_SUCCESS
    assert run("starbug2 -D tests/dat/image-res.fits")==EXIT_SUCCESS
    assert run("starbug2 -fB tests/dat/image-res.fits")==EXIT_SUCCESS
    assert run("starbug2 -fP tests/dat/image-res.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    assert run("starbug2 -fPs GEN_RESIDUAL=1 tests/dat/image-res.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS

    assert run("starbug2 -fSA tests/dat/image.fits")==EXIT_SUCCESS
    clean()


def test_ncores():
    clean()
    os.system("cp tests/dat/image.fits tests/dat/image2.fits")
    assert run("starbug2 tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS
    assert run("starbug2 -n2 tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS
    assert run("starbug2 -vD tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS
    assert run("starbug2 -Dn0 tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS
    assert run("starbug2 -Dn1 tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS
    assert run("starbug2 -Dn2 tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS
    assert run("starbug2 -Dn4 tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS

    assert run("starbug2 -DBP tests/dat/image.fits tests/dat/image2.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    assert run("starbug2 -vDBPn2 tests/dat/image.fits tests/dat/image2.fits -sPSF_FILE=tests/dat/psf.fits")==EXIT_SUCCESS
    assert run("starbug2 -DM tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS
    assert run("starbug2 -DMn2 tests/dat/image.fits tests/dat/image2.fits")==EXIT_SUCCESS

    assert run("starbug2 -D tests/dat/image-ap.fits tests/dat/image.fits")==EXIT_MIXED
    assert run("starbug2 -D bad.fits tests/dat/image.fits")==EXIT_MIXED
    clean()









def clean():
    files=glob.glob("tests/dat/*")
    files.remove("tests/dat/image.fits")
    files.remove("tests/dat/psf.fits")
    for fname in files: os.remove(fname)
    if os.path.exists("starbug.param"): os.remove("starbug.param")


