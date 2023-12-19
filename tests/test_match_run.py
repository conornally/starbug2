import os,glob
import pytest
from starbug2.utils import wget
from starbug2.bin import EXIT_SUCCESS, EXIT_EARLY, EXIT_FAIL
from starbug2.bin.main import starbug_main
from starbug2.bin.match import match_main
run = lambda s:match_main(s.split())

def test_match_start():
    assert run("starbug2-match")==EXIT_FAIL
    assert run("starbug2-match -h")==EXIT_EARLY
    assert run("starbug2-match -vh")==EXIT_EARLY

def test_match_badinput():
    #clean()
    assert run("starbug2-match ")==EXIT_FAIL
    assert run("starbug2-match tests/dat/image.fits")==EXIT_EARLY
    assert run("starbug2-match badinput.fits")==EXIT_FAIL
    assert run("starbug2-match badinput.txt")==EXIT_FAIL
    #assert run("starbug2-match tests/dat/image.fits tests/dat/image.fits")==EXIT_FAIL
    starbug_main("starbug2 -D tests/dat/image.fits".split())
    assert run("starbug2-match tests/dat/image-ap.fits")==EXIT_EARLY
    
    #clean()

def test_match_basicrunthrough():
    #clean()
    starbug_main("starbug2 -Do tests/dat/out1.fits tests/dat/image.fits".split())
    starbug_main("starbug2 -Do tests/dat/out2.fits tests/dat/image.fits".split())
    assert run("starbug2-match tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS
    assert run("starbug2-match -G tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS
    assert run("starbug2-match -C tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS
    #assert run("starbug2-match -D tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS

    assert run("starbug2-match -f tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS
    assert run("starbug2-match -fG tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS
    assert run("starbug2-match -fC tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS
    #assert run("starbug2-match -fD tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS
    #clean()

def test_mask():
    starbug_main("starbug2 -Do tests/dat/out1.fits tests/dat/image.fits".split())
    starbug_main("starbug2 -Do tests/dat/out2.fits tests/dat/image.fits".split())
    assert run("starbug2-match -vmF444W>20 tests/dat/out1-ap.fits tests/dat/out2-ap.fits")==EXIT_SUCCESS



@pytest.fixture(autouse=True)
def init():
    files=glob.glob("tests/dat/*")
    files.remove("tests/dat/image.fits")
    files.remove("tests/dat/psf.fits")
    for fname in files: os.remove(fname)
    if os.path.exists("starbug.param"): os.remove("starbug.param")



