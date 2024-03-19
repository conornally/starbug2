from starbug2.bin.afs import afs_main
from starbug2.bin import EXIT_SUCCESS, EXIT_FAIL
from starbug2.artificialstars import Artificial_StarsIII

run = lambda s : afs_main(["starbug2-afs"]+s.split())

def _test_run():
    assert run("tests/dat/image.fits")==EXIT_SUCCESS
    assert run("nope")==EXIT_FAIL
    
