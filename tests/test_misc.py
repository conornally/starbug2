import os, glob
from starbug2 import filters
from starbug2.misc import *

def xtest_init():
    os.environ["STARBUG_DATDIR"]="/tmp/starbug"
    d=os.getenv("STARBUG_DATDIR")
    init_starbug()

    for f in filters:
        assert glob.glob("%s/*%s*"%(d,f))


