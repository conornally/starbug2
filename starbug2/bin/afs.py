"""StarbugII Artificial Star Testing
usage: starbug2-afs [-vh] [-p file.param] [-s opt=val] image.fits ..
    -h  --help          : show help screen
    -p  --param    file : load a parameter file
    -s  --set    option : set parameter at runtime with syntax "-s KEY=VALUE"
    -v  --verbose       : show verbose stdout output
"""

import os,sys,getopt
import starbug2.bin as scr
from starbug2.starbug import StarbugBase
from starbug2.artificialstars import Artificial_StarsIII
from starbug2.utils import printf,perror,warn,export_table

VERBOSE =0x01
SHOWHELP=0x02
STOPPROC=0x04
KILLPROC=0x08

def afs_parseargv(argv):
    """ Organise the argv line into options, values and arguments """
    options=0
    setopt={}
    cmd,argv = scr.parsecmd(argv)
    opts,args = getopt.gnu_getopt(argv, "hvp:s:o:", ("help","verbose","param=", "set=", "output="))

    for opt,optarg in opts:
        if opt in ("-h","--help"): options |= (SHOWHELP|STOPPROC)
        if opt in ("-p","--param"): setopt["PARAMFILE"]=optarg
        if opt in ("-v","--verbose"): options |= VERBOSE

        if opt in ("-s","--set"): 
            if '=' in optarg:
                key,val=optarg.split('=')
                try: val=float(val)
                except: pass
                setopt[key]=val
            else:
                perror("unable to set parameter, use syntax -s KEY=VALUE\n")
                options|=KILLPROC

        if opt in ("-o","--output"):
            warn("argument \"%s\" not implemented yet\n"%opt)

    return options, setopt, args

def afs_onetimeruns(options, setopt, args):
    """
    Set options, verify run and execute one time functions
    """

    if options&SHOWHELP:
        scr.usage(__doc__,verbose=options&VERBOSE)
        return scr.EXIT_EARLY

    if options & STOPPROC: return scr.EXIT_EARLY
    if options & KILLPROC: 
        perror("..killing process\n")
        return scr.EXIT_FAIL

    return scr.EXIT_SUCCESS

def fn(args):
    return None

def afs_main(argv):
    options, setopt, args= afs_parseargv(argv)
    exit_code=0

    if options or setopt:
        if (exit_code:=afs_onetimeruns(options, setopt, args)):
            return exit_code
    if args:
        setopt["QUIETMODE"]=1
        sb=StarbugBase(args[0], setopt.get("PARAMFILE"), options=setopt)
        afs=Artificial_StarsIII(sb)
        a=afs.auto_run(sb.options.get("NTESTS"), stars_per_test=sb.options.get("NSTARS"))
        export_table(a,"/tmp/out-afs.fits")



    else:
        perror("must include a fits image to work on\n")
        exit_code=scr.EXIT_FAIL

    return exit_code

def afs_mainentry():
    """Command line entry point"""
    return afs_main(sys.argv)
