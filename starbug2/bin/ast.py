"""StarbugII Artificial Star Testing
usage: starbug2-ast [-vh] [-N ntests] [-n ncores] [-p file.param] [-S nstars] [-s opt=val] image.fits ..
    -h  --help          : show help screen
    -N  --ntests    num : number of tests to run
    -n  --ncores  cores : number of cores to split the tests over
    -p  --param    file : load a parameter file
    -S  --nstars    num : number of stars to inject per test
    -s  --set    option : set parameter at runtime with syntax "-s KEY=VALUE"
    -v  --verbose       : show verbose stdout output
"""

import os,sys,getopt
import numpy as np
from multiprocessing import Pool, Process, shared_memory
from itertools import repeat
from time import sleep
from astropy.io import fits

import starbug2.bin as scr
from starbug2.starbug import StarbugBase
from starbug2.artificialstars import Artificial_StarsIII
from starbug2.utils import printf,perror,warn,export_table,tabppend

VERBOSE =0x01
SHOWHELP=0x02
STOPPROC=0x04
KILLPROC=0x08

_c=np.array([0,0,0], dtype=np.int64)
_share=shared_memory.SharedMemory(create=True, size=_c.nbytes)
buf=np.ndarray(_c.shape, dtype=_c.dtype, buffer=_share.buf)

def load(msg="loading"):
    """
    A loading bar that should be run in a subprocess
    It sits and watches the shared memory buffer and periodically
    prints out a progress bar
    """
    while buf[0]<buf[1]:
        sleep(1)
        p=buf[0]/buf[1]
        msg="recovering:%d%%"%buf[2]
        s="\x1b[2K%s|%-40s|%d/%d\r"%(msg, int(p*40)*'=', buf[0], buf[1])

        printf(s)
        sys.stdout.flush()
    printf("\n")

def afs_parseargv(argv):
    """ Organise the argv line into options, values and arguments """
    options=0
    setopt={"NTESTS":100, "NSTARS":10, "QUIETMODE":1}
    cmd,argv = scr.parsecmd(argv)
    opts,args = getopt.gnu_getopt(argv, "hvN:n:p:S:s:o:", ("help","verbose","ncores=","param=", "set=", "output=",
                                                            "ntests=", "nstars"))

    for opt,optarg in opts:
        if opt in ("-h","--help"): options |= (SHOWHELP|STOPPROC)
        if opt in ("-v","--verbose"): options |= VERBOSE
        if opt in ("-p","--param"): setopt["PARAMFILE"]=optarg
        if opt in ("-n","--ncores"): setopt["NCORES"]=int(optarg)
        if opt in ("-o","--output"): setopt["OUTPUT"]=optarg

        if opt in ("-N","--ntests"): setopt["NTESTS"]=int(optarg)
        if opt in ("-S","--nstars"): setopt["NSTARS"]=int(optarg)

        if opt in ("-s","--set"): 
            if '=' in optarg:
                key,val=optarg.split('=')
                try: val=float(val)
                except: pass
                setopt[key]=val
            else:
                perror("unable to set parameter, use syntax -s KEY=VALUE\n")
                options|=KILLPROC

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
    fname, options, setopt, index = args
    out=None
    if os.path.exists(fname):
        sb=StarbugBase(fname, setopt.get("PARAMFILE"), options=setopt)
        opt=sb.options
        afs=Artificial_StarsIII(sb)
        out=afs.auto_run(opt.get("NTESTS"), stars_per_test=opt.get("NSTARS"),
                mag_range=(opt.get("MAX_MAG"),opt.get("MIN_MAG")), loading_buffer=buf)
    return out

def afs_main(argv):
    options, setopt, args= afs_parseargv(argv)
    exit_code=0

    if options or setopt:
        if (exit_code:=afs_onetimeruns(options, setopt, args)):
            _share.unlink()
            return exit_code

    if args:
        fname=args[0]
        _ntests=setopt.get("NTESTS")
        if options & VERBOSE:
            printf("Artificial Stars\n----------------\n")
            printf("-> loading %s\n"%fname)
            printf("-> running %d tests with %d injections per test\n"%(_ntests,setopt.get("NSTARS")))

        buf[0]=0
        buf[1]=_ntests
        loading=Process(target=load, args=("afs",))
        loading.start()

        if (ncores:=setopt.get("NCORES")) is None or ncores==1:
            setopt["NCORES"]=1
            outs=[fn((fname,options,setopt,0)) for fname in args]
        else:
            ncores=min(ncores,_ntests)
            zip_options=np.full(ncores,options,dtype=int)
            for n in range(ncores):
                if n>0: zip_options[n]&=~VERBOSE
            setopt["NTESTS"]=int(np.ceil(_ntests/ncores))


            pool=Pool(processes=ncores)
            outs=pool.map(fn, zip(repeat(fname), zip_options, repeat(setopt), range(1,ncores+1))) 
            pool.close()

        buf[0]=buf[1] #force finish
        loading.join()

        raw=outs[0]
        for res in outs[1:]: raw=tabppend(raw,res)
        if options & VERBOSE:
            printf("-> compiling results\n")
            printf("-> flux recovery: %.2g\n"%(np.nanmean(raw["flux"]/raw["flux_det"])))

        completeness=Artificial_StarsIII.get_completeness(raw)

        sb=StarbugBase(fname, setopt.get("PARAMFILE"), options=setopt)
        spatial_completeness=Artificial_StarsIII.get_spatialcompleteness(raw,sb.image,res=10)

        results=fits.HDUList([fits.PrimaryHDU(),
                fits.BinTableHDU(data=completeness, name="AST"),
                fits.BinTableHDU(data=raw, name="RAW"),
                fits.ImageHDU(data=spatial_completeness,name="CMP")])
        outdir,bname,_=StarbugBase.sort_output_names(fname, param_output=setopt.get("OUTPUT"))
        if options & VERBOSE: printf("--> %s/%s-ast.fits\n"%(outdir,bname))
        results.writeto("%s/%s-ast.fits"%(outdir,bname),overwrite=True)

    else:
        perror("must include a fits image to work on\n")
        exit_code=scr.EXIT_FAIL

    _share.unlink()
    return exit_code

def ast_mainentry():
    """Command line entry point"""
    return afs_main(sys.argv)
