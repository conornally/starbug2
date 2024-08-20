"""StarbugII Artificial Star Testing
usage: starbug2-ast [-vh] [-N ntests] [-n ncores] [-p file.param] [-S nstars] [-s opt=val] image.fits ..
    -h  --help          : show help screen
    -N  --ntests    num : number of tests to run
    -n  --ncores  cores : number of cores to split the tests over
    -o  --output output : output directory or filename to export results to
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
import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

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

        sb=StarbugBase(fname, setopt.get("PARAMFILE"), options=setopt)
        completeness=Artificial_StarsIII.get_completeness(raw)
        _cfit, _compl=Artificial_StarsIII.estim_completeness(completeness)
        head={  "COMPLETE_FN":"F(x)=l/(1+exp(-k(x-xo)))",
                "l":_cfit[0], "k":_cfit[1], "xo":_cfit[2],
                }
                #"COMPLETE 70%": _compl[0],
                #"COMPLETE 50%": _compl[1]}
        #if _compl[0]: printf("-> complete to 90%%: M=%.2f\n"%_compl[0])
        #if _compl[1]: printf("-> complete to 70%%: M=%.2f\n"%_compl[1])
        #if _compl[2]: printf("-> complete to 50%%: M=%.2f\n"%_compl[2])
        for i,frac in enumerate((90,70,50)):
            if _compl[i] and not np.isnan(_compl[i]):
                printf("-> complete to %d%%: %s=%.2f\n"%(frac,sb.filter,_compl[i]))
                head["COMPLETE %d%%"%frac]=_compl[i]

        spatial_completeness=Artificial_StarsIII.get_spatialcompleteness(raw,sb.image,res=10)

        results=fits.HDUList([fits.PrimaryHDU(),
                fits.BinTableHDU(data=completeness, name="AST", header=head),
                fits.BinTableHDU(data=raw, name="RAW"),
                fits.ImageHDU(data=spatial_completeness,name="CMP")])
        outdir,bname,_=StarbugBase.sort_output_names(fname, param_output=setopt.get("OUTPUT"))
        if options & VERBOSE: printf("--> %s/%s-ast.fits\n"%(outdir,bname))
        results.writeto("%s/%s-ast.fits"%(outdir,bname),overwrite=True)

        ## output figure plotting
        if (_fname:=setopt.get("PLOTAST")):
            fig,ax=plt.subplots(1,figsize=(3.5,3),dpi=300)
            ax.scatter(completeness["mag"],completeness["rec"], c='k', lw=0, s=8)
            ax.plot(completeness["mag"],Artificial_StarsIII.scurve(completeness["mag"],*_cfit),c='g',label=r"$f(x)=\frac{%.2f}{1+e^{%.2f(x-%.2f)}}$"%(_cfit[0],-_cfit[1],_cfit[2]))
            ax.axvline(_compl[0], c="seagreen",ls='--', label=("90%%:%.2f"%_compl[0]),lw=0.75)
            ax.axvline(_compl[1], c="seagreen",ls='-.', label=("70%%:%.2f"%_compl[1]),lw=0.75)
            ax.axvline(_compl[2], c="seagreen",ls=':', label=("50%%:%.2f"%_compl[2]),lw=0.75)
            ax.scatter(_compl,(0.9,0.7,0.5),marker='*', c='teal', s=10)
            ax.tick_params(direction="in",top=True,right=True)
            ax.set_title("Artificial Star Test")
            ax.set_xlabel(sb.filter)
            ax.set_ylabel("Fraction Recovered")
            ax.set_yticks([0,.25,.5,.75,1])
            ax.legend(loc="lower left",frameon=False, fontsize=8)
            plt.tight_layout()
            fig.savefig(_fname,dpi=300)
            printf("--> %s\n"%_fname)


    else:
        perror("must include a fits image to work on\n")
        exit_code=scr.EXIT_FAIL

    _share.unlink()
    return exit_code

def ast_mainentry():
    """Command line entry point"""
    return afs_main(sys.argv)
