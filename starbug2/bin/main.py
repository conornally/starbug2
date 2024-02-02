"""StarbugII - JWST PSF photometry
usage: starbug2 [-ABDfGhMPSv] [-b bgdfile] [-d apfile] [-n ncores] [-o ouput] [-p file.param] [-s opt=val] image.fits ...
   -A  --apphot          : run aperture photometry on a source list
   -B  --background      : run background estimation
   -D  --detect          : run source detection
   -G  --geom            : calculate geometric stats on source list
   -M  --match           : match outputs from all input image files
   -P  --psf             : run psf photometry
   -S  --subbgd          : subtract background from image

   -b  --bgdfile         : load background (-bgd.fits) file
   -d  --apfile  ap.fits : load a source detection (-ap.fits) file to skip the source detection step
   -f  --find            : attempt to find associated -ap -bgd files
   -h  --help            : display uasage information
   -n  --ncores      num : number of CPU cores to split process between
   -o  --output      dir : output directory
   -p  --param   a.param : load parameter file
   -s  --set      option : set value in parameter file at runtime (-s SIGSKY=3)
   -v  --verbose         : display verbose outputs

   --> Single run commands
       --init                     : Initialise Starbug (post install)
       --local-param              : Make a local copy of the default parameter file
       --update-param             : Update an out-of-date local parameter file
       --generate-psf             : Generate a single PSF. Set FILTER, DET_NAME, PSF_SIZE with -s
       --generate-region   a.fits : Make a ds9 region file with a detection file
       --generate-run      *.fits : Generate a simple run script
       --version                  : Print starbug2 version

       --apply-zeropint    a.fits : Apply a zeropoint (-s ZP_MAG=1.0) to a.fits
       --calc-instr-zp     a.fits : Calculate and apply an instrumental zero point onto a.fits

   --> typical runs
      $~ starbug2 -vD -p file.param image.fits      //Source detect on image with a parameter file
      $~ starbug2 -vDM -n4 images*.fits             //Source detect and match outputs of a list of images
      $~ starbug2 -vd image-ap.fits -BP image.fits  //PSF photometry on an image with a source file (image-ap.fits)

See https://starbug2.readthedocs.io for more information.

"""
import os,sys,getopt
sys.stdout.write("\x1b[1mlaunching \x1b[36mstarbug\x1b[0m\n")
from starbug2.utils import *
from starbug2 import param
import starbug2.bin as scr

VERBOSE =0x01
KILLPROC=0x02
STOPPROC=0x04
SHOWHELP=0x08

DODETECT=0x100
DOBGDEST=0x200
DOPHOTOM=0x400
FINDFILE=0x800

DOARTIFL=0x1000
DOMATCH =0x2000
DOAPPHOT=0x4000
DOBGDSUB=0x8000
DOGEOM  =0x10000

GENRATPSF =0x100000
GENRATRUN =0x200000
GENRATREG =0x400000
INITSB    =0x800000
UPDATEPRM =0x1000000
DODEBUG   =0x2000000
CALCINSTZP=0x4000000
APPLYZP   =0x8000000

def starbug_parseargv(argv):
    """Organise the sys argv line into options, values and arguments"""
    options=0
    setopt={}

    cmd,argv=scr.parsecmd(argv)
    opts,args=getopt.gnu_getopt(argv,"ABDfGhMPSvb:d:n:o:p:s:",
            (   "apphot","background", "detect", "find", "geom", "help", "match", "psf", "subbgd", "verbose", "xtest",
                "bgdfile=", "apfile=", "ncores=", "output=", "param=", "set=",
                "init", "generate-psf", "local-param", "generate-region=", "version", "generate-run", "update-param",
                "calc-instr-zp=", "apply-zeropoint=",
                "debug", "dev"))
    for opt,optarg in opts:
        if opt in ("-h","--help"):   options|=(SHOWHELP|STOPPROC)
        if opt in ("-p","--param"):  setopt["PARAMFILE"]= optarg
        if opt in ("-v","--verbose"):options|=VERBOSE

        if opt in ("-A","--apphot"):    options |= DOAPPHOT
        if opt in ("-B","--background"):options |=DOBGDEST
        if opt in ("-D","--detect"):    options |= DODETECT
        if opt in ("-G","--geom"):      options |= DOGEOM
        if opt in ("-M","--match"):     options |= DOMATCH
        if opt in ("-P","--psf"):       options |= DOPHOTOM
        if opt in ("-S","--subbgd"):    options |= DOBGDSUB

        if opt == "--dev": options |= DOARTIFL

        if opt in ("-d","--apfile"):
            if os.path.exists(optarg): setopt["AP_FILE"]=optarg
            else: perror("AP_FILE \"%s\" does not exist\n"%optarg)

        if opt in ("-b","--bgdfile"):
            if os.path.exists(optarg): setopt["BGD_FILE"]=optarg
            else: perror("BGD_FILE \"%s\" does not exist\n"%optarg)

        if opt in ("-f","--find"): options|=FINDFILE
        if opt in ("-n","--ncores"): 
            setopt["NCORES"]=max(1,int(optarg))

        if opt in ("-o","--output"):
            output=optarg
            setopt["OUTPUT"]=optarg

        if opt in ("-s","--set"): 
            if '=' in optarg:
                key,val=optarg.split('=')
                try: val=float(val)
                except: pass
                setopt[key]=val
            else:
                perror("unable to set parameter, use syntax -s KEY=VALUE\n")
                options|=KILLPROC

        if opt=="--init": options|=(INITSB|STOPPROC)
        if opt=="--generate-psf": options|=(GENRATPSF|STOPPROC)
        if opt=="--update-param": options|=(UPDATEPRM|STOPPROC)
        if opt=="--generate-run": options|=(GENRATRUN|STOPPROC)
        if opt=="--generate-region":
            setopt["REGION_TAB"]=optarg
            options|=(GENRATREG|STOPPROC)

        if opt=="--local-param":
            param.local_param()
            printf("--> generating starbug.param\n")
            options|=STOPPROC

        if opt=="--version": 
            printf(starbug2.logo%("starbug2-v%s"%get_version()))
            options|=STOPPROC

        if opt=="--calc-instr-zp":
            setopt["ZP_PSF_CAT"]=optarg
            options|=(CALCINSTZP|APPLYZP|STOPPROC)
        if opt=="--apply-zeropoint":
            setopt["ZP_PSF_CAT"]=optarg
            options|=(APPLYZP|STOPPROC)

    return options,setopt,args

def starbug_onetimeruns(options, setopt, args):
    """
    Options set, verify/run one time functions
    """
    from starbug2.misc import init_starbug, generate_psf, generate_runscript, calc_instrumental_zeropint

    if options&SHOWHELP:
        scr.usage(__doc__,verbose=options&VERBOSE)
        return scr.EXIT_EARLY

    ## Load parameter files for onetime runs
    if (pfile:=setopt.get("PARAMFILE"))==None:
        if os.path.exists("./starbug.param"):pfile="starbug.param"
        else: pfile=None

    init_parameters=param.load_params(pfile)

    if options&UPDATEPRM:
        param.update_paramfile(pfile)
        return scr.EXIT_EARLY

    tmp=param.load_default_params()
    if set(tmp.keys())-set(init_parameters.keys()) | set(init_parameters.keys())-set(tmp.keys()):
        warn("Parameter file version mismatch. Run starbug2 --update-param to update\nquitting :(\n")
        return scr.EXIT_FAIL

    init_parameters.update(setopt)
    if (_output:=init_parameters.get("OUTPUT")): output=_output
    else: output='.'

    #########################
    # One time run commands #
    #########################

    if options&INITSB: ## Initialise or update starbug
        init_starbug()

    if options&GENRATPSF: ## Generate a single PSF
        if (fltr:=init_parameters.get("FILTER")):
            detector=init_parameters.get("DET_NAME")
            psf_size=init_parameters.get("PSF_SIZE")
            printf("Generating PSF: %s %s (%d)\n"%(fltr,detector,psf_size))
            psf=generate_psf(fltr, detector=detector, fov_pixels=psf_size)
            if psf: 
                name="%s%s.fits"%(fltr,"" if detector is None else detector)
                printf("--> %s\n"%name)
                psf.writeto(name, overwrite=True)
            else: perror("PSF Generation failed :(\n")
        else: perror("Unable to generate PSF. Set filter with '-s FILTER=FXXX'\n")

    if options&GENRATRUN: ## Generate a run script
        generate_runscript(args, "starbug2 ")
        if not args: perror("no files included to create runscript with\n")

    if options&GENRATREG: ## Generate a region from a table
        fname=setopt.get("REGION_TAB")
        if fname and os.path.exists(fname):
            table=Table.read(fname,format="fits")
            _,name,_=split_fname(fname)
            export_region(table, colour=init_parameters["REGION_COL"], scale_radius=init_parameters["REGION_SCAL"],
                                 region_radius=init_parameters["REGION_RAD"], xcol=init_parameters["REGION_XCOL"],
                                 ycol=init_parameters["REGION_YCOL"], wcs=init_parameters["REGION_WCS"], fname="%s/%s.reg"%(output,name))
            printf("generating region --> %s/%s.reg\n"%(output,name))

    ###########################
    # instrumental zero point #
    ###########################
    if options&APPLYZP:
        _fname=setopt.get("ZP_PSF_CAT")
        _zp=init_parameters.get("ZP_MAG")
        _std=0
        if _fname and os.path.exists(_fname):
            psftable=Table.read(_fname, format="fits")
            with fits.open(_fname) as fp: _header=fp[1].header ##thats a bit rubbish

            if (fltr:=_header.get("FILTER")) is None:
                if (fltr:=setopt.get("FILTER")) is None:
                    perror("Unable to determine table FILTER: set manually with `-s FILTER=F000W`\n")
                    return scr.EXIT_FAIL

            if options&CALCINSTZP:
                _aptable=None
                _apfname=setopt.get("AP_FILE")
                if _apfname and os.path.exists(_apfname):
                    _aptable=Table.read(_apfname, format="fits")

                    if (res:=calc_instrumental_zeropint(psftable, _aptable, fltr=fltr)) is not None:
                        _zp,_std=res

            if _zp is not None:
                psftable.meta["%s ZEROPOINT"%fltr]=_zp
                psftable.meta["%s eZEROPOINT"%fltr]=_std
                psftable[fltr]=psftable[fltr]+_zp

                dname,fname,_=split_fname(_fname)
                printf("--> %s/%s-zp.fits\n"%(dname,fname))
                export_table( psftable, fname="%s/%s-zp.fits"%(dname,fname), header=_header)
            else: perror("Unable to set ZEROPOINT, set it with -sZP_MAG=000\n")
        else: perror("Unable to locate table \"%s\".\n"%_fname)

    if options&STOPPROC: return scr.EXIT_EARLY ## quiet ending the process if required

    if options&KILLPROC:
        perror("..quitting :(\n\n")
        return scr.usage(__doc__, verbose=options&VERBOSE)

    return scr.EXIT_SUCCESS


def starbug_matchoutputs(starbugs, options, setopt):
    """
    Matching output catalogues
    """
    from starbug2.matching import GenericMatch
    if options&VERBOSE: printf("Matching outputs\n")
    params=param.load_params(setopt.get("PARAMFILE"))
    params.update(setopt)

    if (fname:=combine_fnames( [sb.fname for sb in starbugs] )):
        _,name,_=split_fname(os.path.basename(fname))
        fname="%s/%s"%(starbugs[0].outdir, name)
    else: fname="out"

    header=starbugs[0].header
    #colnames=starbug2.match_cols
    #colnames+=[ name for name in params["MATCH_COLS"].split() if name not in colnames]

    match=GenericMatch( threshold= params["MATCH_THRESH"], colnames=None, pfile=setopt.get("PARAMFILE"))

    if options&(DODETECT|DOAPPHOT):
        full=match( [sb.detections for sb in starbugs], join_type="or")
        av =match.finish_matching(full, num_thresh=params["NEXP_THRESH"], zpmag=params["ZP_MAG"])

        printf("-> %s-ap*...\n"%(fname))
        export_table(full, fname="%s-apfull.fits"%(fname), header=header)
        export_table(av, fname="%s-apmatch.fits"%(fname), header=header)

    if options&DOPHOTOM:
        full=match( [sb.psfcatalogue for sb in starbugs], join_type="or")
        av =match.finish_matching(full, num_thresh=params["NEXP_THRESH"], zpmag=params["ZP_MAG"])

        printf("-> %s-psf*...\n"%(fname))
        export_table(full, fname="%s-psffull.fits"%(fname), header=header)
        export_table(av, fname="%s-psfmatch.fits"%(fname), header=header)



#def fn(fname,options=0,setopt={}):
def fn(args):
    from starbug2.starbug import StarbugBase ## Ive put this here because it takes some time
    sb=None
    fname,options,setopt=args
    if os.path.exists(fname):
        dname,bname,ext=split_fname(fname)

        if options&FINDFILE:
            ap="%s/%s-ap.fits"%(dname,bname)
            bgd="%s/%s-bgd.fits"%(dname,bname)
            if os.path.exists(ap)  and not setopt.get("AP_FILE"):  setopt["AP_FILE"]=ap
            if os.path.exists(bgd) and not setopt.get("BGD_FILE"): setopt["BGD_FILE"]=bgd

        ## Sorting out the stdout
        if options&VERBOSE: 
            printf("-> showing starbug stdout for \"%s\"\n"%fname)
            setopt["VERBOSE"]=1
        elif setopt.get("NCORES")>1: printf("-> hiding starbug stdout for \"%s\"\n"%fname)
        else: printf("-> %s\n"%fname)

        if ext==".fits":
            sb=StarbugBase(fname, pfile=setopt.get("PARAMFILE"), options=setopt)
            if sb.verify(): 
                warn("System verification failed\n")
                return None
                #pass
                #_input=input("Continue with warnings y/N:")
                #if _input=="" or _input not in "yY":
                #    return#quit("..quitting :(")

            if options & DODETECT: sb.detect()
            if options & DOBGDEST: sb.bgd_estimate()
            if options & DOBGDSUB: sb.bgd_subtraction()
            if options & DOGEOM: sb.source_geometry()

            if options & DOAPPHOT: sb.aperture_photometry()
            if options & DOPHOTOM: sb.photometry()

            if options & DOARTIFL: sb.artificial_stars()

        else: perror("file must be type '.fits' not %s\n"%ext)
    else: perror("can't access %s\n"%fname)
    return sb




def starbug_main(argv):
    """Command entry"""
    options, setopt, args= starbug_parseargv(argv)

    if options or setopt: 

        if (exit_code:=starbug_onetimeruns(options,setopt, args)):
            return exit_code

    if args:
        import starbug2
        from multiprocessing import Pool
        from itertools import repeat
        puts(starbug2.logo%starbug2.motd)
        exit_code=scr.EXIT_SUCCESS

        if (ncores:=setopt.get("NCORES")) is None or ncores==1 or len(args)==1:
            setopt["NCORES"]=1
            starbugs=[fn((fname,options,setopt)) for fname in args]
        else:

            zip_options=np.full(len(args),options, dtype=int)
            for n in range(len(args)):
                if n>0: zip_options[n]&=~VERBOSE

            pool=Pool(processes=ncores)
            starbugs=pool.map(fn,zip( args,zip_options,repeat(setopt))) ##
            pool.close()

        for n,sb in enumerate(starbugs): 
            if not sb: 
                perror("FAILED: %s\n"%args[n])
                starbugs.remove(sb)
                exit_code=scr.EXIT_MIXED

        if not starbug2: exit_code=EXIT_FAIL

            
        if options&DOMATCH and len(starbugs)>1:
            starbug_matchoutputs(starbugs, options, setopt)
        

    else:
        perror("fits image file must be included\n")
        exit_code=scr.EXIT_FAIL

    return exit_code

def starbug_mainentry():
    """Entry point"""
    return starbug_main(sys.argv)
