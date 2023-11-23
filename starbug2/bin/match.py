"""StarbugII Matching 
usage: starbug2-match [-BCGfhv] [-o output] [-p file.param] [-s KEY=VAL] table.fits ...
    -B  --band               : match in "BAND" mode (does not preserve a column for every frame)
    -C  --cascade            : match in "CASCADE" mode (left justify columns)
    -D  --dither             : match in "DITHER" mode (preserves a column for every frame)
    -f  --full               : export full catalogue
    -G  --generic            : match in "GENERIC" mode
    -h  --help               : show help message
    -o  --output  file.fits  : output matched catalogue
    -p  --param   file.param : load starbug parameter file
    -s  --set     option     : set value in parameter file at runtime (-s MATCH_THRESH=1)

    --> typical runs
       $~ starbug2-match -Gfo outfile.fits tab1.fits tab2.fits
       $~ starbug2-match -sMATCH_THRESH=0.2 -sBRIDGE_COL=F444W -Bo out.fits F*W.fits
"""
import os,sys,getopt
import numpy as np
import pkg_resources
from astropy.table import Table, hstack, vstack
from starbug2 import utils
from starbug2 import matching
import starbug2.bin as scr
import starbug2


VERBOSE =0x01
KILLPROC=0x02
STOPPROC=0x04
SHOWHELP=0x08

BANDMATCH   =0x10
DITHERMATCH =0x20
GENERICMATCH=0x40
CASCADEMATCH=0x80

EXPFULL = 0x100


def match_parsemargv():
    options=0
    setopt={}

    cmd,argv=scr.parsecmd(sys.argv)
    opts,args=getopt.getopt(sys.argv[1:], "BCDfGhvo:p:s:", ("band","cascade","dither","full","generic","help","verbose","output=","param=","set="))
    for opt,optarg in opts:
        if opt in ("-h", "--help"):     options|=(SHOWHELP|STOPPROC)
        if opt in ("-v", "--verbose"):  options|=VERBOSE
        if opt in ("-o", "--output"):   setopt["OUTPUT"]=optarg
        if opt in ("-p", "--param"):    setopt["PARAMFILE"]=optarg

        if opt in ("-f","--full"): options|=EXPFULL
        if opt in ("-s","--set"): 
            if '=' in optarg:
                key,val=optarg.split('=')
                try: val=float(val)
                except: pass
                setopt[key]=val
            else: utils.perror("unable to set parameter, use syntax -s KEY=VALUE\n")

        if opt in ("-B","--band"): options|=BANDMATCH
        if opt in ("-C","--cascade"): options|=CASCADEMATCH
        if opt in ("-D","--dither"): options|=DITHERMATCH
        if opt in ("-G","--generic"): options|=GENERICMATCH
    return options, setopt, args

def match_onetimeruns(options, setopt):
    """
    Options set, one time runs
    """
    if options&SHOWHELP:
        scr.usage(__doc__,verbose=options&VERBOSE)
        return scr.EXIT_EARLY

    return scr.EXIT_SUCCESS

def match_fullbandmatch(tables, parameters):
    tomatch={ starbug2.NIRCAM:[], starbug2.MIRI:[] }
    fname=output if output else "out.fits"
    _colnames=["RA","DEC","flag"]

    for i,tab in enumerate(tables):
        fltr=tab.meta.get("FILTER")
        tomatch[starbug2.filters[fltr].instr].append(tab)
        _colnames+=([fltr,"e%s"%fltr])
    
    if tomatch[starbug2.NIRCAM] and tomatch[starbug2.MIRI]:
        utils.printf("Detected NIRCam to MIRI matching\n")
        nircam_matched=matching.band_match(tomatch[starbug2.NIRCAM], colnames=_colnames)
        miri_matched=matching.band_match(tomatch[starbug2.MIRI], colnames=_colnames)

        load=utils.loading(len(miri_matched), msg="Combining NIRCAM-MIRI(%.2g\")"%dthreshold)
        if (bridgecol:=parameters.get("BRIDGE_COL")):
            mask= np.isnan(nircam_matched[bridgecol])
            utils.printf("-> bridging catalogues with %s\n"%bridgecol)
        else: mask=np.full(len(nircam_matched), False)

        matched,_=matching.generic_match((nircam_matched[~mask],miri_matched), threshold=dthreshold, add_src=True, load=load)
        matched.remove_column("NUM")
        matched=vstack((matched, nircam_matched[mask]))
    else:
        matched=matching.band_match(tables, colnames=_colnames)
        
    return matched

def match_main():
    """Entry Point"""
    options,setopt,args = match_parsemargv()
    if options or setopt:
        if (exit_code:=match_onetimeruns(options,setopt))!=scr.EXIT_SUCCESS:
            return exit_code

    ##########
    # PARAMS #
    ##########
    if not (pfile:=setopt.get("PARAMFILE")):
        if os.path.exists("./starbug.param"): pfile="./starbug.param"
        else: pfile="%s/default.param"%pkg_resources.resource_filename("starbug2", "param/")
    parameters=utils.load_params(pfile)
    parameters.update(setopt)

    #################
    # MAIN ROUTINES #
    #################

    tables=[ utils.import_table(fname, verbose=1) for fname in args]
    if tables:
        colnames=starbug2.match_cols
        colnames+=[ name for name in parameters["MATCH_COLS"].split() if name not in colnames]
        dthreshold=parameters["MATCH_THRESH"]
        nthreshold=parameters["NEXP_THRESH"]
        ## snthresh=parameters["SN_THRESH"]

        ## #################
        ## # SN RATIO CUTS #
        ## #################
        ## if snthresh>0:
        ##     utils.puts("SN Ratio Cuts")
        ##     for i,(tab,fltr) in enumerate(zip(tables, filters)):

        ##         if fltr:
        ##             mask = ((tab[fltr]/tab["e%s"%fltr])<snthresh)
        ##             utils.printf("-> %s: Removing %d sources\n"%(fltr, sum(mask)))
        ##             tables[i].remove_rows(mask)
        ##         else:
        ##             utils.perror("Unable to determine filter of \"%s\"\n"%args[i])

        if options & BANDMATCH:
            av=match_fullbandmatch(tables, parameters)
            full=None

        else:
            if options & DITHERMATCH:  av,full=matching.dither_match(tables, threshold=dthreshold, colnames=colnames)
            elif options & CASCADEMATCH: av,full=matching.cascade_match(tables, threshold=dthreshold, colnames=colnames)
            elif options & GENERICMATCH: av,full=matching.generic_match(tables,threshold=dthreshold, add_src=True, average=True, load=options&VERBOSE)
            else:
                options|=EXPFULL
                av,full=matching.generic_match(tables,threshold=dthreshold, add_src=True, average=True, load=options&VERBOSE)

            if av: 
                av.meta.update(tables[0].meta)
                if nthreshold!=-1:
                    mask=av["NUM"]>=nthreshold
                    av=av[mask]

        output=parameters.get("OUTPUT")
        if output is None or output == '.':
            output=utils.combine_fnames( [ name for name in args] , ntrys=100)
        dname,fname,ext=utils.split_fname(output)

        utils.printf("-> %s/%s*\n"%(dname,fname))
        if options&EXPFULL: utils.export_table(full,fname="%s/%sfull.fits"%(dname,fname))
        if av: utils.export_table(av,"%s/%smatch.fits"%(dname,fname))

        return scr.EXIT_SUCCESS
    else:
        utils.perror("No tables loaded for matching.\n")
        return scr.EXIT_FAIL
