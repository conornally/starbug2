"""StarbugII Matching 
usage: starbug2-match [-BCGfhv] [-e column] [-m mask] [-o output] [-p file.param] [-s KEY=VAL] table.fits ...
    -B  --band               : match in "BAND" mode (does not preserve a column for every frame)
    -C  --cascade            : match in "CASCADE" mode (left justify columns)
    -G  --generic            : match in "GENERIC" mode

    -e  --error   column     : photometric error column ("eflux" or "stdflux")
    -f  --full               : export full catalogue
    -h  --help               : show help message
    -m  --mask    eval       : column evaluation to mask out of matching e.g. -m F444W==nan
    -o  --output  file.fits  : output matched catalogue
    -p  --param   file.param : load starbug parameter file
    -s  --set     option     : set value in parameter file at runtime (-s MATCH_THRESH=1)
    -v  --verbose            : display verbose outputs

        --band-depr          : match in "old" band mode

    --> typical runs
       $~ starbug2-match -Gfo outfile.fits tab1.fits tab2.fits
       $~ starbug2-match -sMATCH_THRESH=0.2 -sBRIDGE_COL=F444W -Bo out.fits F*W.fits
"""
import os,sys,getopt
import numpy as np
from astropy.table import Table, hstack, vstack
from starbug2 import utils
from starbug2.matching import GenericMatch, CascadeMatch, BandMatch, band_match, parse_mask
from starbug2 import param
import starbug2.bin as scr
import starbug2


VERBOSE =0x01
KILLPROC=0x02
STOPPROC=0x04
SHOWHELP=0x08

BANDMATCH   =0x10
BANDDEPR    =0x20
GENERICMATCH=0x40
CASCADEMATCH=0x80

EXPFULL = 0x100


def match_parsemargv(argv):
    options=0
    setopt={}

    cmd,argv=scr.parsecmd(argv)
    opts,args=getopt.gnu_getopt(argv, "BCfGhve:m:o:p:s:", ("band","cascade","dither","full","generic","help","verbose",
                                                    "error=","mask=","output=","param=","set=",
                                                    "band-depr"))
    for opt,optarg in opts:
        if opt in ("-h", "--help"):     options|=(SHOWHELP|STOPPROC)
        if opt in ("-v", "--verbose"):  options|=VERBOSE
        if opt in ("-o", "--output"):   setopt["OUTPUT"]=optarg
        if opt in ("-p", "--param"):    setopt["PARAMFILE"]=optarg

        if opt in ("-e","--error"): setopt["ERRORCOLUMN"]=optarg
        if opt in ("-f","--full"): options|=EXPFULL
        if opt in ("-m","--mask"): setopt["MASKEVAL"]=optarg
        if opt in ("-s","--set"): 
            if '=' in optarg:
                key,val=optarg.split('=')
                try: val=float(val)
                except: pass
                setopt[key]=val
            else: utils.perror("unable to set parameter, use syntax -s KEY=VALUE\n")


        if opt in ("-B","--band"): options|=BANDMATCH
        if opt in ("-C","--cascade"): options|=CASCADEMATCH
        if opt in ("-G","--generic"): options|=GENERICMATCH
        if opt == "--band-depr": options|=BANDDEPR
    return options, setopt, args

def match_onetimeruns(options, setopt):
    """
    Options set, one time runs
    """
    if options&VERBOSE: setopt["VERBOSE"]=1
    if options&SHOWHELP:
        scr.usage(__doc__,verbose=options&VERBOSE)
        return scr.EXIT_EARLY


    return scr.EXIT_SUCCESS

def match_fullbandmatch(tables, parameters):
    utils.perror("THIS NEEDS A TEST\n")
    tomatch={ starbug2.NIRCAM:[], starbug2.MIRI:[] }
    _colnames=["RA","DEC","flag"]
    dthreshold=parameters.get("MATCH_THRESH")

    for i,tab in enumerate(tables):
        fltr=tab.meta.get("FILTER")
        tomatch[starbug2.filters[fltr].instr].append(tab)
        _colnames+=([fltr,"e%s"%fltr])
    
    if tomatch[starbug2.NIRCAM] and tomatch[starbug2.MIRI]:
        utils.printf("Detected NIRCam to MIRI matching\n")
        nircam_matched=band_match(tomatch[starbug2.NIRCAM], colnames=_colnames)
        miri_matched=band_match(tomatch[starbug2.MIRI], colnames=_colnames)

        load=utils.loading(len(miri_matched), msg="Combining NIRCAM-MIRI(%.2g\")"%dthreshold)
        if (bridgecol:=parameters.get("BRIDGE_COL")):
            mask= np.isnan(nircam_matched[bridgecol])
            utils.printf("-> bridging catalogues with %s\n"%bridgecol)
        else: mask=np.full(len(nircam_matched), False)

        m=GenericMatch(threshold=dthreshold, load=load)
        full=m((nircam_matched[~mask],miri_matched))
        matched = m.finish_matching(full)
        #matched,_=matching.generic_match((nircam_matched[~mask],miri_matched), threshold=dthreshold, add_src=True, load=load)
        matched.remove_column("NUM")
        matched=vstack((matched, nircam_matched[mask]))
    else:
        matched=band_match(tables, colnames=_colnames)
        
    return matched

def match_main(argv):
    """"""
    options,setopt,args = match_parsemargv(argv)
    if options or setopt:
        if (exit_code:=match_onetimeruns(options,setopt))!=scr.EXIT_SUCCESS:
            return exit_code

    ##########
    # PARAMS #
    ##########
    if not (pfile:=setopt.get("PARAMFILE")):
        if os.path.exists("./starbug.param"): pfile="./starbug.param"
        else: pfile=None
    parameters=param.load_params(pfile)
    parameters.update(setopt)

    #################
    # MAIN ROUTINES #
    #################
    exit_code=scr.EXIT_SUCCESS

    tables=[ ]
    for fname in args: 
        t=utils.import_table(fname, verbose=1)
        if t is not None: tables.append(t)
    if (raw:=parameters.get("MASKEVAL")):
        masks = [ parse_mask(raw,t) for t in tables ]
        for m in masks:
            try: print(m, sum(m), len(m))
            except: print( m )
    else: masks=None


    if len(tables)>1:
        colnames=starbug2.match_cols
        colnames+=[ name for name in parameters["MATCH_COLS"].split() if name not in colnames]
        dthreshold=parameters["MATCH_THRESH"]
        error_column = setopt.get("ERRORCOLUMN") if setopt.get("ERRORCOLUMN") else "eflux"

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

        if options & BANDDEPR:
            av=match_fullbandmatch(tables, parameters)
            full=None
            options&=~EXPFULL

        else:
            if options & BANDMATCH:
                if isinstance(dthreshold, str):
                    dthreshold=np.array(parameters["MATCH_THRESH"].split(','), float)
                if parameters["FILTER"]!="":
                    fltr=parameters["FILTER"].split(',')
                else: fltr=utils.rmduplicates( [utils.find_filter(t) for t in tables] )
                matcher=BandMatch(threshold=dthreshold, fltr=fltr, verbose=parameters["VERBOSE"])

            elif options & CASCADEMATCH: matcher=CascadeMatch(threshold=dthreshold, colnames=colnames, verbose=parameters["VERBOSE"])
            elif options & GENERICMATCH: matcher=GenericMatch(threshold=dthreshold, colnames=colnames, verbose=parameters["VERBOSE"])
            else: 
                matcher=GenericMatch(threshold=dthreshold, verbose=parameters["VERBOSE"])
                options|=EXPFULL
            print(matcher)
            full= matcher.match( tables, join_type="or", mask=masks )
            av = matcher.finish_matching(full, num_thresh=parameters["NEXP_THRESH"], zpmag=parameters["ZP_MAG"], error_column=error_column)


        output=parameters.get("OUTPUT")
        if output is None or output == '.':
            output=utils.combine_fnames( [ name for name in args] , ntrys=100)
        dname,fname,ext=utils.split_fname(output)

        suffix=""
        if options&EXPFULL: 
            utils.export_table(full,fname="%s/%sfull.fits"%(dname,fname))
            utils.printf("-> %s/%sfull.fits\n"%(dname,fname))
            suffix="match"
        if av:
            utils.export_table(av,"%s/%s%s.fits"%(dname,fname,suffix))
            utils.printf("-> %s/%s%s.fits\n"%(dname,fname,suffix))

        exit_code= scr.EXIT_SUCCESS
    
    elif len(tables)==1: 
        exit_code=scr.EXIT_EARLY
    else:
        utils.perror("No tables loaded for matching.\n")
        exit_code= scr.EXIT_FAIL
    return exit_code

def match_mainentry():
    """StarbugII-match entry"""
    return match_main(sys.argv)
