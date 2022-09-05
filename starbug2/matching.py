"""
Starbug matching functions
Primarily this is the main routines for dither/band/generic matching which are at the core
of starbug2 and starbug2-matc
"""
import numpy as np
import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Column, Table, hstack

import starbug2
from starbug2.utils import *

def exp_info(hdulist):
    """
    Get the exposure information about a hdulist 
    INPUT:  HDUList or ImageHDU or BinTableHDU
    RETURN: dictionary of relevant information:
            >   EXPOSURE, DETECTOR, FILTER
    """
    info={  "FILTER":None,
            "OBSERVTN":0,
            "VISIT":0,
            "EXPOSURE":0,
            "DETECTOR":None
            }

    if type(hdulist) in (fits.ImageHDU, fits.BinTableHDU):
        hdulist=fits.HDUList(hdulist)

    for hdu in hdulist:
        for key in info:
            if key in hdu.header: info[key]=hdu.header[key]
    return info

def sort_exposures(catalogues):
    """
    Given a list of catalogue files, this will return the fitsHDULists as a series of
    nested dictionaries sorted by:
    >   BAND
    >   OBSERVATION ID
    >   VISIT ID
    >   DETECTOR                -- These two have been switched
    >   DITHER (EXPOSURE)       -- These two have been switched
    """
    out={}
    modules=["NRCA1","NRCA2","NRCA3","NRCA4","NRCB1","NRCB2","NRCB3","NRCB4","NRCALONG","NRCBLONG", "UNKNOWN"]
    for cat in catalogues:
        info=exp_info(cat)
        #print(info, cat[0].header["FILENAME"])
        
        if info["FILTER"] not in out.keys(): 
            out[info["FILTER"]]={}

        if info["OBSERVTN"] not in out[info["FILTER"]].keys(): 
            out[info["FILTER"]][info["OBSERVTN"]]={}

        if info["VISIT"] not in out[info["FILTER"]][info["OBSERVTN"]].keys():
            out[info["FILTER"]][info["OBSERVTN"]][info["VISIT"]]={}

        if info["DETECTOR"] not in out[info["FILTER"]][info["OBSERVTN"]][info["VISIT"]].keys():
            out[info["FILTER"]][info["OBSERVTN"]][info["VISIT"]][info["DETECTOR"]]=[]

        #if info["EXPOSURE"] not in out[info["FILTER"]][info["OBSERVTN"]][info["VISIT"]].keys():
            #out[info["FILTER"]][info["OBSERVTN"]][info["VISIT"]][info["EXPOSURE"]]= np.full(len(modules), None)

        out[info["FILTER"]][info["OBSERVTN"]][info["VISIT"]][info["DETECTOR"]].append(cat)
        #index= modules.index(info["DETECTOR"]) if info["DETECTOR"] in modules else -1
        #out[info["FILTER"]][info["OBSERVTN"]][info["VISIT"]][info["EXPOSURE"]][index] = cat
    return out

def _match(cat1, cat2):
    """
    if not len(cat1) and len(cat2): 
        n=len(cat2)
        return np.array((range(n), sys.maxsize*np.ones(n)*u.degree, sys.maxsize*np.ones(n)))
    if not len(cat2) and len(cat1): 
        n=len(cat1)
        return np.array((range(n), sys.maxsize*np.ones(n)*u.degree, sys.maxsize*np.ones(n)))
    """
    _ra_cols= list( name for name in cat1.colnames if "RA" in name)
    _dec_cols= list( name for name in cat1.colnames if "DEC" in name)
    _ra= np.nanmean( tab2array( cat1, colnames=_ra_cols), axis=1) # this still breaks if the source isnt matched in all the columns, the 999 will increase the average
    _dec=np.nanmean( tab2array( cat1, colnames=_dec_cols), axis=1)
    skycoord1=SkyCoord( ra=_ra*u.deg, dec=_dec*u.deg)

    _ra_cols= list( name for name in cat2.colnames if "RA" in name)
    _dec_cols= list( name for name in cat2.colnames if "DEC" in name)
    _ra= np.nanmean( tab2array( cat2, colnames=_ra_cols), axis=1)
    _dec=np.nanmean( tab2array( cat2, colnames=_dec_cols), axis=1)
    skycoord2=SkyCoord( ra=_ra*u.deg, dec=_dec*u.deg)

    return skycoord2.match_to_catalog_3d(skycoord1)

def dither_match(catalogues, threshold, colnames):
    """
    This is the match for when you simultaneously detect sources over
    a series of pipeline stage2 dithers and wich to combine it into a single catalogue
    INPUT:  catalogues: a list of astropy tables
            threshold:  (float) maximum separation between two sources to match
            colnames:   names to include in the output catalogue
    RETURNS: a combined catalogue with paired sources appearing on the same line
    """
    threshold=threshold*u.arcsec
    colnames= list(name for name in colnames if name in catalogues[0].colnames)#list( set(catalogues[0].colnames) & set(colnames) )
    base=Table( None)#, names=colnames )

    for n,cat in enumerate(catalogues,1):
        if not len(base):
            tmp=cat[colnames].copy()
        else:
            idx,d2d,_=_match(base,cat)
            tmp=Table(np.full((len(base),len(colnames)),np.nan), names=colnames)

            for src,IDX,sep in zip(cat,idx,d2d):
                if (sep<=threshold) and (sep==min(d2d[idx==IDX])): ## GOOD MATCH
                    for name in colnames: tmp[IDX][name]=src[name]
                else:   ##BAD MATCH / NEW SOURCE
                    tmp.add_row( src[colnames] )

        tmp.rename_columns( colnames, list("%s_%d"%(name,n) for name in colnames))
        base=hstack((base,tmp))
        base=Table(base,dtype=[float]*len(base.colnames)).filled(np.nan)
    return finish_matching(base, colnames)

def generic_match(catalogues, threshold, colnames):
    """
    match a list of catalogues with RA and DEC columns
    INPUT: 
        catalogues : list( astropy.table.Tables) catalogues to match
        threshold  : max separation (arcsec) for two soures to be considered the same source
        colnames   : column names to include in matches
    RETURN:
        a left aligned catalogue of all the matched values
    """
    threshold*=u.arcsec
    colnames= list( name for name in catalogues[0].colnames if name in colnames)
    ncol=len(colnames)
    base=Table( None, meta=catalogues[0].meta)#, names=colnames )
    load=loading(sum([len(cat) for cat in catalogues[1:]]),"matching")

    for n,cat in enumerate(catalogues,1):
        if n==1:
            tmp=cat[colnames].copy()
        else:
            idx,d2d,_=_match(base,cat)
            #tmp=Table(np.full((len(base),ncol),np.nan), names=colnames)

            ## If the tmp table is larger than cat, then I can copy in the unmatched
            ## sources without add_row
            drow=len(base) 
            mask=(d2d>threshold)
            tmp=Table(np.full((len(base)+sum(mask),ncol),np.nan), names=colnames)

            for src,IDX,sep in zip(cat,idx,d2d):
                load();load.show()
                if (sep<=threshold) and (sep==min(d2d[idx==IDX])): ##It does match
                    for name in colnames: tmp[IDX][name]=src[name]
                else:   ##APPEND
                    if drow<len(tmp): ##This is a time saving idea im trying
                        tmp[drow]=src[colnames]
                        drow+=1 
                    else: 
                        tmp.add_row(src[colnames]) ##i can purely use add_row to simplifiy the code
        tmp.rename_columns(colnames, list("%s_%d"%(name,n) for name in colnames))
        base=hcascade((base,tmp), colnames=colnames)
    return finish_matching(base, colnames)

def band_match(catalogues, threshold, colnames):
    ### ORDER the tables into the correct order (increasing wavelength)
    tables=np.full( len(starbug2.filters), None)
    mask=np.full(len(starbug2.filters), False)
    for tab in catalogues:
        if "FILTER" in tab.meta.keys():
            if tab.meta["FILTER"] in starbug2.filters: 
                ii=list(starbug2.filters.keys()).index(tab.meta["FILTER"])
                tables[ii]=tab
                mask[ii]=True
            else: perror("Unknown filter '%s' (skipping)..\n"%tab.meta["FILTER"])
        else: perror("Cannot find 'FILTER' in table meta (skipping)..\n")
    s="Bands: "
    for fltr,tab in zip(starbug2.filters.keys(),tables):
        if tab: s+="%5s "%fltr
        #else: s+=". "
    puts(s)
    colnames=["RA","DEC","flux","eflux"]
    threshold=threshold*u.arcsec

    ### Match in increasing wavelength order
    base=Table(None)
    load=loading(sum( [len(t) for t in tables[mask][1:]]),"matching", res=10)
    for fltr,tab in zip(starbug2.filters.keys(),tables):
        if not tab: continue
        load.msg="matching:%s"%fltr
        if not len(base): tmp=tab[colnames].copy()
        else:
            idx,d2d,_=_match(base,tab)
            tmp=Table(np.full( (len(base),len(colnames)), np.nan), names=colnames)

            for ii,(src,IDX,sep) in enumerate(zip(tab,idx,d2d)):
                load();load.show()
                if (sep<=threshold) and (sep==min(d2d[idx==IDX])):
                    for name in colnames: tmp[IDX][name]=src[name]
                else:
                    tmp.add_row(src[colnames])

        base=hstack(( base,tmp[["flux","eflux"]] ))
        mag,magerr=flux2ABmag(tmp["flux"], tmp["eflux"],fltr)
        base.add_column(mag,name=fltr)
        base.add_column(magerr,name="e%s"%fltr)

        base.rename_column("flux","%s_flux"%fltr)
        base.rename_column("eflux","%s_eflux"%fltr)
        base=Table(base,dtype=[float]*len(base.colnames)).filled(np.nan)

        ### Only keep the most astromectrically correct position
        if "RA" not in base.colnames: base=hstack(( tmp[["RA","DEC"]], base))
        else:
            _mask=np.logical_and( np.isnan(base["RA"]), tmp["RA"]!=np.nan)
            base["RA"][_mask]=tmp["RA"][_mask]
            base["DEC"][_mask]=tmp["DEC"][_mask]
    export_table(base,"/tmp/tab.fits")
    return base










def finish_matching(tab, colnames):
    """
    Averaging all the values. Combining source flags and building a NUM column
    """
    flags=np.full(len(tab),starbug2.SRC_GOOD, dtype=np.uint32)
    av=Table(np.full((len(tab),len(colnames)),np.nan), names=colnames)
    #if "ap_stdflux" in tab.colnames: 
        #tab.remove_column("ap_stdflux")
        #if "ap_stdflux" in colnames: colnames.remove("ap_stdflux")
        #perror("WARNING, REMOVING DITHER FLUX STDEV COL\n")


    for name in colnames:
        if not (all_cols:=find_colnames(tab,name)): continue
        col=Column(None, name=name)
        ar=tab2array(tab, colnames=all_cols)
        if name=="flux":
            col=Column(np.nanmean(ar,axis=1), name=name)
            median=np.nanmedian(ar,axis=1)
            if "stdflux" not in colnames: av.add_column(Column(np.nanstd(ar,axis=1),name="stdflux")) 
            av.add_column(Column( -2.5*np.log10(col/starbug2.ZP[tab.meta["FILTER"]][0]), name="%s_mag"%tab.meta["FILTER"]))

            ## if median and mean are >5% different, flag as SRC_VAR
            flags[ np.abs(median-col)>(col/5.0)] |= starbug2.SRC_VAR
            
        elif name== "eflux":
            col=Column(np.nansum(ar*ar, axis=1), name=name)
            av.add_column(Column( 2.5*np.log10(1+(av["stdflux"]/av["flux"])), name="%s_emag"%tab.meta["FILTER"] ))
        elif name=="stdflux": 
            col=Column(np.nanmax(ar,axis=1),name=name)

        elif name=="flag":
            col=Column(flags, name=name)
            for fcol in ar.T: col|=fcol.astype(np.uint32)
        else: col=Column(np.nanmedian(ar, axis=1),name=name)
        
        av[name]=col
    narr= np.nansum( np.invert( np.isnan(tab2array(tab,find_colnames(tab,colnames[0])))),axis=1)
    av.add_column(Column(narr, name="NUM"))
    #return hstack((av,tab))
    return (av,tab)

if __name__=="__main__":
    #fnames=(#"/home/conor/dat/NGC346/JWST/stage2-destriped/de-striped_F335M/out.fits",
    fnames=("/home/conor/dat/NGC346/JWST/stage2-destriped/de-striped_F200W/out.fits",
            "/home/conor/dat/NGC346/JWST/stage2-destriped/de-striped_F444W/out.fits",
            )
            #"/home/conor/dat/NGC346/JWST/stage2-destriped/de-striped_F187N/out.fits")
    tables=[Table().read(fname) for fname in fnames]
    band_match(tables,0.5,None)

    """
    fnames=sorted(glob.glob("/home/conor/dat/NGC346/JWST/stage2-destriped/de-striped_F335M/*-apmatch.fits"))
    colnames=("RA","DEC","ap_flux")
    tables=list(Table.read(fname, format="fits") for fname in fnames)
    export_table(generic_match(tables),"out.fits")

    fnames=glob.glob("/home/conor/dat/NGC346/JWST/stage2-destriped/de-striped_F187N/*-ap.fits")
    fps=[ fits.open(fp) for fp in fnames]
    out=sort_exposures(fps)

    print("band ob visit exp det a1 a2 a3 a4 b1 b2 b3 b4 al bl ??")
    for band,obs in out.items():
        for ob,visits in obs.items():
            for visit,exps in visits.items():
                for exp,dets in exps.items():
                    #for det,fp in dets.items():
                    s=":     "
                    for item in dets: s+= "X  " if item else ".. "

                    ar=np.array([Table(fp[1].data._get_raw_data()) if fp else None for fp in dets])
                    out=match_detectormodule( ar )
                    print(band,ob,visit,exp, s)
    """
