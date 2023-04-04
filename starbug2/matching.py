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

def generic_match(catalogues, threshold=0.25, add_src=True):
    """
    """
    threshold=threshold*u.arcsec
    base=Table(None)

    for n,cat in enumerate(catalogues,1):
        if "Catalogue_Number" in cat.colnames: cat.remove_column("Catalogue_Number")
        if not len(base):
            tmp=cat.copy()
        else:
            tmp=Table(tmp,dtype=[float]*len(tmp.colnames)).filled(np.nan) ## fill empty values with null
            idx,d2d,_=_match(base,cat)
            tmp=Table(np.full( (len(base),len(cat.colnames)), np.nan), names=cat.colnames)

            for src,IDX,sep in zip(cat,idx,d2d):
                if (sep<=threshold) and (sep==min(d2d[idx==IDX])): ## GOOD MATCH
                    tmp[IDX]=src
                elif add_src:   ##BAD MATCH / NEW SOURCE
                    tmp.add_row( src )
        tmp.rename_columns( tmp.colnames, ["%s_%d"%(name,n) for name in tmp.colnames] )
        base=hstack((base,tmp))
        #base=Table(base,dtype=[float]*len(base.colnames)).filled(np.nan) ## fill empty values with null
    return base



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

def cascade_match(catalogues, threshold, colnames):
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
    base=Table(base,dtype=[float]*len(base.colnames)).filled(np.nan)
    return finish_matching(base, colnames)

def band_match(catalogues, threshold, colnames):
    """
    Given a list of catalogues (with filter names in the meta data), match them
    in order of decreasing astrometric accuracy. 
    If F115W, F444W, F770W are input, the F115W centroid positions will be 
    taken as "correct". If a source is not resolved in this band, the next most 
    astrometrically accurate position is taken, i.e. F444W
    """
    threshold=threshold*u.arcsec
    colnames= list( name for name in catalogues[0].colnames if name in colnames)
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

def stage_match(stage2, stage3, threshold):
    """
    Match together a stage 2 and stage 3 catalogue
    if the star is resolved in just the stage 3, it takes on that value
    if the star is resolved in both, the stage3 is ignored
    """
    
    if stage2.meta.get("CALIBLEVEL")!=2: perror("WARNING: stage2 catalogue CALIBLEVEL=%s does not match\n"%stage2.meta.get("CALIB_LEVEL"))
    if stage3.meta.get("CALIBLEVEL")!=3: perror("WARNING: stage3 catalogue CALIBLEVEL=%s does not match\n"%stage2.meta.get("CALIB_LEVEL"))
    #av,tab=dither_match((stage2,stage3),threshold,stage2.colnames)
    idx,d2d,_=_match(stage2,stage3)
    
    tmp=Table(np.full((len(stage2),len(stage3.colnames)), np.nan), names=stage3.colnames)
    for src,IDX,sep in zip(stage3,idx,d2d):
        if (sep<threshold) and (sep==min(d2d[idx==IDX])):
            tmp[IDX]=src
        else: tmp.add_row(src)

    mask=np.isfinite(tmp["RA"])
    tmp.rename_columns(tmp.colnames, list("%s_stage3"%name for name in tmp.colnames))
    tab=hstack( (stage2.copy(),tmp) )
    tab.add_column(mask, name="MASK")

    print(tab)
    return tab 



def finish_matching(tab, colnames):
    """
    Averaging all the values. Combining source flags and building a NUM column
    """
    flags=np.full(len(tab),starbug2.SRC_GOOD, dtype=np.uint16)
    av=Table(np.full((len(tab),len(colnames)),np.nan), names=colnames)

    for name in colnames:
        all_cols=find_colnames(tab,name)
        #if not (all_cols:=find_colnames(tab,name)): continue
        if not all_cols: continue
        col=Column(None, name=name)
        ar=tab2array(tab, colnames=all_cols)
        if name=="flux":
            col=Column(np.nanmean(ar,axis=1), name=name)
            median=np.nanmedian(ar,axis=1)
            if "stdflux" not in colnames: av.add_column(Column(np.nanstd(ar,axis=1),name="stdflux")) 
            ## if median and mean are >5% different, flag as SRC_VAR
            flags[ np.abs(median-col)>(col/5.0)] |= starbug2.SRC_VAR
        elif name== "eflux":
            col=Column(np.sqrt(np.nansum(ar*ar, axis=1)), name=name)
        elif name=="stdflux": 
            col=Column(np.nanmax(ar,axis=1),name=name)
        elif name=="flag":
            col=Column(flags, name=name)
            for fcol in ar.T: col|=fcol.astype(np.uint16)
        else: col=Column(np.nanmedian(ar, axis=1),name=name)
        
        av[name]=col
    mag,magerr=flux2ABmag(av["flux"],av["eflux"], tab.meta["FILTER"])
    av.add_column(mag,name=tab.meta["FILTER"])
    av.add_column(magerr,name="e%s"%tab.meta["FILTER"])
    narr= np.nansum( np.invert( np.isnan(tab2array(tab,find_colnames(tab,colnames[0])))),axis=1)
    av.add_column(Column(narr, name="NUM"))
    return (av,tab)


def remove_NUM(tab, N):
    """
    Remove sources from the list of tab if they have >N non matches
    """
    if "NUM" in tab.colnames:
        pass#mask=tab["NUM"]


if __name__=="__main__":
    stage2=Table.read("/dat/1zw18/jwst/stage2/1zw18-apstage2.fits")
    stage3=Table.read("/dat/1zw18/jwst/stage3/1Zw18-apstage3.fits")


    tab=stage_match(stage2,stage3, 0.25*u.arcsec)
    export_table(tab,fname="/tmp/tab.fits")
