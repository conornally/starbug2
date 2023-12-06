"""
Starbug matching functions
Primarily this is the main routines for dither/band/generic matching which are at the core
of starbug2 and starbug2-matc
"""
import os
import numpy as np
import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Column, Table, hstack, vstack

import starbug2
from starbug2.utils import *
from starbug2.param import load_params

class Matcher(object):
    """
    Base matching class

    Parameters
    ----------
    catalogues :

    colnames : list
        List of str column names to include in the matching.
        Everything else will be discarded

    fltr : str
        Specifically set the filter of the catalogues 

    threshold : float
        Separation threshold in arcseconds

    pfile : str
        Parameter filename

    """
    method="Generic Matching"
    def __init__(self, threshold=None, colnames=None, fltr=None, verbose=None, pfile=None):
        options=load_params(pfile)
        self.threshold  =options.get("MATCH_THRESH")
        self.filter     =options.get("FILTER")
        self.verbose    =options.get("VERBOSE")

        if threshold is not None: self.threshold=threshold
        self.threshold *= u.arcsec

        if fltr is not None: self.filter=fltr
        if verbose is not None: self.verbose=verbose

        self.colnames=colnames 
        self.load=loading(1)
        
    def log(self,msg):
        if self.verbose: printf(msg)

    def __str__(self):
        s=[ "%s:"%self.method,
            "Filter: %s"%self.filter,
            "Colnames: %s"%self.colnames,
            "Threshold: %s\""%self.threshold]
        return "\n".join(s)

    def __call__(self, *args, **kwargs):
        return self.match(*args, **kwargs)

    def init_catalogues(self, catalogues):
        """
        """
        ## Must copy here maybe?
        if len(catalogues)>=2:
            self.load=loading( sum( len(cat) for cat in catalogues[1:]), msg="matching")

        if self.colnames is None:
            self.colnames=[]
            for cat in catalogues:
                self.colnames+=cat.colnames
        self.colnames = rmduplicates(self.colnames)
        if "Catalogue_Number" in self.colnames: self.colnames.remove("Catalogue_Number")

        for n,catalogue in enumerate(catalogues):
            keep=set(catalogue.colnames)&set(self.colnames)
            keep=sorted( keep, key= lambda s:self.colnames.index(s))
            catalogues[n]=catalogue[keep]
            self.colnames=keep  # This maybe wants to go somewhere else but it ensures that colnames doesnt contain anything not in any tables

        if not self.filter:
            if (fltr:=catalogues[0].meta.get("FILTER")) is None:
                fltr="MAG"
            self.filter=fltr

        return catalogues


    def match(self, catalogues, join_type="or", **kwargs):
        """
        This matching works as a basic match. Everything is included and the column
        names have _N appended to the end. 
        Parameters
        ----------
        join_type : str
            Joing method
            "or" include sources in any catalogue
            "and" only include sources in all catalogues
        
        Returns
        -------
        .
        """
        catalogues=self.init_catalogues(catalogues)
        base=Table(None)
        if join_type=="and": perror("join_type 'and' not fully implemented\n")

        for n,cat in enumerate(catalogues,1):
            if not len(base):
                tmp=cat.copy()
            else:
                idx,d2d,_=self._match(base,cat)
                tmp=Table(np.full( (len(base),len(cat.colnames)), np.nan), names=cat.colnames, dtype=cat[self.colnames].dtype)

                for src,IDX,sep in zip(cat,idx,d2d):
                    self.load()
                    if self.verbose: self.load.show()

                    if (sep<=self.threshold) and (sep==min(d2d[idx==IDX])): ## GOOD MATCH
                        tmp[IDX]=src
                    elif join_type=="or":##BAD MATCH / NEW SOURCE
                        tmp.add_row( src )

            tmp.rename_columns( tmp.colnames, ["%s_%d"%(name,n) for name in tmp.colnames] )
            base=fill_nan(hstack((base,tmp)))
        return base

    @staticmethod
    def _match(cat1, cat2):
        """
        Base matching function between two catalogues
        
        Parameters
        ----------
        cat1 and cat2 : `astropy.table.Table`
            two astropy tables containing columns with "RA/DEC" in the column names.
            This could be RA_1, RA_2 .. 
            If several columns are located, they will be nanmeaned together

        Returns
        -------
            idx,d2d,d3d : the same as SkyCoord.match_to_catalog_3d
        """
        cat1=fill_nan(cat1)
        cat2=fill_nan(cat2)

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


    def finish_matching(self, tab, error_column="eflux", num_thresh=-1, discard_outliers=False, zpmag=0, colnames=None):
        """
        Averaging all the values. Combining source flags and building a NUM column

        Parameters
        ----------
        tab : `astropy.table.Table`
            Table to work on

        error_column : str
            Column containing resultant photometric errors to be used to calculate the magnitude error
            "eflux" - use the eflux column (the normal photometric error column)
            "stdflux"   - use "stdflux" column as error on flux

        num_thresh : int
            Minimum number of matches a source must have.
            NUM values smaller than this will be removed from the table.
            If num_thresh<=0, no cropping will happen

        discard_outliers : bool
            Choose whether to remove outling values from the averaging. 
            This may be usful if some flux values are wildly different from others in the set

        zpmag : float
            Zero point (Magnitude) to be applied to the magnitude after it is calculated

        colnames : list
            List of colnames to include in the averaging. If None, use self.colnames instead
        """
        flags=np.full(len(tab),starbug2.SRC_GOOD, dtype=np.uint16)
        av=Table(np.full((len(tab),len(self.colnames)),np.nan), names=self.colnames)
        
        if colnames is None: colnames=self.colnames 
        for name in colnames:
            if (all_cols:=find_colnames(tab,name)):
                col=Column(None, name=name)
                ar=tab2array(tab, colnames=all_cols)
                if name=="flux":
                    col=Column(np.nanmedian(ar,axis=1), name=name)
                    mean=np.nanmean(ar,axis=1)
                    if "stdflux" not in self.colnames: av.add_column(Column(np.nanstd(ar,axis=1),name="stdflux")) 
                    ## if median and mean are >5% different, flag as SRC_VAR
                    flags[ np.abs(mean-col)>(col/5.0)] |= starbug2.SRC_VAR
                elif name== "eflux":
                    col=Column(np.sqrt(np.nansum(ar*ar, axis=1)), name=name)
                elif name=="stdflux": 
                    col=Column(np.nanmedian(ar,axis=1),name=name)
                elif name=="flag":
                    col=Column(flags, name=name)
                    for fcol in ar.T: flags|=fcol.astype(np.uint16)
                elif name=="NUM":
                    col=Column(np.nansum(ar, axis=1), name=name)
                else: col=Column(np.nanmedian(ar, axis=1),name=name)
                
                av[name]=col

        av["flag"]=Column(flags,name="flag")
        if "flux" in av.colnames:
            ecol=av[error_column] if error_column in av.colnames else None
            mag,magerr=flux2mag(av["flux"], fluxerr=ecol)
            mag+=zpmag

            if self.filter in av.colnames: av.remove_column(self.filter)
            if "e%s"%self.filter in av.colnames: av.remove_column("e%s"%self.filter)
            av.add_column(mag,name=self.filter)
            av.add_column(magerr,name="e%s"%self.filter)
            #perror("There was no zero point added here!\n")

        if "NUM" not in av.colnames:
            narr= np.nansum( np.invert( np.isnan(tab2array(tab,find_colnames(tab,self.colnames[0])))),axis=1)
            av.add_column(Column(narr, name="NUM"))

            if num_thresh>0:
                av.remove_rows( av["NUM"]<num_thresh)
        return av
    
    def build_meta(self,catalogues):
        pass
        

class CascadeMatch(Matcher):
    method="Cascade Matching"
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def match(self, catalogues, **kwargs):
        """
        match a list of catalogues with RA and DEC columns
        INPUT: 
            catalogues : list( astropy.table.Tables) catalogues to match
            threshold  : max separation (arcsec) for two soures to be considered the same source
            colnames   : column names to include in matches
        RETURN:
            a left aligned catalogue of all the matched values
        """
        catalogues=self.init_catalogues(catalogues)
        base=Table( None, meta=catalogues[0].meta)

        for n,cat in enumerate(catalogues,1):
            if n==1:
                tmp=cat.copy()
            else:
                idx,d2d,_=self._match(base,cat)
                #tmp=Table(np.full((len(base),ncol),np.nan), names=colnames)

                ## If the tmp table is larger than cat, then I can copy in the unmatched
                ## sources without add_row
                drow=len(base) 
                mask=(d2d>self.threshold)
                tmp=Table(np.full((len(base)+sum(mask),len(self.colnames)),np.nan), names=self.colnames, dtype=cat[self.colnames].dtype)

                for src,IDX,sep in zip(cat,idx,d2d):
                    self.load()
                    if self.verbose: self.load.show()

                    if (sep<=self.threshold) and (sep==min(d2d[idx==IDX])): ##It does match
                        tmp[IDX]=src
                    else:   ##APPEND
                        if drow<len(tmp): ##This is a time saving idea im trying
                            tmp[drow]=src
                            drow+=1 
                        else: 
                            tmp.add_row(src[self.colnames]) ##i can purely use add_row to simplifiy the code
            tmp.rename_columns( tmp.colnames, ["%s_%d"%(name,n) for name in tmp.colnames] )
            base=hcascade((base,tmp), colnames=self.colnames)
        base=fill_nan(base)
        return base

class DitherMatch(Matcher):
    def __init__(self, catalogues, pfile=None):
        super(self,DitherMatch).__init__(catalogues, pfile)

    def match(self, **kwargs):
        return None

class BandMatch(Matcher):
    method="Band Matching"
    def __init__(self, **kwargs):

        if "fltr" in kwargs:
            if not isinstance(kwargs["fltr"], list):
                warn("fltr input should be a list, there may be unexpected behaviour\n")

        if "threshold" in kwargs:
            if isinstance(kwargs["threshold"],list):
                kwargs["threshold"]=np.array(kwargs["threshold"])

        super().__init__(**kwargs)

    def order_catalogues(self,catalogues):
        """
        Reorder catalogue list into increasing wavelength size
        This only works for JWST bands. Unrecognised filters will be left unchanged.
        The function should also set the self.filter variable if possible

        Parameter
        ---------
        catalogues : list
            List of `astropy.table.Table` with meta keys "FILTER"

        Returns
        -------
        The same list reordered
        """

        status=-1
        _ii=None
        sorters = [ lambda t: list(starbug2.filters.keys()).index( t.meta.get("FILTER")),   ## META in JWST filters
                    lambda t: list(starbug2.filters.keys()).index( (set(t.colnames)&set(starbug2.filters.keys())).pop()), ## colnames in JWST filters
                    lambda t: self.filter.index( t.meta.get("FILTER")),                     ## META in self.filters
                    lambda t: self.filter.index( (set(t.colnames)&set(self.filter)).pop() ) ## colnames in JWST filters
                    ]

        for n,fn in enumerate(sorters):
            try:
                catalogues.sort(key=fn)
                _ii=map(fn,catalogues)
                status=n
                break
            except Exception as e:
                pass

        if status<0:
            perror("Unable to reorder catalogues, leaving input order untouched.\n")
        elif status<=1 and (_ii is not None): ## JWST filters
            self.filter=[list(starbug2.filters.keys())[i] for i in _ii]

        self.load=loading(sum(len(c) for c in catalogues[1:]))

        return catalogues

    def jwst_order(self,catalogues):
        pass

    def match(self, catalogues, method="first", **kwargs):
        """
        Given a list of catalogues, it will reorder them into increasing wavelength
        or to match the fltr= keyword in the initialiser. 
        The matching then uses the shortest wavelegnth availables position.
        I.e If F115W, F444W, F770W are input, the F115W centroid positions will be 
        taken as "correct". If a source is not resolved in this band, the next most 
        astrometrically accurate position is taken, i.e. F444W

        Parameters
        ----------
        catalogues : list
            List of `astropy.table.Table` objects containing the meta item "FILTER=XXX"

        method: str
            Centroid method:
            -   "first" :   Use the position corresponding to the earliest appearance of the source
            -   "last"  :   Use the position corresponding to the latest appearance of the source
            -   "bootsrap": ..
            -   "average" : ..

        Returns
        -------
        """
        catalogues=self.order_catalogues(catalogues)

        if isinstance(self.filter,list) and len(self.filter)==len(catalogues):
            printf("Bands: %s\n"%', '.join(self.filter))
        else: printf("Bands: Unknown\n")

        if type(self.threshold.value) in (list,np.ndarray):# and len(self.threshold)==(len(catalogues)-1):
            if len(self.threshold)!=(len(catalogues)-1):
                warn("Threshold values must be scalar or list with length 1 less than the catalogue list. The final element is being ignored.\n")
                self.threshold = self.threshold[:-1]
        else: 
            self.threshold=np.full(len(catalogues)-1, self.threshold)*u.arcsec
        printf("Thresholds: %s\n"%", ".join(["%g\""%g for g in self.threshold.value]))

        if self.colnames is None: self.colnames=["RA","DEC", "flag", "NUM", *self.filter, *["e%s"%f for f in self.filter]]
        printf("Columns: %s\n"%", ".join(self.colnames))

        if method not in ("first","last","bootstrap"): method="first"

        #########
        # Begin #
        #########
        
        base=Table(None)
        for n,tab in enumerate(catalogues):
            colnames= [ name for name in self.colnames if name in tab.colnames]
            self.load.msg="%s (%g\")"%(self.filter[n], self.threshold[n-1].value)

            if not len(base): 
                tmp=tab[colnames].copy()
            else:
                idx,d2d,_=self._match(base,tab)
                tmp=Table(np.full( (len(base),len(colnames)), np.nan), names=colnames, dtype=tab[colnames].dtype)

                for ii,(src,IDX,sep) in enumerate(zip(tab,idx,d2d)):
                    self.load();self.load.show()
                    if (sep<=self.threshold[n-1]) and (sep==min(d2d[idx==IDX])):
                        tmp[IDX]=src[colnames]
                    else:
                        tmp.add_row(src[colnames])

            colnames.remove("RA")
            colnames.remove("DEC")
            base=fill_nan(hstack((base, tmp[colnames])))
            base.rename_columns(colnames, ["%s_%d"%(name,n+1) for name in colnames])

            if "RA" not in base.colnames: base=fill_nan( hstack((tmp["RA","DEC"], base)) )
            elif method=="first":
                _mask=np.logical_and( np.isnan(base["RA"]), tmp["RA"]!=np.nan)
                base["RA"][_mask]=tmp["RA"][_mask]
                base["DEC"][_mask]=tmp["DEC"][_mask]
            elif method=="last":
                _mask= ~np.isnan(tmp["RA"])
                base["RA"][_mask]=tmp["RA"][_mask]
                base["DEC"][_mask]=tmp["DEC"][_mask]
            elif method=="bootstrap":
                _mask= ~np.isnan(tmp["RA"])
                base.rename_columns(("RA","DEC"),("_RA_%d"%n, "_DEC_%d"%n))
                base=hstack( (base,tmp[["RA","DEC"]]) )

        ####################
        # Fix column names #
        for name in self.colnames:
            all_cols=find_colnames(base,name)
            if len(all_cols)==1:
                base.rename_column(all_cols.pop(), name)

        ################################
        # Finalise NUM and flag column #
        tmp=self.finish_matching( base, colnames=["NUM", "flag"] )
        base.remove_columns( (*find_colnames(base,"NUM"), *find_colnames(base,"flag")) )
        base.add_column(tmp["NUM"], index=2)
        base.add_column(tmp["flag"],index=3)

        return base


def band_match(catalogues, colnames=("RA","DEC")):
    """
    Given a list of catalogues (with filter names in the meta data), match them
    in order of decreasing astrometric accuracy. 
    If F115W, F444W, F770W are input, the F115W centroid positions will be 
    taken as "correct". If a source is not resolved in this band, the next most 
    astrometrically accurate position is taken, i.e. F444W
    """
    #threshold=threshold*u.arcsec
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
        elif (_tmp:=set(starbug2.filters.keys()) & set(tab.colnames)):
            ii=list(starbug2.filters.keys()).index(_tmp.pop())
            tables[ii]=tab
            mask[ii]=True
        else: perror("Cannot find 'FILTER' in table meta (skipping)..\n")
    s="Bands: "
    for fltr,tab in zip(starbug2.filters.keys(),tables):
        if tab: s+="%5s "%fltr
        #else: s+=". "
    puts(s)

    ### Match in increasing wavelength order
    base=Table(None)
    load=loading(sum( [len(t) for t in tables[mask][1:]]),"matching", res=100)
    for fltr,tab in zip(starbug2.filters.keys(),tables):
        if not tab: continue
        tab.remove_rows( np.isnan(tab[fltr]) ) ## removing empty magnitude rows
        load.msg="matching:%s"%fltr
        _colnames= list( name for name in tab.colnames if name in colnames)
        if not len(base): 
            tmp=tab[_colnames].copy()
        else:
            idx,d2d,_=Matcher._match(base,tab)
            tmp=Table(np.full( (len(base),len(_colnames)), np.nan), names=_colnames)
            
            ###################################
            # Hard coding separations for now #
            separation=0.06
            f_id=list(starbug2.filters.keys()).index(fltr)
            if f_id >= list(starbug2.filters.keys()).index("F277W"): separation=0.10
            if f_id >= list(starbug2.filters.keys()).index("F560W"): separation=0.15
            if f_id >= list(starbug2.filters.keys()).index("F1000W"): separation=0.20
            if f_id >= list(starbug2.filters.keys()).index("F1500W"): separation=0.25


            for ii,(src,IDX,sep) in enumerate(zip(tab,idx,d2d)):
                load.msg="matching:%s(%.2g\")"%(fltr,separation)
                load();load.show()
                if (sep<=separation*u.arcsec) and (sep==min(d2d[idx==IDX])):
                    for name in _colnames: tmp[IDX][name]=src[name]
                else:
                    tmp.add_row(src[_colnames])
            #tmp=generic_match((base,tab), threshold=threshold, add_src=True, load=load)

        #base=hstack(( base,tmp[["flux","eflux"]] ))
        #mag,magerr=flux2ABmag(tmp["flux"], tmp["eflux"],fltr)
        #base.add_column(mag,name=fltr)
        #base.add_column(magerr,name="e%s"%fltr)

        #base.rename_column("flux","%s_flux"%fltr)
        #base.rename_column("eflux","%s_eflux"%fltr)

        tmp.rename_column("flag","flag_%s"%fltr)
        base=hstack(( base,tmp[[fltr,"e%s"%fltr,"flag_%s"%fltr]] ))#.filled(np.nan)
        base=Table(base,dtype=[float]*len(base.colnames)).filled(np.nan)

        ### Only keep the most astromectrically correct position
        if "RA" not in base.colnames: base=hstack(( tmp[["RA","DEC"]], base))
        else:
            _mask=np.logical_and( np.isnan(base["RA"]), tmp["RA"]!=np.nan)
            base["RA"][_mask]=tmp["RA"][_mask]
            base["DEC"][_mask]=tmp["DEC"][_mask]

    ## Sort out flags
    flag=np.zeros(len(base),dtype=np.uint16)
    for fcol in find_colnames(base,"flag"):
        flag|=base[fcol].value.astype(np.uint16)
        base.remove_column(fcol)
    base.add_column(flag,name="flag")
    
    return base.filled(np.nan)






####################################
# Disgarded functions and what not #
####################################


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

def generic_match():pass
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

def bootstrap_match(catalogues):
    """
    A more generalised band matching routine
    A - B
    B - C
    C - D
    """
    print("THIS ISNT READY")

    output=Table()

    THRESHHOLD=0.3
    filters= [c.meta.get("FILTER") for c in catalogues]
    base=None

    if len( catalogues ) >1:
        base=catalogues[0]["RA","DEC","%s"%filters[0],"e%s"%filters[0]]

        for n,cat in enumerate(catalogues[1:],1):
            f=filters[n]
            cat=cat["RA","DEC","%s"%f,"e%s"%f]
            base=generic_match( (base,cat), add_src=True, average=False)

            _f=filters[n-1]   #fix base colnames. Lower filter columns get finalised
            base.rename_columns(("RA_1","DEC_1","%s_1"%_f,"e%s_1"%_f), 
                                ("RA_%s"%_f,"DEC_%s"%_f,"%s"%_f,"e%s"%_f) ) 
            for cn in base.colnames:
                if cn[-2:]=="_1": base.rename_column(cn,cn[:-2])

            #Get larger wavelength ready for next match
            base.rename_columns( ("RA_2","DEC_2","%s_2"%f,"e%s_2"%f),
                                 ("RA","DEC","%s"%f,"e%s"%f) )

            print(n,len(catalogues))
            if n==len(catalogues)-1: ##Final match
                base.rename_columns(("RA","DEC"),("RA_%s"%f,"DEC_%s"%f) )

    else:
        perror("Must include more than one catalogue.\n")
    return base


