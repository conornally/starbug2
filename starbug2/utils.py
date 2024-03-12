import time
import os, sys, numpy as np
from parse import parse
import pkg_resources
from astropy.table import Table,hstack,Column,MaskedColumn
from astropy.io import fits
from astropy.wcs import WCS
import starbug2
import requests

printf=sys.stdout.write
perror=sys.stderr.write
puts=lambda s:printf("%s\n"%s)
sbold=lambda s:"\x1b[1m%s\x1b[0m"%s
warn=lambda s:perror("%s%s"%(sbold("Warning: "), s))

def strnktn(s,n,c):
    for _ in range(n): s+=c
    return s

def nputch(n,c): printf(strnktn("",n,c))

def split_fname(fname):
    dname,bname=os.path.split(fname)
    base,ext=os.path.splitext(bname)
    if(not dname): dname='.'
    return dname,base,ext

class loading(object):
    bar=40
    n=0
    N=1
    msg=""

    def __init__(self, N, msg="", res=1):
        self.setlen(N)
        self.msg=msg
        self.startime=time.time()
        self.res=int(res)

    def setlen(self,n): 
        self.N=abs(n)

    def __call__(self):
        self.n+=1
        return self.n<=self.N

    def show(self):
        dec=self.n/self.N
        if (dec==1) or (not self.n%self.res): ## only show once per self.res loads
            out="%s|"%self.msg
            for i in range(self.bar+0):
                out+= ('=' if (i<(self.bar*dec)) else ' ')
            out+="|%.0f%%"%(100*dec)
            
            if self.n:
                etc=(time.time()-self.startime)*(self.N-self.n)/self.n
                nhrs=etc//3600
                nmins= (etc-(nhrs*3600))//60
                nsecs= (etc-(nhrs*3600)-(nmins*60))
                stime=""
                if nhrs: stime+="%dh"%int(nhrs)
                if nmins: stime+="%dm"%int(nmins)
                stime+="%ds"%int(nsecs)
                out+= " ETC:%s"%stime

            printf("\x1b[2K%s\r"%out)
            sys.stdout.flush()
        if(dec==1): printf("\n")

def tabppend(base, tab):
    """
    Is this the same as vstack?
    """
    if(not base): base=tab
    else:
        for line in tab: base.add_row(line)
    return base

def export_region(tab, colour="green", scale_radius=1, region_radius=3, xcol="RA", ycol="DEC", wcs=1, fname="/tmp/out.reg"):
    """
    A handy function to convert the detections in a DS9 region file
    """
    if xcol not in tab.colnames:
        xcols= list(filter(lambda s: 'x'==s[0],tab.colnames))
        if xcols:
            xcol=xcols[0]
            printf("Using '%s' as x position column\n"%sbold(xcol))
            wcs=0

    if ycol not in tab.colnames:
        ycols= list(filter(lambda s: 'y'==s[0],tab.colnames))
        if ycols:
            ycol=ycols[0]
            printf("Using '%s' as y position column\n"%sbold(ycol))
            wcs=0

    if "flux" in tab.colnames and scale_radius: 
        r= (-40.0/np.log10(tab["flux"]))
        r[r<region_radius]=region_radius
        r[np.isnan(r)]=region_radius
    else: r=np.ones(len(tab))*region_radius

    prefix="fk5;" if wcs else ""

    with open(fname, 'w') as fp:
        fp.write("global color=%s width=2\n"%(colour))
        if tab:
            for src, ri in zip(tab,r[r>0]):
                #fp.write("circle %f %f %f;"%(1+src[xcol], 1+src[ycol], ri))
                fp.write("%scircle %f %f %fi\n"%(prefix,src[xcol], src[ycol], ri))
        else:
            perror("unable to open %f\n"%fname)

def parse_unit(raw):
    """
    Take a value with the ability to be cast into several units and parse it
    i.e. 123p -> 123 'pixels'
    """
    recognised={'p':starbug2.PIX, 's':starbug2.ARCSEC, 'm':starbug2.ARCMIN, 'd':starbug2.DEG}
    value=None
    unit=None
    if raw:
        try: 
            value=float(raw)
            unit=None
        except:
            try:
                value=float(raw[:-1])
                unit=recognised.get(raw[-1])
            except:
                perror("unable to parse '%s'\n"%raw)
    return value,unit

def tab2array(tab,colnames=None):
    """
    Returns the contents of the table as a notmal 2D numpy array
    NB: this is different from Table.asarray(), which returns an array of numpy.voids

    if colnames not None, return the subset of the table corresponding to this list
    """
    if not colnames: colnames=tab.colnames
    else: colnames=list( set(colnames)&set(tab.colnames) )
    return np.array( tab[colnames].as_array().tolist() )

def collapse_header(header):
    """
    Convert a dictionary to a Header. 
    Parameters in PARAMFILES have keys longer than 8 chars
    which can cause issues in the fits format. This function turns
    those to comment cards.

    Parameters
    ----------
    header : dict , `fits.Header`
        Header or dictionary to convert to collapse header

    Returns
    -------

    `fits.Header`
        Collapsed Header
    """
    out=fits.Header()
    for key,value in header.items():
        if len(key)>8: 
            out["comment"]=":".join([key,str(value)])
        else: out[key]=value
    return out


def export_table(table, fname=None, header=None):
    """ Export table with correct dtypes """
    dtypes=[]
    table=reindex(table)
    for name in table.colnames:
        if name=="Catalogue_Number": dtypes.append(str)
        elif name=="flag": dtypes.append(np.uint16)
        else: dtypes.append(float)
    table=fill_nan(Table(table,dtype=dtypes))

    if not fname: fname="/tmp/starbug.fits"
    btab=fits.BinTableHDU(data=table, header=header).writeto(fname, overwrite=True, output_verify="fix")

def import_table(fname, verbose=0):
    """
    Slight tweak to `astropy.table.Table.read`. This makes sure that the 
    proper column dtypes are maintained

    Parameters
    ----------
    fname : str
        Path to binary fits table file

    verbose : bool
        Display verbose information
    """
    tab=None
    if os.path.exists(fname):
        if os.path.splitext(fname)[1]==".fits":
            tab=fill_nan(Table.read(fname,format="fits"))
            #tmp=Table()

            #names=[]
            #for index in range(len(tab.dtype.names)):
            #    if tab.dtype[index].kind!='f':
            #        name=tab.colnames[index]
            #        names.append(name)
            #        tmp.add_column( tab[name] )
            #tab.remove_columns( names )

            #tab=tab.filled(np.nan)
            #if tmp: tab=hstack((tab,tmp))

            if not tab.meta.get("FILTER"):
                if (fltr:=find_filter(tab)): 
                    tab.meta["FILTER"]=fltr
            if verbose: printf("-> loaded %s (%s:%d)\n"%(fname,tab.meta.get("FILTER"), len(tab)))
                
        else: perror("Table must fits format\n")
    else: perror("Unable to locate \"%s\"\n"%fname)
    return tab

def fill_nan(table):
    """
    Fill empty values in table with nans
    This is useful for tables that have columns that
    dont support nans (e.g. starbug flag). These will be set to zero instead
    """
    for i,name in enumerate(table.colnames):
        fill_val=np.nan if table.dtype[i].kind=='f' else 0
        if type(table[name])==MaskedColumn: table[name]=table[name].filled(fill_val)

    return table

def find_colnames(tab, basename):
    """
    find substring (basename) within the table colnames
    """
    return [colname for colname in tab.colnames if colname[:len(basename)]==basename]

def combine_fnames(fnames, ntrys=10):
    """
    when matching catalogues, combines the file names into an appropriate combination
    of all the inputs

    Parameters
    ----------
    fnames : list of file names
    
    ntrys : int
        The number of mismatched characters it will allow
    """
    trys=0
    fname=""
    dname,_,ext=split_fname(fnames[0])
    fnames= [split_fname(name)[1] for name in fnames]
    
    for i in range(len(fnames[0])):
        chars= [name[i] for name in fnames if len(name)>i]
        if len(set(chars))==1: fname+=chars[0]
        else: 
            fname+="(%s)"%"".join(sorted(set(chars)))
            trys+=1
        if trys>ntrys: return None
    while ")(" in fname: fname=fname.replace(")(","")
    return "%s/%s%s"%(dname,fname,ext)


def hcascade(tables, colnames=None):
    """
    Similar use as hstack
    Except rather than adding a full new column, the inserted value
    is placed into the leftmost empty column

    Parameters
    ----------
    tables: list of Tables 
        Table to hcascade

    colnames: list of str
        List of column names to include in the stacking.
        If colnames=None, use all possible columns
    """
    tab=fill_nan(hstack(tables))

    if not colnames: colnames=tables[0].colnames
    for name in colnames:
        cols=find_colnames(tab,name)
        if not cols: continue
        move=1
        while move:
            move=0
            for n in range(len(cols)-1,0,-1):
                currmask= np.invert( np.isnan( tab[cols[n]] ) ) ##everything that has a value
                leftmask= np.isnan(tab[cols[n-1]])              ##everything empty in left neighbouring column
                mask=np.logical_and(currmask,leftmask)          ##cur has value and left is empty

                tab[cols[n]]=MaskedColumn(tab[cols[n]])
                tab[cols[n-1]][mask] = tab[cols[n]][mask]
                tab[cols[n]][mask]=tab[cols[n]].info.mask_val
                        #try: tab[cols[n]][mask]*=np.nan
                        #except: tab[cols[n]][mask]*=0
                if sum(mask): move=1
                #tab[cols[n]]=MaskedColumn(np.ma.array(tab[cols[n]], mask=mask), mask=mask)

                #print(tab[cols[n]].info)

        #if name != "flag": ## this is a bodge because it was removing the column if all the stars were good
        #if tab[cols].dtype.kind=='f': # I suspect this could be done with masked bad_vals

            #empty= ( np.nansum( tab2array( tab[cols] ),axis=0 ) ==0)
            #if any(empty): tab.remove_columns(np.array(cols)[empty])

            #for col in cols:
            #    if sum(np.isnan(tab[col]))==len(col):
            #        tab.remove_columns(col)

        cols=find_colnames(tab,name)#[ colname for colname in tab.colnames if name in colname]
        if cols: tab.rename_columns(cols, ["%s_%d"%(name,i+1) for i in range(len(cols))])

    for name in tab.colnames:
        col=tab[name]
        try:
            if col.info.n_bad==col.info.length:
                tab.remove_column(col)
        except: pass
    return tab

def extnames(hdulist): return list(ext.name for ext in hdulist)

def flux2mag(flux,fluxerr=None, zp=1):

    ## sort any type issues in FLUX
    if type(flux)!=np.array: flux=np.array(flux)
    if not flux.shape: flux=np.array([flux])

    # sort type issues in FLUXERR
    if fluxerr is None: fluxerr=np.zeros(len(flux))
    if type(fluxerr)!=np.array: fluxerr=np.array(fluxerr)
    if not fluxerr.shape: fluxerr=np.array([fluxerr])

    mag=np.full( len(flux), np.nan )
    magerr=np.full( len(flux), np.nan )

    maskflux = (flux>0)
    maskferr = (fluxerr>=0)
    mask= np.logical_and( maskflux, maskferr)

    mag[maskflux]= -2.5*np.log10(flux[maskflux]/zp)
    magerr[mask] = 2.5*np.log10( 1.0+( fluxerr[mask]/flux[mask]) )

    return mag,magerr


def flux2ABmag(flux,fluxerr=None):
    return flux2mag( flux, fluxerr, zp=3631.0)


def wget(address, fname=None):
    """
    A really simple "implementation" of wget
    """
    r=requests.get(address)
    if r.status_code==200:
        fname=fname if fname else os.path.basename(address)
        with open(fname,"wb") as fp:
            for chunk in r.iter_content(chunk_size=128):
                fp.write(chunk)
    else: perror("Unable to download \"%s\"\n"%address)

def reindex(table):
    """
    Add indexes into a table
    """
    if "Catalogue_Number" in table.colnames: table.remove_column("Catalogue_Number")
    column=Column(["CN%d"%i for i in range(len(table))], name="Catalogue_Number")
    table.add_column(column,index=0)
    return table

def colour_index(table,keys):
    """
    Allow table indexing with A-B
    """
    out=Table()
    for key in keys:
        if key in table.colnames: out.add_column(table[key])
        elif '-' in key:
            a,b=key.split('-')
            out.add_column(table[a]-table[b],name=key)
    return out

def get_MJysr2Jy_scalefactor(ext):
    """
    Find the unit scale factor to convert an image from MJy/sr to Jy
    """
    scalefactor=1
    if ext.header.get("BUNIT")=="MJy/sr":
        if "PIXAR_SR" in ext.header:
            scalefactor=1e6*float(ext.header["PIXAR_SR"])
    return scalefactor

def find_filter(table):
    """
    Attempt to identify filter for a table from the meta data or column names

    Parameters
    ----------
    table : `astropy.table.Table`
        Table to work on

    Returns
    -------
    Identified filter value, otherwise None
    """
    fltr=None
    if not (fltr:=table.meta.get("FILTER")):
        lst=(set(table.colnames)&set(starbug2.filters.keys()))
        if lst: fltr= lst.pop()
    return fltr

def get_version():
    """
    Try to determine the installed starbug version on the system
    """
    try: version=pkg_resources.get_distribution("starbug2").version 
    except: version="UNKNOWN" ## Github pytest work around for now
    return version

def rmduplicates(seq):
    """
    Take a sequence and rm its duplicates while preserving the order
    of the input

    Parameters
    ----------
    seq : list
        Input list to work on

    Returns
    -------
    A copy of the list with the duplicate elements removed
    """
    seen = set()
    return [x for x in seq if not (x in seen or seen.add(x))]

def cropHDU(hdu, xlim=None, ylim=None):
    """
    Crop an image with multiple extensions

    Parameters:
    -----------
    hdu : fits.HDUList
        A multi frame fits HDUList

    xlim : list
        Pixel X bounds to crop image between

    ylim : list
        Pixel Y bounds to crop image between
    """
    if xlim is None or ylim is None: return None

    for ext in hdu:
        if type(ext) not in (fits.PrimaryHDU, fits.ImageHDU): continue
        if not ext.header["NAXIS"]: continue
        
        ctype=ext.header.get("CTYPE") 
        #if ctype and "-SIP" in ctype:
        ext.header["CTYPE"]="%s-SIP"%ctype
            #ext.header["CTYPE"]=ctype.replace("-SIP","")

        w=WCS(ext.header,relax=False)
        ext.data=ext.data[ xlim[0]:xlim[1], ylim[0]:ylim[1] ]
        ext.header.update( w[ xlim[0]:xlim[1], ylim[0]:ylim[1]].to_header())
    return hdu
        


if __name__ == "__main__":
    print(parse_unit(""))
    print(parse_unit("10p"))
    print(parse_unit("10 p"))
    print(parse_unit("10 D"))
    print(parse_unit("10"))
    print(parse_unit("p10"))


