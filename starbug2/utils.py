import time
import os, sys, numpy as np
from parse import parse
from scipy.optimize import minimize
from astropy.table import Table,hstack
from astropy.io import fits
import starbug2

printf=sys.stdout.write
perror=sys.stderr.write
puts=lambda s:printf("%s\n"%s)

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
        self.res=res

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
    if(not base): base=tab
    else:
        for line in tab: base.add_row(line)
    return base

def export_region(tab, fname="/tmp/out.reg"):
    """
    A handy function to convert the detections in a DS( region file
    """
    xcol="xcentroid" if "xcentroid" in tab.colnames else "x_0"
    ycol="ycentroid" if "ycentroid" in tab.colnames else "y_0"
    xcol="RA"
    ycol="DEC"
    
    if "peak" in tab.colnames: r= np.sqrt(tab["peak"]**0.7)*2.5 if "peak" in tab.colnames else 5
    else: r=np.ones(len(tab))*10

    with open(fname, 'w') as fp:
        fp.write("global color=cyan\n")
        if tab:
            for src, ri in zip(tab,r[r>0]):
                #fp.write("circle %f %f %f;"%(1+src[xcol], 1+src[ycol], ri))
                fp.write("fk5;circle %f %f %fi\n"%(src[xcol], src[ycol], ri))
        else:
            perror("unable to open %f\n"%fname)

def load_table(fname):
    table=None
    if os.path.exists(fname):
        fp=fits.open(fname)
        table=Table(fp[1].data._get_raw_data())
        fp.close()
    else: perror("unable to open %f\n"%fname)
    return table

def load_params(fname):
    """
    Convert a parameter file into a dictionary of options
    INPUT:  fname=path/to/file.param
    RETURN: dictionary of options
    """
    config={}
    if os.path.exists(fname):
        with open(fname, "r") as fp:
            for line in fp.readlines():
                if line[0] in "# \t\n":
                    continue
                key,value,_=parse("{}={}//{}\n",line)
                key=key.strip()
                value=value.strip()
                try:
                    value=float(value)
                except:
                    pass

                ## Special case values
                if key=="PSFDIR": value=os.path.expandvars(value)
                config[key]=value
    else:
        perror("config file \"%s\" does not exist\n"%fname)
        config=None
    return config

def tab2array(tab,colnames=None):
    """
    Returns the contents of the table as a notmal 2D numpy array
    NB: this is different from Table.asarray(), which returns an array of numpy.voids

    if colnames not None, return the subset of the table corresponding to this list
    """
    if not colnames: colnames=tab.colnames
    else: colnames=list( set(colnames)&set(tab.colnames) )
    return np.array( tab[colnames].as_array().tolist() )


def export_table(table, fname=None, header=None):
    if table:
        if not fname: fname="/tmp/starbug.fits"
        fits.BinTableHDU(data=table, header=header).writeto(fname, overwrite=True)

def find_colnames(tab, basename):
    """
    find substring (basename) within the table colnames
    """
    return [colname for colname in tab.colnames if colname[:len(basename)]==basename]

def combine_fnames(fnames, ntrys=10):
    """
    when matching catalogues, combines the file names into an appropriate combination
    of all the inputs
    INPUT:  fnames=list of file names
            ntrys=the number of mismatched characters it will allow
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
    return "%s/%s%s"%(dname,fname,ext)


def hcascade(tables, colnames=None):
    """
    Similar use as hstack
    Except rather than adding a full new column, the inserted value
    is placed into the leftmost empty column
    INPUT:
        tables: list of Tables to hstack
        colnames: colnames to use
    """
    tab=hstack(tables)
    tab=Table(tab, dtype=[float]*len(tab.colnames)).filled(np.nan)

    if not colnames: colnames=tables[0].colnames
    for name in colnames:
        cols=find_colnames(tab,name)
        move=1
        while move:
            move=0
            for n in range(len(cols)-1,0,-1):
                currmask= np.invert( np.isnan( tab[cols[n]] ) ) ##everything that has a value
                leftmask= np.isnan(tab[cols[n-1]])              ##everything empty in left neighbouring column
                mask=np.logical_and(currmask,leftmask)          ##cur has value and left is empty

                tab[cols[n-1]][mask] = tab[cols[n]][mask]
                tab[cols[n]][mask]*=np.nan
                if sum(mask): move=1

        empty= ( np.nansum( tab2array( tab[cols] ),axis=0 ) ==0)
        if any(empty): tab.remove_columns(np.array(cols)[empty])

        cols=find_colnames(tab,name)#[ colname for colname in tab.colnames if name in colname]
        tab.rename_columns(cols, ["%s_%d"%(name,i+1) for i in range(len(cols))])
    return tab

def extnames(hdulist): return list(ext.name for ext in hdulist)

def flux2ABmag(flux,fluxerr=None, filter=None):
    zp=3631.0
    if filter and filter in list(starbug2.ZP.keys()): zp=starbug2.ZP[filter][0]
    mag=-2.5*np.log10(flux/zp)
    magerr=2.5*np.log10(1.0+(fluxerr/flux)) if any(fluxerr) else None
    return(mag,magerr)


if __name__ == "__main__":
    import glob
    fnames=glob.glob("/home/conor/dat/NGC346/JWST/stage2-destriped/de-striped_F335M/de-striped_jw01227002001_02101_0000*_nrc*long_cal.fits")
    fname=combine_fnames(fnames)
    print(fname)
