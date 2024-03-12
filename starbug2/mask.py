import getopt
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import Polygon
from astropy.table import Table
from starbug2.utils import tab2array, colour_index, fill_nan

class Mask(object):
    colour='k'
    def __init__(self, bounds, keys, label=None, **kwargs):
        self.path=Path(bounds)
        if len(keys)==2:
            self.keys=keys
        else: raise Exception
        self.label=label
        
        if "colour" in kwargs: self.colour=kwargs.get("colour")


    def apply(self, dat):
        d=fill_nan(colour_index(dat,self.keys))
        return self.path.contains_points(tab2array(d))
    
    @staticmethod
    def from_file(fname):
        with open(fname) as fp:
            return Mask.from_string(fp.readline())
    @staticmethod
    def from_string(string):
        """
        String Format:
            [-x XCOL] [-y YCOL] [-l Label] : x1 y1 x2 y2 x3 y3 ...
        """
        label=None
        keys=[None,None]
        colour='k'
        _opts,_coords=string.split(':')
        opts,args=getopt.getopt(_opts.split(' '), "c:l:x:y:")
        for opt,optarg in opts:
            if opt=="-x": keys[0]=optarg
            if opt=="-y": keys[1]=optarg
            if opt=="-l": label=optarg.replace('_',' ')
            if opt=="-c": colour=optarg
        coords=_coords.strip().rstrip().split(' ')
        points=np.array(coords, dtype=float).reshape((int(len(coords)/2),2))
        return Mask(points,keys,label=label, colour=colour)

    def plot(self, ax, **kwargs):
        patch=Polygon(self.path._vertices, label=self.label.replace('_',' '), fill=False, edgecolor=self.colour, **kwargs)
        ax.add_patch(patch)
        



if __name__=="__main__":
    s="-yF115W -xF115W-F200W -lTestCut 0 20 1 21 1 24 0 24"
    t=Table.read("/home/conor/sci/proj/ngc6822/paper1/dat/ngc6822.fits",format="fits").filled(np.nan)
    m=Mask.from_string(s)
    mask=m.apply(t)
    import matplotlib.pyplot as plt
    tt=colour_index(t,("F115W-F200W","F115W"))
    plt.scatter(tt["F115W-F200W"], tt["F115W"], c='k', lw=0, s=1)
    #plt.scatter(tt["F115W-F200W"][mask], tt["F115W"][mask], c='r', lw=0, s=1, label=m.label)
    m.plot( plt.gca(), fill=False, edgecolor="blue", label="test")
    plt.legend()
    plt.show()

