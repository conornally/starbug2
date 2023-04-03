import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

tab=Table.read("tests/dat/image-ap.fits", format="fits")

sky=SkyCoord(x=tab["xcentroid"],y=tab["ycentroid"],z=np.zeros(len(tab)), representation_type='cartesian')

idx,d2d,d3d=sky.match_to_catalog_3d(sky)
for i in range(len(idx)):
    print(idx[i],d2d[i],d3d[i])

tab.add_column(tab["xcentroid"])
