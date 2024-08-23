******************
Catalogue Matching
******************
..    include:: <isopub.txt>

*Starbug* offers a selection of generalised matching routines, to accommodate most basic matching requirements. These are accessed through the separate :code:`starbug2-match` executable.
This program takes an arbitrary number of *fits* tables and uses one of three combination methods to match between the RA/DEC columns. 

The basic usage can be displayed with ::

    StarbugII Matching 
    usage: starbug2-match [-BCGfhv] [-e column] [-m mask] [-o output] [-p file.param] [-s KEY=VAL] table.fits ...
        -B  --band               : match in "BAND" mode (does not preserve a column for every frame)
        -C  --cascade            : match in "CASCADE" mode (left justify columns)
        -G  --generic            : match in "GENERIC" mode

        -e  --error   column     : photometric error column ("eflux" or "stdflux")
        -f  --full               : export full catalogue
        -h  --help               : show help message
        -m  --mask    eval       : column evaluation to mask out of matching e.g. -m"~np.isnan(F444W)"
        -o  --output  file.fits  : output matched catalogue
        -p  --param   file.param : load starbug parameter file
        -s  --set     option     : set value in parameter file at runtime (-s MATCH_THRESH=1)
        -v  --verbose            : display verbose outputs

            --band-depr          : match in "old" band mode

        --> typical runs
           $~ starbug2-match -Gfo outfile.fits tab1.fits tab2.fits
           $~ starbug2-match -sMATCH_THRESH=0.2 -sBRIDGE_COL=F444W -Bo out.fits F*W.fits

.. important::
   Matching is a diverse topic, with many intricacies edge cases. :code:`starbug2-match` provides a lot of useful functionality but may not do exactly what the user may require for that specific case. 
   See the soure code documentation for :doc:`Generic Match <../api>` which can be used to write more complex and specialised matching routines if required.


Generic Matching
----------------

The Generic Matching mode, set by :code:`-D` or :code:`--generic` is this simplest combination mode. It will stack the tables horizontally, with every matched source appearing on the same row, and every non-matched source appending to the end. I.e : [*flux_1*,*flux_2*,*flux_3*, ...]

This mode will work through all input tables, building onto the master matched catalogue. This means that a source present in table 2 and 3 but not in 1, will still match to each other, even though they were present in the first.

:code:`starbug2-match` will try to average the appropriate columns by taking either a mean (most columns) or a median (flux/flux_err columns). These will be inserted into the left side of the table. The full result of matching three tables, containing two unique columns will take the form:

+---------+---------+---+---------+---------+---------+---------+---------+---------+
|     A   |     B   |NUM|    A_1  |    B_1  |    A_2  |    B_2  |     A_3 |    B_3  |
+---------+---------+---+---------+---------+---------+---------+---------+---------+
| |check| | |check| | 2 | |check| | |check| | |check| | |check| |         |         |
+---------+---------+---+---------+---------+---------+---------+---------+---------+
| |check| | |check| | 1 | |check| | |check| |         |         |         |         |
+---------+---------+---+---------+---------+---------+---------+---------+---------+
| |check| | |check| | 2 |         |         | |check| | |check| | |check| | |check| |
+---------+---------+---+---------+---------+---------+---------+---------+---------+



Cascade Matching
----------------

Cascade Matching mode, set with :code:`-C` or :code:`--cascade` is generally the same as the Generic mode, however the table gets left justified in between each new catalogue match. The effect of this, is to reduce the overall table size, which may be desired when matching very large numbers of catalogues together. In this mode, the positioning of the columns is considered not important.

A use case may be where the user is matching a large array of (not necessarily overlapping) exposures, covering a large area. In this case, there may be hundreds of tables being combined, which would result in a very large, generally empty final table.

The example table above, would be resolved to take the form :

+---------+---------+---+---------+---------+---------+---------+
|     A   |     B   |NUM|    A_1  |    B_1  |    A_2  |    B_2  |
+---------+---------+---+---------+---------+---------+---------+
| |check| | |check| | 2 | |check| | |check| | |check| | |check| |
+---------+---------+---+---------+---------+---------+---------+
| |check| | |check| | 1 | |check| | |check| |         |         |
+---------+---------+---+---------+---------+---------+---------+
| |check| | |check| | 2 | |check| | |check| | |check| | |check| |
+---------+---------+---+---------+---------+---------+---------+

It is important to note, that the left justification of the final catalogue is a computationally expensive process and slows the execution time.

Band Matching
-------------

Band Matching is used when combining catalogues that are from different photometric bands. In JWST the large spread in FWHM makes is inappropriate to treat the short wavelength data and the long wavelength data the same.

This routine orders the input tables based on increasing PSF FWHM. At every stage, rather than averaging all astrometric positions, it take the shortest wavelength filter possible and uses the position measured in that band. 




Dither Matching
---------------

This is another Generic Match mode that is conducted immediately after the execution of the :code:`starbug2` executable. It is initiated with :code:`-M` or :code:`--match` and simply matches all the final photometric catalogues that were produced in that run. For example ::

    $~ starbug2 -DM image1.fits image2.fits
    --> image1-ap.fits
    --> image2-ap.fits
    --> image(1,2)-apmatch.fits     //This has matched 1 and 2 together

A Typical Run
-------------

In this typical run, we will imagine the scenario where we have four photometric band data: F200W, F444W, F770W, F1000W. That is two NIRCam (one short and one long wavelength) and two MIRI bands. We have conducted photometry on all the individual exposures independently and now wish to match all data into a single catalogue::

    $~ ls
    F200W-expo01-ap.fits  F200W-expo02-ap.fits  F200W-expo03-ap.fits  F200W-expo04-ap.fits
    F444W-expo01-ap.fits  F444W-expo02-ap.fits  F444W-expo03-ap.fits  F444W-expo04-ap.fits
    F770W-expo01-ap.fits  F770W-expo02-ap.fits  F770W-expo03-ap.fits  F770W-expo04-ap.fits
    F1000W-expo01-ap.fits F1000W-expo02-ap.fits F1000W-expo03-ap.fits F1000W-expo04-ap.fits

    // Make a combined for each photometric band individually
    $~ starbug2-match -Go F200W-ap.fits -sMATCH_THRESH=0.1 F200W-expo*.fits 
    $~ starbug2-match -Go F444W-ap.fits -sMATCH_THRESH=0.1 F444W-expo*.fits 
    $~ starbug2-match -Go F770W-ap.fits -sMATCH_THRESH=0.2 F770W-expo*.fits 
    $~ starbug2-match -Go F1000W-ap.fits -sMATCH_THRESH=0.2 F1000W-expo*.fits 

    // Combine all the catalogues together, with an increasing matching threshold
    $~ starbug2-match -Bo final.fits -sMATCH_THRESH=0.1,0.15,0.2 F200W-ap.fits F444W-ap.fits F770W-ap.fits F1000W-ap.fits







Extra Options
*************

Full Catalogue : :code:`-f` or :code:`--full`
    The matching routines will create two different final catalogues:
    - A raw data one, with all the unique columns 
    - A simplified one where all the raw columns have bee averaged.

    The latter simplified version is usually the desired output however if the user needs to keep all the raw columns, set :code:`-f` to output this as well.

Error Column : :code:`-e` or :code:`--error`
    This value sets the column to calculate the error on the magnitude. It is largely deprecated however, the user may wish to set whether to measure an error from "*eflux*" (default) or "*stdflux*" the standard deviation in the flux distribution. 

Mask Logic : :code:`-m` or :code:`--mask`
    **This option is still under development and may change in the future.**
    Add a logic based mask to the matching. This will **exclude** sources from the matching, that resolve "False" to the mask, however will still append them on to the end of the final table. This can be useful for quality controlling the matching so that mismatching is less frequent. An example to remove faint source might be ::

        $~ starbug2-match -G table1.fits table2.fits -m"MAG<24"

    Or a more complex example to match NIRCam and MIRI catalogues together, but ensuring that you only match to NIRCam sources that contain a long wavelength detection, to reduce mismatching chances ::
        
        $~ starbug2-match -G nircam.fits miri.fits -m"~np.isnan(F444W)"


Output : :code:`-o` or :code:`--output`
    :code:`starbug2-match` tries to naively build a combined filename for the matched output, based on the input filenames. This works well if the input filenames are simple but results in horrendous names very quickly. Set the output filename with this option. If paired with :code:`-f` tag, then "*match*" and "*full*" will be appended on the end of this given name.

**MATCH_THRESH**: :code:`-sMATCH_THRESH=`
    Set the maximum separation threshold in arcseconds here. This can be a single value or a comma separated list of multiple values. 
    
    In the latter case, the first pair of tables will take the first value, and the second pair the second value.. and so on. Therefor the list of matching thresholds should be one shorter than the list of input tables. Use this when matching been multiple photometric bands which may have different FWHMs to consider. An example may look like ::

        $~ starbug2-match F115W.fits F200W.fits F356W.fits -sMATCH_THRESH="0.1,0.2"

**MATCH_COLS**: :code:`-sMATCH_COLS=`
    Set any columns to include in the final catalogue here, with a comma separated list. This can be used to drop no longer wanted columns from the inputs.

**NEXP_THRESH**: :code:`-sNEXP_THRESH=`
    Set the minimum number of matches a source must have to be kept in the final catalogue. When matching between multiple exposures of the same area of sky, this can be useful for dropping artefacts like cosmic rays that have been accidentally detected as point sources. By setting a minimum value, a single source with no matches will not persist in the final catalogue
