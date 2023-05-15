#import numpy as np
#from astropy.io import fits
#from astropy.convolution import convolve, RickerWavelet2DKernel, Gaussian2DKernel
#from astropy.visualization import ZScaleInterval
#import matplotlib.pyplot as plt
#from scipy.signal import ricker
#from scipy.ndimage import convolve
#from skimage.feature import match_template
#from starbug2.routines import Detection_Routine
#import os
#import webbpsf
#
#instr=webbpsf.NIRCam()
#instr.filter="F444W"
#instr.calc_psf().writeto("/tmp/psf.fits",overwrite=True)
##os.system("ds9 -multiframe /tmp/test.fits")
import time, os

print(os.system("starbug2 -h"))



#import multiprocessing
#
#def fnwait(t):
#    print("%d waiting:%d"%(os.getpid(),t))
#    time.sleep(t)
#    print("done waiting:%d"%t)
#
#    return 0
#
#with multiprocessing.Pool(3) as p:
#    p.map(fnwait,(10,11,12))
#    print('oh?')






#dat=fits.open("/home/conor/dat/NGC346/JWST/forced/F444W/jw01227002001_02105_00001_nrcalong_destrip_jhat.fits")
##dat=fits.open("/home/conor/dat/NGC346/JWST/stage3/F770W/jw01227-c1002_t005_miri_f770w_i2d.fits")
##plt.imshow(ZScaleInterval()(dat[1].data))
##plt.show()
#
#kernel = Gaussian2DKernel(x_stddev=2)
#kernel = RickerWavelet2DKernel(1)
#
#mask=[np.where( dat[1].data==0)]
#image=dat[1].data
#mask=np.where(np.isnan(image))
#image[mask]*=0
#conv=convolve(image, kernel)
#corr=match_template(conv/np.amax(conv), kernel.array)
#offset=kernel.array.shape[0]//2
#
#
#
#tab=Detection_Routine(sig_src=3, sharplo=0.1, sharphi=1, verbose=1).find_stars(corr)
#
#fig,(ax1,ax2, ax3)=plt.subplots(1,3, figsize=(30,10), sharex=True, sharey=True)
#
#ax1.imshow(ZScaleInterval()(image),interpolation=None)
#ax2.imshow(ZScaleInterval()(conv),interpolation=None)
#ax3.imshow(ZScaleInterval()(corr),interpolation=None)
#ax2.scatter( tab["xcentroid"]+offset, tab["ycentroid"]+offset, s=80, facecolors='none', edgecolors='r')
#plt.show()



#def RickerWavelet (data, thres_val, mask):
#    #Create a RickerWavelet2dKernel ('Mexican hat') and convolve the data with it.
#    kernel = RickerWavelet2DKernel(1)
#    conv = convolve(data, kernel)
#
#    #Calculate the sigma clipped statistics and use source detection using the RickerWavelet2DKernel
#    mean, median, std = sigma_clipped_stats(conv, sigma=3.0, mask=mask)
#    print ('mean, median and std of convolved image =', mean, median, std)
#    threshold = median + (thres_val * std)
#    tbl_conv = find_peaks(conv, threshold, box_size=11)
#    positions_conv = np.transpose((tbl_conv['x_peak'], tbl_conv['y_peak']))
#    conv_apertures = CircularAperture(positions_conv, r=2.0)
#
#    #Match the data with the template of the RickerWavelet 2D kernel and use this for source detection.
#    corr = match_template(conv/np.amax(conv), kernel.array)
#
#    #conv and corr have a slight size difference due to the match templating.
#    #The mask must be adjusted accordingly.
#    ny, nx = conv.shape #Size of conv image
#    y, x = corr.shape #Size of corr image
#    size_dif_y = (ny - y)//2
#    size_dif_x = (nx - x)//2
#
#    #Mask for corr image
#    mask_corr = mask[size_dif_y:ny-size_dif_y, size_dif_x:nx-size_dif_x]
#
#    offset = kernel.array.shape[0] // 2
#    mean, median, std = sigma_clipped_stats(corr, sigma=3.0 , mask=mask_corr)
#    print ('mean, median and std of "corr" image =', mean, median, std)
#    threshold = median + (thres_val * std) #Arbitrary threshold
#    tbl_corr = find_peaks(corr, threshold, box_size=11)
#    positions_corr = np.transpose((tbl_corr['x_peak']+offset, tbl_corr['y_peak']+offset))
#    corr_apertures = CircularAperture(positions_corr, r=2.0)
#
#    return conv, corr, positions_corr, corr_apertures, positions_conv, conv_apertures
#
#
#
#
