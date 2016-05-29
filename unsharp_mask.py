#!/usr/bin/env python

import pyfits, os, sys, numpy
import scipy, scipy.ndimage


if __name__ == "__main__":

    for infile in sys.argv[1:]:

        hdulist = pyfits.open(infile)

        # for now, assume it's a swarped file
        for ext in hdulist:
            data = ext.data
            if (type(data) != numpy.ndarray):
                continue
            if (data.ndim != 2):
                continue

            print infile, ext.name, ext.data.shape

            # smoothed = scipy.ndimage.filters.gaussian_filter(
            #     input=data, 
            #     sigma=10, 
            #     order=0, output=None, mode='constant', cval=0.0
            #     )
            # leftover = data - smoothed

            median = scipy.ndimage.filters.median_filter(
                input=data, 
                size=21, 
                output=None, mode='constant', cval=0.0
                )
            leftover = data - median

            ext.data = leftover


        hdulist.writeto(infile[:-5]+".unsharped.fits", clobber=True)



        
