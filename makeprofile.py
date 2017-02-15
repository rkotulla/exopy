#!/usr/bin/env python

import os, sys, numpy, pyfits
import scipy, scipy.interpolate

def compute_profile(data, fn):

    iy,ix = numpy.indices(data.shape, dtype=numpy.float)
    iy -= data.shape[0]/2.
    ix -= data.shape[1]/2.

    r = numpy.hypot(ix,iy)

    valid = numpy.isfinite(data)

    profile = numpy.append(r[valid].reshape((-1,1)), 
                           data[valid].reshape((-1,1)),
                           axis=1)
    numpy.savetxt(fn[:-5]+".profile", profile[::5])

    #
    # compute median in radial bins
    #
    max_width = data.shape[1]/3
    n_bins = 100
    bin_width = max_width / n_bins

    profile_median = numpy.zeros((n_bins,2))
    profile_median[:,:] = numpy.NaN

    for ibin in range(n_bins):
        min_r = ibin*bin_width
        max_r = min_r + bin_width

        in_this_bin = (profile[:,0] > min_r) & (profile[:,0] <= max_r)
        if (numpy.sum(in_this_bin) <= 0):
            continue

        median = numpy.median(profile[in_this_bin][:,1])

        profile_median[ibin] = [ 0.5*(max_r+min_r), median]

    numpy.savetxt(fn[:-5]+".profile2", profile_median)

    good_profile = profile_median[numpy.isfinite(profile_median[:,0])]

    #
    # Now compute an interpolation function
    #
    interpol = scipy.interpolate.interp1d(
        good_profile[:,0],
        good_profile[:,1],
        bounds_error=False)

    scale_2d = interpol(r)
    pyfits.PrimaryHDU(data=scale_2d).writeto(fn[:-5]+".scale.fits", clobber=True)

    return good_profile, scale_2d

if __name__ == "__main__":

    for fn in sys.argv[1:]:
        hdu = pyfits.open(fn)

        for ext in hdu:
            if (ext.header['NAXIS'] == 2):
                data = ext.data
                break
            
            good_profile, scale_2d = compute_profile(data, fn)
            
