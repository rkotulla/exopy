#!/usr/bin/env python

import os, sys
import pyfits
from astLib import astWCS
import numpy
import scipy
import scipy.ndimage

from extract_psf import rotate_me

sys.path.append("/work/podi_devel")
from podi_definitions import *


def subtract_psf(psf, fn, hdu, ra, dec, superscale=1., normalize_zp=30.):

    cutout_size = 700
    
    # hdu = pyfits.open(fn)
    data = hdu[1].data
    
    wcs = astWCS.WCS(hdu[1].header, mode='pyfits')
    xy = wcs.wcs2pix(ra, dec)
    print xy
    # resulting pixels are 0-based (not 1, as usual FITS coords)

    skylevel = hdu[0].header['SKYLEVEL']

    #cutout_size = 700

    x1 = int(xy[0] - cutout_size)
    x2 = int(xy[0] + cutout_size)
    y1 = int(xy[1] - cutout_size)
    y2 = int(xy[1] + cutout_size)
    print x1,x2,y1,y2

    cutout = data[y1:y2, x1:x2]

    #
    # Now rotate the PSF image and insert it into the full frame
    #

    # rotator angle
    rotator_angle = hdu[0].header['ROTSTART']
    if (rotator_angle > 180.): rotator_angle -= 360.
    if (rotator_angle < -180.): rotator_angle += 360.
    elevation = numpy.degrees(float(hdu[0].header['ELMAP']))
    print rotator_angle, -1*rotator_angle, elevation

    rotate_angle = (rotator_angle - elevation)

    print "running geometric_transform"
    output_shape = (3*cutout_size, 3*cutout_size)
    setup = {
            'in_center_x': psf.shape[1]/2., 
            'in_center_y': psf.shape[0]/2.,
            'sin_rot': numpy.sin(numpy.radians(rotate_angle)), 
            'cos_rot': numpy.cos(numpy.radians(rotate_angle)),
            'out_center_x': cutout.shape[1]/2., 
            'out_center_y': cutout.shape[0]/2., 
            'scale': 1./superscale,
            }

    out = scipy.ndimage.interpolation.geometric_transform(
        input=psf,
        mapping=rotate_me,
        output_shape=cutout.shape,
        order=1,
        mode='constant',
        cval=numpy.NaN,
        prefilter=False,
        extra_keywords=setup
        )
    print "done!"

    # pyfits.PrimaryHDU(data=out).writeto(fn[:-5]+".psf2sub.fits", clobber=True)
    # out_bpm = scipy.ndimage.interpolation.geometric_transform(
    #     input=bpm,
    #     mapping=rotate_me,
    #     output_shape=output_shape,
    #     order=1,
    #     mode='constant',
    #     cval=numpy.NaN,
    #     prefilter=False,
    #     extra_keywords=setup
    #     )
    print "done!"
    # out[out_bpm > 0.1] = numpy.NaN


    # # make bad pixel masks, and grow a bit
    # cutout[cutout > 50000] = numpy.NaN
    
    # bpm = numpy.zeros(cutout.shape)
    # bpm[numpy.isnan(cutout)] = 1.
    
    # bpm_conv = scipy.ndimage.filters.convolve(
    #     input=bpm, 
    #     weights=numpy.ones((5,5)), 
    #     output=None, 
    #     mode='reflect', 
    #     cval=0.0)
    # bpm = bpm_conv

    # cutout[numpy.isnan(cutout)] = 0.
    # cutout -= skylevel


    # rotated = scipy.ndimage.interpolation.rotate(
    #     input=cutout,
    #     angle=-1*rotator_angle,
    #     order=3,
    #     )
    # bpm_rotated = scipy.ndimage.interpolation.rotate(
    #     input=bpm,
    #     angle=-1*rotator_angle,
    #     order=3,
    #     )


    # rotated[bpm_rotated > 0.1] = numpy.NaN

    # pyfits.PrimaryHDU(data=cutout, header=hdu[0].header).writeto("psf.fits", clobber=True)
    # pyfits.PrimaryHDU(data=rotated, header=hdu[0].header).writeto("psf_rotated.fits", clobber=True)

    # pyfits.PrimaryHDU(data=bpm_rotated, header=hdu[0].header).writeto("psf_rotated_bpm.fits", clobber=True)

    # pyfits.PrimaryHDU(data=out, header=hdu[0].header).writeto("psf_geom.fits", clobber=True)

    #
    # Apply the zeropoint scaling
    #
    magzero = hdu[0].header['PHOTZP_X']
    scale = numpy.power(10., -0.4*(magzero - 30.))
    print magzero, scale
    print "#####", fn, magzero, normalize_zp, scale

    #out = out / scale

    #cutout -= out
    
    #data[:,:] = 0.
    data[y1:y2, x1:x2] -= (out / scale)
    #cutout = out

    img_hdu = pyfits.ImageHDU(data=data, header=hdu[1].header)

    return img_hdu


if __name__ == "__main__":

    psf_fn = sys.argv[1]
    psf_hdu = pyfits.open(psf_fn)
    psf = psf_hdu[1].data
    print psf.shape

    ra = float(sys.argv[2])
    dec = float(sys.argv[3])
    
    fn = sys.argv[4]
    hdu = pyfits.open(fn)

    try:
        out_fn = sys.argv[5]
    except:
        out_fn = fn[:-5]+".psfsub.fits"

    superscale = 1.0
    psfsub_hdu = subtract_psf(psf, fn, hdu, ra, dec)


    phdu = pyfits.PrimaryHDU(header=hdu[0].header)
    phdu.header['SPRSCALE'] = superscale

    hdulist = pyfits.HDUList([phdu, psfsub_hdu])
    clobberfile(out_fn)
    hdulist.writeto(out_fn, clobber=True)
    
