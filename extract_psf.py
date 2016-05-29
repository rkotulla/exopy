#!/usr/bin/env python

import os, sys
import pyfits
from astLib import astWCS
import numpy
import scipy
import scipy.ndimage

sys.path.append("/work/podi_devel/")
from podi_definitions import *

def rotate_me(output_coords, 
              in_center_x, in_center_y,
              sin_rot, cos_rot,
              out_center_x, out_center_y, scale,
):
    
    # in_center_x = extra[0]
    # in_center_y = extra[1]
    # rot = extra[2]

    # out_center_x = extra[3]
    # out_center_y = extra[4]
    # scale = extra[5]

    rel_x = (output_coords[0] - out_center_x) / scale
    rel_y = (output_coords[1] - out_center_y) / scale

    in_x = cos_rot*rel_x - sin_rot*rel_y + in_center_x
    in_y = sin_rot*rel_x + cos_rot*rel_y + in_center_y
    return (in_x, in_y)



def extract_psf(fn, ra, dec,
                cutout_size=700,
                superscale=2.,
                normalize_zp=30.):

    hdu = pyfits.open(fn)
    hdu.writeto("dump2.fits", clobber=True)

    data = hdu[1].data

    #pyfits.PrimaryHDU(data=data).writeto("dump.fits", clobber=True)
    #os._exit(0)

    #return data

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
    #return cutout

    # make bad pixel masks, and grow a bit
    cutout[cutout > 50000] = numpy.NaN
    #return cutout

    bpm = numpy.zeros(cutout.shape)
    bpm[numpy.isnan(cutout)] = 1.
    
    bpm_conv = scipy.ndimage.filters.convolve(
        input=bpm, 
        weights=numpy.ones((5,5)), 
        output=None, 
        mode='reflect', 
        cval=0.0)
    bpm = bpm_conv

    #return bpm

    cutout[numpy.isnan(cutout)] = 0.
    cutout -= skylevel

    # rotator angle
    rotator_angle = hdu[0].header['ROTSTART']
    if (rotator_angle > 180.): rotator_angle -= 360.
    if (rotator_angle < -180.): rotator_angle += 360.
    elevation = numpy.degrees(float(hdu[0].header['ELMAP']))
    print rotator_angle, -1*rotator_angle, elevation

    derotate_angle = -1. * (rotator_angle - elevation)

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

    print "running geometric_transform"
    output_shape = (3*cutout_size, 3*cutout_size)
    setup = {
            'in_center_x': cutout.shape[1]/2., 
            'in_center_y': cutout.shape[0]/2.,
            'sin_rot': numpy.sin(numpy.radians(derotate_angle)), 
            'cos_rot': numpy.cos(numpy.radians(derotate_angle)),
            'out_center_x': output_shape[1]/2., 
            'out_center_y': output_shape[0]/2., 
            'scale': superscale,
            }
    out = scipy.ndimage.interpolation.geometric_transform(
        input=cutout,
        mapping=rotate_me,
        output_shape=output_shape,
        order=1,
        mode='constant',
        cval=numpy.NaN,
        prefilter=False,
        extra_keywords=setup
        )
    print "done!"
    out_bpm = scipy.ndimage.interpolation.geometric_transform(
        input=bpm,
        mapping=rotate_me,
        output_shape=output_shape,
        order=1,
        mode='constant',
        cval=numpy.NaN,
        prefilter=False,
        extra_keywords=setup
        )
    print "done!"
    out[out_bpm > 0.1] = numpy.NaN

    # rotated[bpm_rotated > 0.1] = numpy.NaN


    # pyfits.PrimaryHDU(data=cutout, header=hdu[0].header).writeto("psf.fits", clobber=True)
    # pyfits.PrimaryHDU(data=rotated, header=hdu[0].header).writeto("psf_rotated.fits", clobber=True)

    # pyfits.PrimaryHDU(data=bpm_rotated, header=hdu[0].header).writeto("psf_rotated_bpm.fits", clobber=True)

    # pyfits.PrimaryHDU(data=out, header=hdu[0].header).writeto("psf_geom.fits", clobber=True)

    #
    # Figure out an appropriate scaling factor to normalize PSF extractions
    #
    magzero = hdu[0].header['PHOTZP_X']
    scale = numpy.power(10., -0.4*(magzero - normalize_zp))
    out *= scale

    print "#####", fn, magzero, normalize_zp, scale
    return out

if __name__ == "__main__":
    
    ra = float(sys.argv[1])
    dec = float(sys.argv[2])
    
    fn = sys.argv[3]
    
    try:
        out_fn = sys.argv[4]
    except:
        out_fn = fn[:-5]+".psf.fits"

    superscale = 1.0
    psf = extract_psf(fn, ra, dec, superscale=1.)

    phdu = pyfits.PrimaryHDU(data=psf)
    phdu.header['SPRSCALE'] = superscale

    clobberfile(out_fn)
    phdu.writeto(out_fn, clobber=True)
    
