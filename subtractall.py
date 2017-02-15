#!/usr/bin/env python

import os, sys

import extract_psf
import subtract_psf
import makeprofile

sys.path.append("/work/podi_devel/")
from podi_definitions import *
import podi_imcombine


if __name__ == "__main__":

    task = sys.argv[1]

    superscale = 1.0
    magnitude = 0.


    if (task == "create"):

        ra = float(sys.argv[2])
        dec = float(sys.argv[3])

        shapes = []
        files = sys.argv[4:]
        shape_list = {}
        for fn in files:

            print "\n"*3, "extracting PSF from %s" % (fn)
            out_fn = fn[:-5]+".psf.fits"
            psf = extract_psf.extract_psf(fn, ra, dec, superscale=1., magnitude=magnitude)

            phdu = pyfits.PrimaryHDU(data=psf)
            phdu.header['SPRSCALE'] = superscale

            clobberfile(out_fn)
            phdu.writeto(out_fn, clobber=True)

            #
            # Compute profile
            #
            good_profile, scale_2d = makeprofile.compute_profile(data=psf, fn=fn)

            #
            # Divide PSF by profile to create shape image
            #
            shape_img = psf / scale_2d
            shapes.append(shape_img)
            #shape_list[fn] = shape_img

            out_fn = fn[:-5]+".shape.fits"
            clobberfile(out_fn)
            phdu = pyfits.PrimaryHDU(data=shape_img)
            phdu.writeto(out_fn)


    elif (task == "combine"):
        #
        # Now we have all the shape images, run imcombine
        #
        shape_combined = podi_imcombine.imcombine_data(shapes, operation="sigmaclipmean")
        mastershape_fn  = "%s_mastershape.fits" % (object_name)
        pyfits.PrimaryHDU(data=shape_combined).writeto(mastershape_fn, clobber=True)

    elif (task == "subtract"):

        mastershape_fn = sys.argv[2]
        ra = float(sys.argv[3])
        dec = float(sys.argv[4])

        mastershape_hdu = pyfits.open(mastershape_fn)
        mastershape = mastershape_hdu[1].data

        files = sys.argv[5:]

        # 
        # now scale the psf shape image back up, and subtract from image
        #
        for fn in files:

            psf = extract_psf.extract_psf(fn, ra, dec, superscale=1., 
                                          magnitude=magnitude, rotate=False)
            good_profile, scale_2d = makeprofile.compute_profile(data=psf, fn=fn)


            psf = mastershape * scale_2d

            hdu = pyfits.open(fn)
            psfsub_hdu = subtract_psf.subtract_psf(psf, fn, hdu, ra, dec)

            out_fn = fn[:-5]+".psfsub.fits"
            phdu = pyfits.PrimaryHDU(header=hdu[0].header)
            hdulist = pyfits.HDUList([phdu, psfsub_hdu])
            clobberfile(out_fn)
            hdulist.writeto(out_fn, clobber=True)
