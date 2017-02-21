#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
import math


if __name__ == "__main__":

    fn = sys.argv[1]

    cx = float(sys.argv[2])
    cy = float(sys.argv[3])

    seeing_fwhm = float(sys.argv[4])
    mag = float(sys.argv[5])
    out_dir = sys.argv[6]

    if (not os.path.isdir(out_dir)):
        os.makedirs(out_dir)

    n_total = 1000.
    n_per_image = 100.
    max_position_offset = 300

    print fn, cx, cy
    print "seeing in arcsec", seeing_fwhm

    hdulist = pyfits.open(fn)
    data = hdulist[0].data

    n_chunks = int(math.ceil(n_total / float(n_per_image)))

    n_sigma = 5
    pixelscale = hdulist[0].header['PIXLSCAL']
    seeing_sigma_pixels = seeing_fwhm / pixelscale / 2.355
    print "sigma in pixels", seeing_sigma_pixels

    star_box_size = 2*n_sigma*seeing_fwhm/pixelscale
    starmodel = numpy.zeros((star_box_size, star_box_size))

    sy, sx = numpy.indices(starmodel.shape, dtype=numpy.float)
    sy -= star_box_size/2.
    sx -= star_box_size/2.
    sr = numpy.hypot(sx,sy)

    starmodel = numpy.exp( -sr**2 / (2*seeing_sigma_pixels**2))
    total_flux = numpy.sum(starmodel)
    zeropoint = hdulist[0].header['MAGZERO']
    total_mag = -2.5* math.log10(total_flux) + zeropoint
    print "mag/zp:",total_mag, zeropoint

    delta_mag = total_mag - mag
    print "mag offset:", delta_mag
    scaling = numpy.power(10., 0.4*(delta_mag))
    print scaling
    starmodel *= scaling

    print "calibrated mag:", -2.5*numpy.log10(numpy.sum(starmodel)) + zeropoint
    pyfits.PrimaryHDU(data=starmodel).writeto("starmodel.fits", clobber=True)

    gain = hdulist[0].header['GAIN']
    skylevel = 0.#hdulist[0].header['SKYLEVEL']
    starmodel_photons = (starmodel + skylevel) * gain
    starmodel_sigma = numpy.sqrt(starmodel_photons)

    #
    # Now we have a valid star-model, add it to a number of frames
    #
    stars = [pyfits.PrimaryHDU()]
    for frame_number in range(n_chunks):

        _outimg = data.copy()

        center_pos = numpy.round(
            (numpy.random.random((n_per_image, 2)) - 0.5) * 2 * \
                     max_position_offset + [cx,cy], 0).astype(numpy.int)
        r = numpy.hypot(center_pos[:,0]+1.-cx, center_pos[:,1]+1.-cy)

        # write a catalog of positions so we can go back and check which
        # stars we are able to recover.
        cat_fn = "%s/mag%.1f__chunk%02d.cat" % (
            out_dir, mag, frame_number+1
        )
        # when writing, add 1 to x & y to account for 1-indexing in FITS
        numpy.savetxt(
            cat_fn,
            numpy.append(center_pos,r.reshape((-1,1)),axis=1)+[1.,1.,0],
            "%.2f"
        )

        for i_star, starpos in enumerate(center_pos):
            # simulate a noisy version of the starmodel
            star_noise = numpy.random.normal(
                loc=0, scale=starmodel_sigma, size=starmodel_sigma.shape
            ) / gain
            star = starmodel + star_noise
            stars.append(pyfits.ImageHDU(data=star))

            # xy1 = starpos - 0.5*star.shape
            # xy2 = xy1 + star.shape

            x1 = starpos[0]-int(0.5*star.shape[1])
            y1 = starpos[1]-int(0.5*star.shape[0])
            x2 = x1 + star.shape[1]
            y2 = y1 + star.shape[0]

            _outimg[int(y1):int(y2),
                    int(x1):int(x2)] += star


        out_fn = "%s/mag%.1f__chunk%02d.fits" % (
            out_dir, mag, frame_number+1
        )
        hdulist[0].data = _outimg
        hdulist[0].header['_NINSERT'] = n_per_image
        hdulist[0].header['_STARMAG'] = mag
        hdulist[0].header['_CENTERX'] = cx
        hdulist[0].header['_CENTERY'] = cy
        hdulist.writeto(out_fn, clobber=True)


    pyfits.HDUList(stars).writeto("stars.fits", clobber=True)