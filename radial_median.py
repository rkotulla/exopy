#!/usr/bin/env python

import os
import sys
import pyfits
import numpy



if __name__ == "__main__":

    fn = sys.argv[1]

    cx = float(sys.argv[2])
    cy = float(sys.argv[3])

    print fn, cx, cy

    hdulist = pyfits.open(fn)
    data = hdulist[0].data

    iy,ix = numpy.indices(data.shape, dtype=numpy.float)
    print iy[:5,:5]
    print ix[:5,:5]

    # compute relative X and Y coordinates
    rel_y = iy - cy
    rel_x = ix - cx

    # based on that, compute polar coordinates r and phi
    print("computing polar coordinates")
    r = numpy.hypot(rel_x, rel_y)
    phi = numpy.degrees(numpy.arctan2(rel_x, rel_y))

    # Now just for practice, extract a radial ring
    max_r = 300
    min_r = 50
    in_ring = (r > min_r) & (r < max_r)
    data[~in_ring] = numpy.NaN #numpy.degrees(phi)[~in_ring]
    #pyfits.PrimaryHDU(data=data).writeto("ring.fits", clobber=True)

    #
    # apply cutoff at r=r_max
    #
    use_data = in_ring
    workdata = data[use_data]
    work_r = r[use_data]
    work_phi = phi[use_data]

    #print workdata.shape

    polar_median = numpy.zeros_like(workdata)
    polar_count = numpy.zeros_like(workdata)

    dr = 15
    dphi = 1
    dpix = 4
    for pixel in range(workdata.shape[0]):

        if ((pixel % 1000) == 0):
            print pixel,"of",workdata.shape[0]

        _r = work_r[pixel]
        _phi = work_phi[pixel]


        #
        # To speed up processing, select neighboring pixels in 2 steps,
        # first in phi and then in r
        #
        dphi2 = dphi #numpy.degrees(numpy.arcsin(dpix/(_r-dr)))
        #print dphi2
        right_angle = (work_phi > _phi-dphi2) & (work_phi < _phi+dphi2)

        this_r = work_r[right_angle]
        this_phi = work_phi[right_angle]
        this_data = workdata[right_angle]

        right_r = (this_r > (_r-dr)) & (this_r  < (_r+dr))
        this_r = this_r[right_r]
        this_phi = this_phi[right_r]
        this_data = this_data[right_r]

        # Now make sure we also fulfill the criteria to be part of the
        # selection "box"
        delta_phi = this_phi - _phi
        pixel_offset = numpy.abs(numpy.sin(numpy.radians(delta_phi))*this_r)
        final_select = (pixel_offset < dpix)

        polar_median[pixel] = numpy.nanmedian(this_data[final_select])
        polar_count[pixel] = numpy.sum(final_select)

        # neighbors = (work_r > _r-dr) & (work_r < _r+dr) & \
        #             (work_phi > _phi-dphi) & (work_phi < _phi+dphi)
        # polar_median[pixel] = numpy.nanmedian(workdata[neighbors])

        # numpy.savetxt("pixel_%d" % (pixel+1),
        #               numpy.array([this_r[final_select],
        #                            this_phi[final_select],
        #                            this_data[final_select]]).T)
        # if (pixel > 25):
        #     break

    data[use_data] = polar_median
    pyfits.PrimaryHDU(data=data).writeto("smoothring.fits", clobber=True)

    data[use_data] = polar_count
    pyfits.PrimaryHDU(data=data).writeto("smoothring_count.fits", clobber=True)
