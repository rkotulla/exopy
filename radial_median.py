#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
import multiprocessing

def parallel_worker(jobqueue,
                    data, r, phi,
                    dr, dphi, dpix,
                    return_queue,):

    while (True):

        pixel = jobqueue.get()
        if (pixel is None):
            jobqueue.task_done()
            return

        if ((pixel % 1000) == 0):
            print pixel,"of",data.shape[0]

        _r = r[pixel]
        _phi = phi[pixel]

        median = numpy.NaN
        count = 0

        try:
            #
            # To speed up processing, select neighboring pixels in 2 steps,
            # first in phi and then in r
            #
            dphi2 = dphi #numpy.degrees(numpy.arcsin(dpix/(_r-dr)))
            #print dphi2
            right_angle = (phi > _phi - dphi2) & (phi < _phi + dphi2) & numpy.isfinite(data)

            this_r = r[right_angle]
            this_phi = phi[right_angle]
            this_data = data[right_angle]

            right_r = (this_r > (_r-dr)) & (this_r  < (_r+dr))
            this_r = this_r[right_r]
            this_phi = this_phi[right_r]
            this_data = this_data[right_r]

            # Now make sure we also fulfill the criteria to be part of the
            # selection "box"
            delta_phi = this_phi - _phi
            pixel_offset = numpy.abs(numpy.sin(numpy.radians(delta_phi))*this_r)
            final_select = (pixel_offset < dpix)

            # polar_median[pixel] = numpy.nanmedian(this_data[final_select])
            # polar_count[pixel] = numpy.sum(final_select)
            median = numpy.nanmedian(this_data[final_select])
            count = numpy.sum(final_select)
        except:
            pass

        return_queue.put([pixel, median, count])
        jobqueue.task_done()

        # neighbors = (work_r > _r-dr) & (work_r < _r+dr) & \
        #             (work_phi > _phi-dphi) & (work_phi < _phi+dphi)
        # polar_median[pixel] = numpy.nanmedian(workdata[neighbors])

        # numpy.savetxt("pixel_%d" % (pixel+1),
        #               numpy.array([this_r[final_select],
        #                            this_phi[final_select],
        #                            this_data[final_select]]).T)
        # if (pixel > 25):
        #     break



if __name__ == "__main__":

    fn = sys.argv[1]

    cx = float(sys.argv[2])
    cy = float(sys.argv[3])

    print fn, cx, cy

    hdulist = pyfits.open(fn)
    # check what extension to use
    try:
        extname = sys.argv[4]
    except:
        extname = "PRIMARY"
    data = hdulist[extname].data

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
    max_r = 700
    min_r = 50
    in_ring = (r > min_r) & (r < max_r)
    data[~in_ring] = numpy.NaN #numpy.degrees(phi)[~in_ring]
    #pyfits.PrimaryHDU(data=data).writeto("ring.fits", clobber=True)

    #
    # apply cutoff at r=r_max
    #
    use_data = in_ring
    work_data = data[use_data]
    work_r = r[use_data]
    work_phi = phi[use_data]

    #print workdata.shape

    polar_median = numpy.zeros_like(work_data)
    polar_count = numpy.zeros_like(work_data)

    dr = 15
    dphi = 1
    dpix = 4

    jobqueue = multiprocessing.JoinableQueue()
    resultqueue = multiprocessing.Queue()

    for pixel in range(work_data.shape[0]):
        jobqueue.put(pixel)

    n_workers = multiprocessing.cpu_count()
    processes = []
    for i in range(n_workers):
        jobqueue.put(None)
    for i in range(n_workers):
        p = multiprocessing.Process(
            target=parallel_worker,
            kwargs=dict(
                jobqueue=jobqueue,
                data=work_data, r=work_r, phi=work_phi,
                dr=dr, dphi=dphi, dpix=dpix,
                return_queue=resultqueue,
            )
        )
        p.start()
        processes.append(p)

    # read all results
    for i in range(work_data.shape[0]):

        (pixel,median,count) = resultqueue.get()
        polar_median[pixel] = median
        polar_count[pixel] = count


        # if ((pixel % 1000) == 0):
        #     print pixel,"of",data.shape[0]
        #
        # _r = r[pixel]
        # _phi = phi[pixel]
        #
        #
        # #
        # # To speed up processing, select neighboring pixels in 2 steps,
        # # first in phi and then in r
        # #
        # dphi2 = dphi #numpy.degrees(numpy.arcsin(dpix/(_r-dr)))
        # #print dphi2
        # right_angle = (phi > _phi - dphi2) & (phi < _phi + dphi2)
        #
        # this_r = r[right_angle]
        # this_phi = phi[right_angle]
        # this_data = data[right_angle]
        #
        # right_r = (this_r > (_r-dr)) & (this_r  < (_r+dr))
        # this_r = this_r[right_r]
        # this_phi = this_phi[right_r]
        # this_data = this_data[right_r]
        #
        # # Now make sure we also fulfill the criteria to be part of the
        # # selection "box"
        # delta_phi = this_phi - _phi
        # pixel_offset = numpy.abs(numpy.sin(numpy.radians(delta_phi))*this_r)
        # final_select = (pixel_offset < dpix)
        #
        # polar_median[pixel] = numpy.nanmedian(this_data[final_select])
        # polar_count[pixel] = numpy.sum(final_select)
        #
        # # neighbors = (work_r > _r-dr) & (work_r < _r+dr) & \
        # #             (work_phi > _phi-dphi) & (work_phi < _phi+dphi)
        # # polar_median[pixel] = numpy.nanmedian(workdata[neighbors])
        #
        # # numpy.savetxt("pixel_%d" % (pixel+1),
        # #               numpy.array([this_r[final_select],
        # #                            this_phi[final_select],
        # #                            this_data[final_select]]).T)
        # # if (pixel > 25):
        # #     break

    _data = data.copy()
    _data[use_data] = polar_median
    pyfits.PrimaryHDU(data=_data).writeto("smoothring.fits", clobber=True)

    _data = data.copy()
    _data[use_data] = polar_count
    pyfits.PrimaryHDU(data=_data).writeto("smoothring_count.fits", clobber=True)

    #
    # Finally, subtract the smoothed model from the input data
    #
    data[use_data] -= polar_median
    out_fn = fn[:-5]+".radialbgsub.fits"
    hdulist[extname].data = data
    # pyfits.PrimaryHDU(data=data).writeto(out_fn, clobber=True)
    hdulist.writeto(out_fn, clobber=True)
