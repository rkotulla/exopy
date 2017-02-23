#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
import multiprocessing
import time
import math
from optparse import OptionParser



def parallel_worker(jobqueue,
                    data, r, phi,
                    dr, dphi, dpix,
                    xy_list,
                    return_queue,
                    chunksize=1000,):

    pixelnumber = 0
    max_box_size = int(math.ceil(numpy.hypot(dr, dpix)))

    while (True):

        chunk = jobqueue.get()
        if (chunk is None):
            jobqueue.task_done()
            print "Shutting down worker"
            break

        c1 = chunk*chunksize
        c2 = c1 + chunksize
        if (c2 > xy_list.shape[0]): c2=xy_list.shape[0]
        pixel_list = xy_list[c1:c2,:]

        median = numpy.zeros((c2-c1))
        median[:] = numpy.NaN

        count = numpy.zeros((c2-c1))

        print chunk, median.shape, count.shape,

        for i_pixel, pixel in enumerate(pixel_list):

            #print pixel
            x,y = pixel[0], pixel[1]

            _r = r[y,x]
            _phi = phi[y,x]
            _data = data[y,x]
            #print pixel, x, y, _r, _phi, _data

            #
            # Work out the maximum dimension of a box around the current position
            # that contains all pixels we need to work on
            #
            box_size_y = int(math.ceil(max_box_size * numpy.fabs(numpy.cos(numpy.radians(_phi))))) + dpix
            box_size_x = int(math.ceil(max_box_size * numpy.fabs(numpy.sin(numpy.radians(_phi))))) + dpix

            #median[i_pixel] = box_size_y #numpy.nanmedian(this_data[final_select])
            #count[i_pixel] = box_size_x #numpy.sum(final_select)
            #continue

            y1 = y - box_size_y
            y2 = y + box_size_y
            x1 = x - box_size_x
            x2 = x + box_size_x

            block_data = data[y1:y2, x1:x2]
            block_phi = phi[y1:y2, x1:x2]
            block_r = r[y1:y2, x1:x2]


            pixelnumber += 1
            if ((pixelnumber % 100) == 0):
                sys.stdout.write(".")
                # sys.stdout.write("\r%6d %4d %4d" % (pixelnumber, x,y))
                sys.stdout.flush()

            #_r = r[pixel]
            #_phi = phi[pixel]

            if (True): #try:
                #
                # To speed up processing, select neighboring pixels in 2 steps,
                # first in phi and then in r
                #
                #dphi2 = dphi #numpy.degrees(numpy.arcsin(dpix/(_r-dr)))
                #print dphi2

                right_angle = (block_phi > _phi - dphi) & \
                              (block_phi < _phi + dphi) & \
                              (block_r > (_r-dr)) & \
                              (block_r < (_r+dr)) & \
                              numpy.isfinite(block_data)

                this_r = block_r[right_angle]
                this_phi = block_phi[right_angle]
                this_data = block_data[right_angle]

                # right_r = (this_r > (_r-dr)) & (this_r  < (_r+dr))
                # this_r = this_r[right_r]
                # this_phi = this_phi[right_r]
                # this_data = this_data[right_r]

                # Now make sure we also fulfill the criteria to be part of the
                # selection "box"
                delta_phi = this_phi - _phi
                pixel_offset = numpy.abs(numpy.sin(numpy.radians(delta_phi))*this_r)
                final_select = (pixel_offset < dpix)

                # polar_median[pixel] = numpy.nanmedian(this_data[final_select])
                # polar_count[pixel] = numpy.sum(final_select)
                median[i_pixel] = numpy.nanmedian(this_data[final_select])
                count[i_pixel] = numpy.sum(final_select)
            # except:
            #     pass

        return_queue.put((chunk, c1, c2, median, count))
        jobqueue.task_done()
        print
            #break


            # neighbors = (work_r > _r-dr) & (work_r < _r+dr) & \
            #             (work_phi > _phi-dphi) & (work_phi < _phi+dphi)
            # polar_median[pixel] = numpy.nanmedian(workdata[neighbors])

            # numpy.savetxt("pixel_%d" % (pixel+1),
            #               numpy.array([this_r[final_select],
            #                            this_phi[final_select],
            #                            this_data[final_select]]).T)
            # if (pixel > 25):
            #     break

    print "terminating worker"


if __name__ == "__main__":


    parser = OptionParser()
    parser.add_option("", "--maxr", dest="maxr",
                      default="700", type=int)
    parser.add_option("", "--dr", dest="dr",
                      default=15, type=float)
    parser.add_option("", "--dpix", dest="dpix",
                       default=4., type=float)
    (options, cmdline_args) = parser.parse_args()

    fn = cmdline_args[0]
    cx = float(cmdline_args[1])
    cy = float(cmdline_args[2])

    print fn, cx, cy

    hdulist = pyfits.open(fn)
    # check what extension to use
    try:
        extname = sys.argv[4]
    except:
        extname = "PRIMARY"
    data = hdulist[extname].data

    iy,ix = numpy.indices(data.shape)
    #print iy[:5,:5]
    #print ix[:5,:5]

    # compute relative X and Y coordinates
    rel_y = iy.astype(numpy.float) - cy
    rel_x = ix.astype(numpy.float) - cx

    # based on that, compute polar coordinates r and phi
    print("computing polar coordinates")
    r = numpy.hypot(rel_x, rel_y)
    phi = numpy.degrees(numpy.arctan2(rel_x, rel_y))

    # Now just for practice, extract a radial ring
    max_r = 1000
    min_r = 50
    in_ring = (r > min_r) & (r < max_r)
    # data[~in_ring] = numpy.NaN #numpy.degrees(phi)[~in_ring]
    #pyfits.PrimaryHDU(data=data).writeto("ring.fits", clobber=True)

    x_work = ix[in_ring]
    y_work = iy[in_ring]

    print x_work.shape, y_work.shape
    xy_work = numpy.array([x_work, y_work]).T
    print xy_work.shape


    #
    #
    # os._exit(0)
    #
    #
    # #
    # # apply cutoff at r=r_max
    # #
    use_data = in_ring
    # work_data = data[use_data]
    # work_r = r[use_data]
    # work_phi = phi[use_data]
    #
    # #print workdata.shape

    polar_median = numpy.zeros((xy_work.shape[0]))
    polar_count = numpy.zeros((xy_work.shape[0]))

    dr = 15
    dphi = 1
    dpix = 4

    jobqueue = multiprocessing.JoinableQueue()
    resultqueue = multiprocessing.Queue()

    chunksize = 5000
    chunk_start = 0
    chunk_end = 0
    n_chunks = int(math.ceil(float(xy_work.shape[0])/float(chunksize)))
    print xy_work.shape[0], n_chunks
    print "Feeding worker queue",n_chunks
    for n in range(n_chunks):
        jobqueue.put(n)
        # while (chunk_end < xy_work.shape[0]):
        #     chunk_end = chunk_start + chunksize
        #     if (chunk_end >= xy_work.shape[0]):
        #         chunk_end = xy_work.shape[0]
        #     jobqueue.put(xy_work[chunk_start:chunk_end, :])
        #     chunk_start = chunk_end
        sys.stdout.write(".")
        sys.stdout.flush()

    print "starting worker"
    n_workers = multiprocessing.cpu_count()
    for i in range(n_workers):
        jobqueue.put(None)

    processes = []
    for i in range(n_workers):
        p = multiprocessing.Process(
            target=parallel_worker,
            kwargs=dict(
                jobqueue=jobqueue,
                data=data, r=r, phi=phi,
                dr=dr, dphi=dphi, dpix=dpix,
                xy_list=xy_work, chunksize=chunksize,
                return_queue=resultqueue,
            )
        )
        p.daemon = True
        p.start()
        processes.append(p)

    # for i,pixel in enumerate(xy_work): #range(xy_work.shape[0]):
    #     jobqueue.put(pixel)
    #
    # print "starting to wait beforre quitting"
    # time.sleep(3)
    # os._exit(0)

    # read all results
    #radial_median = numpy.zeros((xy_work.shape[0]))
    #radial_count = numpy.zeros((xy_work.shape[0]))

    # time.sleep(5)
    #
    # print "joining processes"
    # print processes
    # for p in processes:
    #     p.join()
    #     # if (p.is_alive()):
    #     #     print "terminating"
    #     #     p.terminate()
    #     #     p.join(timeout=0.01)
    #
    # print processes


    for i in range(n_chunks):

        print "reading results"
        (chunk, c1, c2, median, count) = resultqueue.get()
        #print "Received chunk %d" % (chunk), median.shape, count.shape

        polar_median[c1:c2] = median[:]
        polar_count[c1:c2] = count[:]



        # (pixel,median,count) = resultqueue.get()
        # polar_median[pixel] = median
        # polar_count[pixel] = count


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

    print use_data.shape
    print (data[use_data]).shape, polar_median.shape

    print "writing smoothring"
    _data = data.copy()
    _data[use_data] = polar_median
    pyfits.PrimaryHDU(data=_data).writeto("smoothring.fits", clobber=True)

    print "writing smoothring_count"
    _data = data.copy()
    _data[use_data] = polar_count
    pyfits.PrimaryHDU(data=_data).writeto("smoothring_count.fits", clobber=True)

    #
    # Finally, subtract the smoothed model from the input data
    #
    print "writing radialbgsub"
    data[use_data] -= polar_median
    out_fn = fn[:-5]+".radialbgsub.fits"
    hdulist[extname].data = data
    # pyfits.PrimaryHDU(data=data).writeto(out_fn, clobber=True)
    hdulist.writeto(out_fn, clobber=True)
