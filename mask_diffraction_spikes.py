#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
import multiprocessing
import time
import math
from optparse import OptionParser

from astLib import astWCS
import ephem

def mask_diffraction_spikes(data, hdr, x, y, width, length, write_region=None):

    # rotator angle
    rotator_angle = hdr['ROTSTART']
    if (rotator_angle > 180.): rotator_angle -= 360.
    if (rotator_angle < -180.): rotator_angle += 360.
    elevation = numpy.degrees(float(hdr['ELMAP']))
    # print rotator_angle, -1*rotator_angle, elevation

    rotate_angle = (rotator_angle - elevation)

    #
    # Select a region to speed up things by avoiding having to compare too
    # many pixels
    #
    ix1 = int(math.floor(x-length))
    ix2 = int(math.ceil(x+length))
    iy1 = int(math.floor(y-length))
    iy2 = int(math.ceil(y+length))
    iy,ix = numpy.indices(data.shape, dtype=numpy.float)

    _ix = ix[iy1:iy2, ix1:ix2] - x
    _iy = iy[iy1:iy2, ix1:ix2] - y
    _phi = numpy.arctan2(_ix, _iy)
    _phi[_phi<0] += 2*math.pi
    _r = numpy.hypot(_ix, _iy)

    mask = numpy.zeros(data.shape, dtype=numpy.bool)

    #
    # Create a ds9 region file drawing lines where the diffraction spikes are
    #
    if (write_region is not None):
        ds9 = open(write_region, "w")
        print >>ds9, """\
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image"""

    angles = numpy.arange(4) * 90. + 45. - rotate_angle
    print "Diffraction spikes at %s" % (",".join(["%.1f" % (a) for a in angles]))

    for i, angle in enumerate(numpy.fmod(numpy.radians(angles), 2*math.pi)):

        # print angle, numpy.degrees(angle)

        dx = length * numpy.sin(angle)
        dy = length * numpy.cos(angle)
        if (write_region is not None):
            print >>ds9, "line(%.2f,%.2f,%.2f,%.2f) # line=0 0 text={%.2f}" % (
                x,y,x+dx,y+dy,numpy.degrees(angle),
            )

        angle_rad = numpy.radians(angle)
        # dist_x = _r * numpy.cos(math.pi/2-(angle-_phi))
        dist_x = _r * numpy.sin(angle-_phi)

        _data = data.copy()
        #_data[iy1:iy2, ix1:ix2][part_of_spike] = numpy.NaN
        #d_angle = numpy.fmod((_phi-angle+2*math.pi), 2*math.pi)

        d_angle = (_phi-angle)
        d_angle[d_angle > math.pi] -= 2*math.pi
        d_angle[d_angle < -1*math.pi] += 2*math.pi
        d_angle = numpy.fabs(d_angle)

        part_of_spike = (numpy.fabs(dist_x) < width) & (d_angle < math.pi/2) & (_r < length)
            #(numpy.fmod(numpy.fabs(_phi-angle),2*math.pi) < math.pi/2)

        # data[iy1:iy2, ix1:ix2][part_of_spike] = numpy.NaN
        mask[iy1:iy2, ix1:ix2][part_of_spike] = True

        #_data[iy1:iy2, ix1:ix2] = d_angle #(_phi-angle) #numpy.fmod(numpy.fabs(_phi-angle), 2*math.pi)

        # data[iy1:iy2, ix1:ix2] = numpy.fabs(dist_x)
        #line(5:19:12.062,+40:06:35.516,5:19:08.873,+40:05:48.712) # line=0 0

        # hdulist[extname].data = _data
        # hdulist.writeto(fn[:-5]+".spikemask%d.fits" % (i+1), clobber=True)

    return mask
    #

if __name__ == "__main__":


    parser = OptionParser()
    parser.add_option("", "--maxr", dest="maxr",
                      default="700", type=int)
    parser.add_option("", "--dr", dest="dr",
                      default=15, type=float)
    parser.add_option("", "--dpix", dest="dpix",
                       default=4., type=float)
    parser.add_option("", "--dmask", dest="dmask",
                       default=None, type=float)
    parser.add_option("", "--pixel", dest="pixelcoords",
                      action="store_true", default=False)
    parser.add_option("", "--extname", dest="extname",
                      default='PRIMARY', type=str)

    (options, cmdline_args) = parser.parse_args()


    x_ra = cmdline_args[0]
    y_dec = cmdline_args[1]

    print options.extname

    for fn in cmdline_args[2:]:
        print "\n\nWorking on %s\n\n" % (fn)

        hdulist = pyfits.open(fn)
        ext = hdulist[options.extname]
        data = ext.data
        hdr = hdulist[0].header

        write_region = fn[:-5]+".diffspikes.reg"

        mask = mask_diffraction_spikes(
            data=data, hdr=hdr,
            # fn=fn, extname=options.extname,
            x=float(x_ra), y=float(y_dec),
            write_region=write_region,
            width=10, length=750,
        )

        data[mask] = numpy.NaN
        hdulist[options.extname].data = data
        hdulist.writeto(fn[:-5]+".spikemask.fits", clobber=True)