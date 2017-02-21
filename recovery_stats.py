#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
import math
import subprocess
import scipy.spatial

if __name__ == "__main__":

    out_dir = sys.argv[1]
    mag = float(sys.argv[2])

    inserted_all = None
    recovered_all = None

    for frame_number in range(1000):
        filebase = "%s/mag%.1f__chunk%02d" % (
            out_dir, mag, frame_number+1
        )

        inserted_fn = filebase+".cat"
        recovered_fn = filebase+".recover"

        if (not os.path.isfile(inserted_fn) or
            not os.path.isfile(recovered_fn)):
            break

        inserted = numpy.loadtxt(inserted_fn)
        inserted_all = inserted if inserted_all is None else \
            numpy.append(inserted_all, inserted, axis=0)

        recovered = numpy.loadtxt(recovered_fn)
        recovered_all = recovered if recovered_all is None else \
            numpy.append(recovered_all, recovered, axis=0)

    numpy.savetxt("inserted.all", inserted_all)
    numpy.savetxt("recovered.all", recovered_all)

    bins = numpy.linspace(0, 400, 1)
    print bins

    count_inserted,_ = numpy.histogram(
        inserted_all[:,2], bins=bins,
    )
    count_recovered,_ = numpy.histogram(
        recovered_all[:,2], bins=bins,
    )
    bin_center = 0.5*(bins[0:-1] + bins[1:])
    bin_width = numpy.diff(bins)

    numpy.savetxt("analysis_mag%.1f" % (mag),
        numpy.array([bin_center, bin_width, count_inserted,
                     count_recovered]).T)