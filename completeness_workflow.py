#!/usr/bin/env python

import os
import sys
import pyfits
import numpy
import math


import completeness
import recover_stars
import recovery_stats

if __name__ == "__main__":

    fn = sys.argv[1]

    cx = float(sys.argv[2])
    cy = float(sys.argv[3])

    seeing_fwhm = float(sys.argv[4])
    mag = [float(m) for m in sys.argv[5].split(",")]
    out_dir = sys.argv[6]

    sextractor_conf = sys.argv[7]
    sextractor_params = sys.argv[8]
    matching_radius_arcsec = 1.

    analysis_basename = sys.argv[9]

    try:
        n_total = int(sys.argv[10])
    except:
        n_total = 1000.

    n_per_image = 100.
    max_position_offset = 300

    for cur_mag in mag:

        # insert sources
        completeness.insert_artificial_sources(
            fn=fn,
            cx=cx, cy=cy,
            seeing_fwhm=seeing_fwhm,
            mag=cur_mag,
            out_dir=out_dir,
            n_total=n_total,
            n_per_image=n_per_image,
            max_position_offset=max_position_offset,
        )

        # do recovery
        recover_stars.recover_stars(
            out_dir=out_dir,
            mag=cur_mag,
            sextractor_conf=sextractor_conf,
            sextractor_params=sextractor_params,
            matching_radius_arcsec=matching_radius_arcsec,
        )

        # run analysis
        results = recovery_stats.recovery_stats(out_dir, cur_mag)
        numpy.savetxt("%s_mag%.1f" % (analysis_basename, cur_mag), results)

