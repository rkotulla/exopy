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

    sextractor_conf = sys.argv[3]
    sextractor_params = sys.argv[4]

    matching_radius_arcsec = 1.
    for frame_number in range(1000):
        filebase = "%s/mag%.1f__chunk%02d" % (
            out_dir, mag, frame_number+1
        )

        fits_fn = filebase+".fits"
        cat_fn = filebase+".cat"

        if (not os.path.isfile(fits_fn)):
            break

        #
        # Now run source extractor
        #
        sex_cat_fn = filebase+".sexsrc"
        sex_cmd = "sex -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s %s" % (
            sextractor_conf, sex_cat_fn, sextractor_params, fits_fn
        )
        print sex_cmd

        ret = subprocess.Popen(sex_cmd.split(),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        sextractor_pid = ret.pid
        (sex_stdout, sex_stderr) = ret.communicate()
        if (ret.returncode != 0):
            print("Sextractor might have a problem, check the log")
            print("Stdout=\n"+sex_stdout)
            print("Stderr=\n"+sex_stderr)
            continue


        #
        # Load the catalog of reference positions
        #
        inserted_positions = numpy.loadtxt(cat_fn)

        #
        # Also load the list of sources found by Sextractor
        #
        sex_cat = numpy.loadtxt(sex_cat_fn)

        #
        # Do some prep work before matching catalogs
        #
        hdu = pyfits.open(fits_fn)
        pixelscale = hdu[0].header['PIXLSCAL']
        magzero = hdu[0].header['MAGZERO']
        n_insert = hdu[0].header['_NINSERT']
        starmag = hdu[0].header['_STARMAG']
        matching_radius_pixels = matching_radius_arcsec / pixelscale
        cx = hdu[0].header['_CENTERX']
        cy = hdu[0].header['_CENTERY']

        #
        # Match the two catalogs using KD-Trees
        #
        ref_tree = scipy.spatial.cKDTree(sex_cat[:,0:2])
        d,i = ref_tree.query(
            x=inserted_positions[:, 0:2],
            k=1, # max 5 matches
            p=2,
            distance_upper_bound=matching_radius_pixels,
        )
        # print i
        print d.shape
        matched = numpy.isfinite(d)
        print numpy.sum(matched)

        matched_ref = inserted_positions[matched]
        matched_src = sex_cat[i[matched]]
        print matched_ref.shape, matched_src.shape

        matched_combined = numpy.append(matched_ref, matched_src, axis=1)
        print matched_combined.shape

        # apply photometric zeropoint
        matched_combined[:,5] += magzero

        #
        # compute additional columns that might be relevant
        #
        delta_mag = matched_combined[:,5]-starmag
        r_center = numpy.hypot(matched_combined[:,0]-cx,
                               matched_combined[:,1]-cy)
        extra_columns = numpy.array([delta_mag, r_center]).T
        matched_combined = numpy.append(
            matched_combined, extra_columns,
            axis=1
        )

        # true match: magnitude within 2*phot error
        # matched_combined = numpy.append(matched_combined,
        #                                 (matched_combined[:,
        #                                  4]-starmag).reshape((-1,1)), axis=1)
        combined_fn = filebase+".recover"
        numpy.savetxt(combined_fn, matched_combined)

        # true_match = numpy.fabs(matched_combined[:,4] - starmag) < \
        #              2*matched_combined[:,5]
        # valid_matched = matched_combined[true_match]
        # print valid_matched.shape
