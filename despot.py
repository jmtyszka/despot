#!/usr/bin/env python3
"""
Perform DESPOT modeling of combined bSSFP, SPGR and IR-SPGR images
- models T1, T2, rho, omega, B1 from multiflip imaging
- expects BIDS format data

Usage
----
despot.py -i <BIDS source directory>
despot.py -h

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2017-01-25 JMT From scratch

License
----
This file is part of atlaskit.

    atlaskit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atlaskit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with atlaskit.  If not, see <http://www.gnu.org/licenses/>.

Copyright
----
2017 California Institute of Technology.
"""

__version__ = '0.2.0'

import os
import sys
import argparse
import json
import nibabel as nib
import numpy as np
from scipy import optimize
from bids.grabbids import BIDSLayout


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='DESPOT reconstruction of bSSFP, SPGR, IR-SPGR dataset')
    parser.add_argument('-i','--bidsdir', required=True, help='BIDS format source data')

    # Parse command line arguments
    args = parser.parse_args()

    bidsdir = args.bidsdir

    layout = BIDSLayout(bidsdir)

    for sub in layout.get_subjects():

        for ses in layout.get_sessions():

            print('Subject %s : Session %s' % (sub, ses))

            spgr_tr, spgr_te, spgr_alpha, spgr_img = load_spgr(layout, sub, ses)

    # Clean exit
    sys.exit(0)


def load_spgr(layout, sub, ses):

    spgr_objs = layout.get(subject=sub, session=ses, type='^SPGR$', extensions='.nii.gz')
    n_spgr = len(spgr_objs)

    print('Found %d SPGR images' % n_spgr)

    # Allocate SPGR parameters
    spgr_tr = np.zeros(n_spgr)
    spgr_te = np.zeros(n_spgr)
    spgr_alpha = np.zeros(n_spgr)

    # Load SPGR data
    # Note use of anchors in regex to force exact matching ('^expr$')
    for i, spgr_obj in enumerate(spgr_objs):

        # Use in-house metadata function until pybids updated for session subdirs
        md_dict = get_metadata(spgr_obj)

        # Imaging parameters required for DESPOT model
        spgr_tr[i] = np.float(md_dict['RepetitionTime']) * 1000.0
        spgr_te[i] = np.float(md_dict['EchoTime']) * 1000.0
        spgr_alpha[i] = np.float(md_dict['FlipAngle'])

        print('(TR, TE, alpha) : (%0.2f, %0.2f, %0.1f)' % (spgr_tr[i], spgr_te[i], spgr_alpha[i]))

        # Load SPGR image data


def get_metadata(img_obj):

    json_fname = img_obj.filename.replace('.nii.gz','.json').replace('.nii','.json')
    if os.path.exists(json_fname):
        md_dict = json.load(open(json_fname, "r"))

    return md_dict


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
