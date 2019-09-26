#!/usr/bin/env python
"""
This file reproduces the same basic functionality of find_particles.py
in the yt.lagrangian_volume repository, but it allows one to track
multiple haloes simultaneously and spit the outputs to a file.

Usage
-----
python find_particle_list.py <pf_later> <pf_early> <input_text> <output_text>

where <input_text> is a file where each line is a halo id from <pf_later> to
be tracked and given a bounding box output to <output_text>.
"""

from yt.mods import *
from yt.analysis_modules.halo_finding.api import *
from os import environ
environ['CFLAGS'] = "-I"+np.get_include()

import pyximport; pyximport.install()
import particle_ops
import argparse
import find_particles as fp

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("later", type=str)
    parser.add_argument("earlier", type=str)
    parser.add_argument("input_file", type=str)
    parser.add_argument("output_file", type=str)
    args = parser.parse_args()
    pf1 = load(args.later)
    pf2 = load(args.earlier)
    f_in = open(args.input_file, 'r') # each line should be a halo num
    f_out = open(args.output_file, 'w')
    if pf1 is None or pf2 is None:
        print "Sorry, couldn't load your parameter files!"
        sys.exit(1)
    halos = FOFHaloFinder(pf1)
    for halo in f_in:
        halo_id = int(halo)
        ind = halos[halo_id]["particle_index"]
        min_pos, max_pos = fp.correlate_particles(ind, pf2)

        # Now calculate some info for output
        center_later = halos[halo_id].center_of_mass()
        box_size = pf2.domain_right_edge - pf2.domain_left_edge
        box_center = (pf2.domain_right_edge + pf2.domain_left_edge) / 2.0
        new_min_pos = (min_pos - center_later) * box_size + box_center
        new_max_pos = (max_pos - center_later) * box_size + box_center
        bound_volume = 1.
        for i in range(3):
            bound_volume *= (max_pos[i] - min_pos[i])

        f_out.write("Halo id: %i\n" % halo_id)
        f_out.write("The bounding box:\n")
        f_out.write("    %0.8e %0.8e %0.8e\n" % tuple(min_pos))
        f_out.write("    %0.8e %0.8e %0.8e\n" % tuple(max_pos))
        f_out.write("\n")
        f_out.write("centered on COM location of halo in later pf:\n")
        f_out.write("    %0.8e %0.8e %0.8e\n" % tuple(center_later))
        f_out.write("\n")
        f_out.write("which yields a boundary box in the new frame:\n")
        f_out.write("    %0.8e %0.8e %0.8e\n" % tuple(new_min_pos))
        f_out.write("    %0.8e %0.8e %0.8e\n" % tuple(new_max_pos))
        f_out.write("\n")
        f_out.write("with a total bound volume of %0.8e\n" % bound_volume)
        f_out.write("---------------------------------------------\n")
        f_out.write("\n")
        f_out.write("\n")
