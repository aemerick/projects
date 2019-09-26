"""
Particle operations for Lagrangian Volume

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Author: Cameron Hummels <chummels@gmail.com>
Affiliation: University of Arizona
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.mods import *
from yt.analysis_modules.halo_finding.api import *
from os import environ
environ['CFLAGS'] = "-I"+na.get_include()

import pyximport; pyximport.install()
import particle_ops
import argparse

def correlate_particles(particle_ids, pf):
    # First we identify all the particles and find their maximum extent
    min_pos = pf.domain_right_edge.copy().astype(float)
    max_pos = pf.domain_left_edge.copy().astype(float)

    # We simultaneously compute the maximum extent of particles that are
    # shifted over by 1/2 a box length to take into account periodic 
    # boundary conditions
    min_periodic_pos = pf.domain_right_edge.copy().astype(float)
    max_periodic_pos = pf.domain_left_edge.copy().astype(float)
    box_size = pf.domain_right_edge - pf.domain_left_edge

    mask_to_get = na.zeros(particle_ids.size, dtype='int32')
    for g in pf.h.grids:
        if g.NumberOfParticles == 0: continue
        found_any, mask = particle_ops.mask_particles(
            particle_ids, g["particle_index"], mask_to_get)
        if found_any == 0: continue
        mask = mask.astype("bool")
        for i, ax in enumerate('xyz'):
            min_pos[i] = min(g["particle_position_%s" % ax][mask].min(),
                             min_pos[i])
            max_pos[i] = max(g["particle_position_%s" % ax][mask].max(),
                             max_pos[i])

            # now do the same calculation with the position shifted over 0.5
            # box length to the right
            g["particle_position_%s" % ax][mask] += (box_size[i] / 2.0) 
            overshoot = g["particle_position_%s" % ax] >= pf.domain_right_edge[i]
            g["particle_position_%s" % ax][overshoot] -= box_size[i]
            min_periodic_pos[i] = min(g["particle_position_%s" % ax][mask].min(),
                             min_periodic_pos[i])
            max_periodic_pos[i] = max(g["particle_position_%s" % ax][mask].max(),
                             max_periodic_pos[i])
    
    # calculate the extent of each bounding box for non-periodic and
    # periodic frames before transforming periodic back to original frame
    extent = max_pos - min_pos
    periodic_extent = max_periodic_pos - min_periodic_pos

    # shift periodic extrema back over half a box length to the left
    min_periodic_pos -= (box_size / 2.0)
    max_periodic_pos -= (box_size / 2.0)
    for i in range(3):
        if min_periodic_pos[i] < pf.domain_left_edge[i]:
            min_periodic_pos[i] += box_size[i]
#        if max_periodic_pos[i] < pf.domain_left_edge[i]:
        max_periodic_pos[i] += box_size[i]

#    # Print out what periodic and non-periodic bounds are
#    print "Min: %s" % min_pos
#    print "Max: %s" % max_pos
#    print "Min Periodic: %s" % min_periodic_pos
#    print "Max Periodic: %s" % max_periodic_pos

    # So which is the correct set of bounds?  The non-shifted or the shifted?
    # It depends on which one gives a smaller extent in their frame.
    for i in range(3):
        if periodic_extent[i] < extent[i]:
            min_pos[i] = min_periodic_pos[i]
            max_pos[i] = max_periodic_pos[i]
    
    return min_pos, max_pos

if __name__ == "__main__":
    """
    find_particles accepts two parameter files (early and late) and a halo 
    ID number from the later pf.  It runs the halo finder on the later pf and 
    then traces back all of the particles from the desired halo ID
    back to their positions in the early pf.  It determines the bounding box
    which encompasses all of these particles in the early pf and the center
    of this range.

    Usage
    -----
    $ find_particles <pf_later> <pf_earlier> <halo_id_later>

    Example
    -------
    To trace the most massive halo from a late redshift run to its
    location at high redshift:

    $ find_particles DD0063/DD0063 DD0000/DD0000 0
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("later", type=str)
    parser.add_argument("earlier", type=str)
    parser.add_argument("halo_id", type=int)
    args = parser.parse_args()
    pf1 = load(args.later)
    pf2 = load(args.earlier)
    if pf1 is None or pf2 is None:
        print "Sorry, couldn't load your parameter files!"
        sys.exit(1)
    halos = HaloFinder(pf1)
    #halos = FOFHaloFinder(pf1)
    ind = halos[args.halo_id]["particle_index"]
    min_pos, max_pos = correlate_particles(ind, pf2)

    # Now calculate some info for output
    center_later = halos[args.halo_id].center_of_mass()
    box_size = pf2.domain_right_edge - pf2.domain_left_edge
    box_center = (pf2.domain_right_edge + pf2.domain_left_edge) / 2.0
    new_min_pos = (min_pos - center_later) * box_size + box_center
    new_max_pos = (max_pos - center_later) * box_size + box_center

    print
    print "The bounding box:"
    print "    %0.8e %0.8e %0.8e" % tuple(min_pos)
    print "    %0.8e %0.8e %0.8e" % tuple(max_pos)
    print
    print "centered on COM location of halo in later pf:"
    print "    %0.8e %0.8e %0.8e" % tuple(center_later)
    print
    print "which yields a boundary box in the new frame:"
    print "    %0.8e %0.8e %0.8e" % tuple(new_min_pos)
    print "    %0.8e %0.8e %0.8e" % tuple(new_max_pos)
