# -*- coding: utf-8 -*-
##    Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
##
##    This file is part of Nascent.
##
##    Nascent is free software: you can redistribute it and/or modify
##    it under the terms of the GNU Affero General Public License as
##    published by the Free Software Foundation, either version 3 of the
##    License, or (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU Affero General Public License for more details.
##
##    You should have received a copy of the GNU Affero General Public License
##    along with this program. If not, see <http://www.gnu.org/licenses/>.


"""
Refs:
* http://colorbrewer2.com/
* https://vis4.net/blog/posts/avoid-equidistant-hsv-colors/
* http://zzolo.org/colortools/
* http://scholarship.claremont.edu/cgi/viewcontent.cgi?article=1881&context=cmc_theses
* http://alumni.media.mit.edu/~wad/color/palette.html
* http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors


"""

from random import random
import colorsys

from .debug import printd, pprintd


class StrandColor():
    def __init__(self, randomize=False):
        self.st_spec_colors = {}
        self.last_hue = 0
        self.hue_delta = 1/6  # rygcbm
        self.used = []
        self.randomize = randomize
        self.rounds = 0
        self.mindist = 0.1

    def strand_color(self, strand):
        """ Return a specie-specific strand color. """
        N_MAX_ABORT = 200
        n_tries = 0
        if strand.name not in self.st_spec_colors:
            printd("last_hue: %s, self.mindist=%s" % (self.last_hue, self.mindist))
            hue = random() if self.randomize else self.last_hue
            while any(abs(hue-h2) < self.mindist for h2 in self.used):
                n_tries += 1
                hue += self.hue_delta
                if hue >= 1:
                    self.hue_delta /= 2  #1/(1/self.hue_delta + 1)
                    #self.hue_delta /= (1/self.hue_delta + 1)
                    self.rounds += 1
                    self.mindist /= 2
                    #self.mindist = 1/self.rounds
                    hue %= 1
                if n_tries > N_MAX_ABORT:
                    hue = random()
                    print("UNABLE TO FIND NEW HUE (self.rounds = %s, self.mindist=%s), RETURNING RANDOM HUE: %s"
                          % (self.rounds, self.mindist, hue))
                    break
            self.last_hue = hue
            self.st_spec_colors[strand.name] = hue
            self.used.append(hue)
            printd("Hue for strand %s: %s" % (strand, hue))
            # For 0 < frac < 0.5, use frac*2 as value, max out at 1.
            # e.g. frac=1/4 (3 domains), val=0.5; frac=1/8 (7 domains), val=1/4, etc.
        else:
            hue = self.st_spec_colors[strand.name]
        return hue

    def domain_hsv(self, strand, frac):
        """
        Frac is (domain_idx+1)/(len(strand.domains)+1)
        So that for strand with 1 domain, frac=0.5,
        strand with 2 domains, frac=1/3 or 2/3.
        3 domains: frac = 1/4, 2/4, 3/4,
        etc.
        """
        assert 0 < frac < 1
        hue = self.strand_color(strand)
        val = min((frac*2, 1.0))
        sat = min((2-frac*2, 1.0))
        printd("Strand %s with frac %s yields hsv = %s" % (strand, frac, (hue, sat, val)))
        return hue, sat, val

    def domain_rgb(self, strand, frac):
        """
        """
        hsv = self.domain_hsv(strand, frac)
        rgb = colorsys.hsv_to_rgb(*hsv)
        printd("Strand %s with frac %s yields rgb = %s" % (strand, frac, rgb))
        return rgb


class StrandColorRandom(StrandColor):
    def strand_color(self, strand):
        """ Return a specie-specific strand color. """
        # Generate a large random sample, pick the number that is furthest away from existing values:
        if strand.name in self.st_spec_colors:
            hue = self.st_spec_colors[strand.name]
        else:
            if not self.used:
                hue = 0
            else:
                candidates = [(i, random()) for i in range(1000)]
                bestidx, hue = max(candidates, key=lambda tup: min(abs(tup[1]-h2) for h2 in self.used))
            self.st_spec_colors[strand.name] = hue
            self.used.append(hue)
        printd("Hue for strand %s: %s" % (strand, hue))
        return hue
