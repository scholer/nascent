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

# pylint: disable=C0103,W0142

"""

Module for analysing sequences.

"""

import os
import math
import string

strand1 = "GGCACGACAGGTTTCC"
strand2 = "GGAAACCTGTCGTGCC"
min_domain_length = 1
slen = len(strand1)

assert len(strand1) == len(strand2)
assert len(strand1) in [2**i for i in range(10)]

WC = dict(zip("ATGC", "TACG"))
def rcompl(seq):
    return "".join(WC[n] for n in seq[::-1])

fnfmt = "duplex_{slen}bp-d{ndoms}.txt"
headers = "Strand	Domains	Sequence".split("\t")
fieldsep = "\t"
add_ndoms_postfix = True

print("pwd '.':", os.path.abspath('.'))

for n in range(int(math.log(len(strand1), 2)) - min_domain_length + 2):
    # n = number of times we break in half
    seq = strand1
    ndoms = 2**n
    postfix = "/%s" % ndoms if add_ndoms_postfix else ""
    strands = [{'Strand': 's1'+postfix}, {'Strand': 's2'+postfix}]
    dlen = int(slen / (2**n))
    dseqs1 = [seq[start:start+dlen] for start in range(0, slen, dlen)]
    assert len(dseqs1) == ndoms
    dseqs2 = [rcompl(s) for s in reversed(dseqs1)]
    dnames1 = [c+postfix for c in string.ascii_lowercase[:ndoms]]
    dnames2 = [c.upper() for c in reversed(dnames1)]
    strands[0]['Domains'] = ",".join(dnames1)
    strands[1]['Domains'] = ",".join(dnames2)
    strands[0]['Sequence'] = ",".join(dseqs1)
    strands[1]['Sequence'] = ",".join(dseqs2)
    with open(fnfmt.format(n=n, slen=slen, ndoms=ndoms), 'w') as fp:
        fp.write(fieldsep.join(headers) + "\n")
        for strand in strands:
            fp.write(fieldsep.join(strand[field] for field in headers) + "\n")
