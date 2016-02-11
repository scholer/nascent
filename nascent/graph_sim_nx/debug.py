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

#pylint: disable=C0103,C0111,W0613


from __future__ import absolute_import, print_function, division
from pprint import pprint

do_print = False


def print_debug(*args, **kwargs):
    """ Change module-level do_print variable to toggle behaviour. """
    if 'origin' in kwargs:
        del kwargs['origin']
    if do_print:
        print(*args, **kwargs)

def pprint_debug(*args, **kwargs):
    if do_print:
        pprint(*args, **kwargs)


pprintd = pprint_debug
printd = print_debug
