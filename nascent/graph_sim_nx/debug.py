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

from pprint import pprint


def print_debug(*args, origin=None, **kwargs):
    print(*args, **kwargs)

def mute(*args, **kwargs):
    pass

# do_print = False
do_print = True

if do_print:
    pprintd = pprint
    printd = print_debug
else:
    pprintd = mute
    printd = mute
