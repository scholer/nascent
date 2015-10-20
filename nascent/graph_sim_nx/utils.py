# -*- coding: utf-8 -*-
##    Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.





def sequential_number_generator(start=0, exclude=None):
    uid = start
    while True:
        if not (exclude and uid in exclude):
            # if we have a set of excluded ids and uid is in that set, do not return uid, continue to next number
            yield uid
        uid += 1

sequential_uuid_gen = sequential_number_generator()


def touch(filepath):
    with open(filepath, 'a'):
        os.utime(filepath)


def find_unique_fn(filename):
    if not os.path.exists(filename):
        return filename
    regex_prog = re.compile(r"(.+?)(\d*)")
    fnroot, ext = os.path.splitext(filename)
    match = regex_prog.fullmatch(fnroot) # python 3.4+ only. Appends ^ before and $ after the regex.
    fnbase, num = match.groups()
    num = int(num) if num else 0
    fnfmt = "%s_%s%s"
    while os.path.exists(os.path.join(filename)):
        num += 1
        filename = fnfmt % (fnbase, num, ext)
    return filename
