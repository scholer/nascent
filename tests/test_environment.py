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

# pylint: disable=C0103


"""

Module for testing the current runtime environment.
This includes:
* Disabling of random seeding for hashing.
    This requires setting the environment variable export PYTHONHASHSEED=0
    On Windows, do this by invoking ```set PYTHONHASHSEED=0``` on your command line before starting python.
    On Unix/OSX, invoke ```export PYTHONHASHSEED=0``` or directly while starting python with:
        ```PYTHONHASHSEED=0 python myscript.py```


Other things to consider:

    Windows:
    * If using the old command prompt, make sure unicode is enabled by:
        (a) Selecting a proper TrueType font, e.g. Lucida Console.
        (b) Invoking ```chcp 65001``` in your terminal.  (UNLESS you are using Python 2.7 or pypy)



"""


import sys
import os


def test_hash_seed_disabled():
    """
    Assert that random hash seeding has been disabled on the current system.
    Random seed hashing is enabled by default for python 3.3 and above, but
    can be disabled by setting the environment variable PYTHONHASHSEED to 0.

    This may vary depending on python version and implementation!
    hash seed randomization was not enabled by default until python 3.3
    For alternative hashing methods that should guarantee invariant results, see:
        * https://github.com/flier/pyfasthash

    """
    if sys.version_info >= (3, 3):
        assert os.environ['PYTHONHASHSEED'] == '0'
        assert hash("abc") == 4596069200710135518
    else:
        # Checked with CPython 2.6, 2.7,
        assert hash("abc") == -1600925533
        if not os.environ.get('PYTHONHASHSEED') == '0':
            print("WARNING: PYTHONHASHSEED environment variable not set, "
                  "but that is OK because you are also using an old python version.")


if __name__ == "__main__":
    ## Test 2:
    test_hash_seed_disabled()
    print(" - test_hash_seed_disabled OK")
