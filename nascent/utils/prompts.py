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


def prompt_yes_no(prompt, default=None):
    """ Prompt user for a yes/no answer, defaulting to :default:.
    If default is None, this function will continue to ask until a clear yes or no has been given."""
    answer = ""
    while answer not in ('y', 'n'):
        answer = input(prompt).strip().lower()
        if not answer:
            if default:
                answer = default
        else:
            answer = answer[0].lower()
        if answer not in ('y', 'n'):
            print("Answer must start with Y/y or N/n. Please try again.")
    return answer


def start_interact(local=None, readfunc=None, banner=None):
    try:
        import rlcompleter
        import readline
        readline.set_completer(rlcompleter.Completer(locals()).complete)
        readline.parse_and_bind("tab: complete")
    except ImportError:
        pass
    from importlib import reload
    import pdb
    import code #from code import interact
    if local is None:
        local = locals()
    else:
        local['reload'] = reload
        local['pdb'] = pdb
    code.interact(local=local, readfunc=readfunc, banner=banner)
