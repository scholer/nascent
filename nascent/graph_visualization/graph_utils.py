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


def directed_for_all_edges(edges):
    """
    Returns a 'directed' argument for edges.
    If all edges have the same directed value, only this value is used.
    """
    directed = [edge.get('directed') for edge in edges]
    if len(set(directed)) < 2:
        # If all edges have the same 'directed' value, we can optimize by only providing one value:
        directed = directed[0]
    return directed
