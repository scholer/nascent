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

Formats:
* text/csv
* json
* yaml
* msgpack


"""

import pdb
import yaml
from collections import Counter
try:
    import msgpack
except ImportError as e:
    print("Error importing msgpack library: %s" % e)
    print(" - msgpack output not available to stats writer.")


from nascent.graph_sim_nx.reactionmgr import ReactionAttrs


def simplify(data, _list=list, _set=set):
    if isinstance(data, dict):
        # cannot simplify if k is frozenset... :(
        return {simplify(k, _list=tuple, _set=frozenset): simplify(v, _list=_list, _set=_set)
                         for k, v in data.items()}
    if isinstance(data, (tuple, list)):
        return _list(simplify(v, _list=_list, _set=_set) for v in data)
    if isinstance(data, (set, frozenset)):
        return _set(simplify(v, _list=tuple, _set=frozenset) for v in data)
    return str(data)


class StatsWriter():
    """
    TODO: Implement support for msgpack format.
    """

    def __init__(self, sysmgr, simulator, config):

        self.sysmgr = sysmgr
        self.simulator = simulator
        self.config = config
        self.open_files = []


        self.stats_field_sep = config.get('stats_field_sep', '\t')
        self.stats_field_sep2 = config.get('stats_field_sep2', ', ')

        # Use standard fields,
        self.stats_total_fields = config.get('stats_total_fields',
            ['tau', 'sim_system_time', 'temperature',
             'N_domains', 'n_hybridizable_domains', 'n_hybridized_domains',
             'N_strands', 'n_partially_hybridized_strands', 'n_fully_hybridized_strands',
             'n_stacked_ends',
             'n_complexes'])
        # Or use a custom line format:
        # self.stats_total_fmt = config.get('stats_total_fmt')

        self.stats_per_domain_sysfields = config.get('stats_per_domain_sysfields',
                                                     ['tau', 'sim_system_time', 'temperature'])
        self.stats_per_domain_fields = config.get('stats_per_domain_fields',
                                                  ['n_total', 'n_hybridized'])
        self.stats_per_domain_species = config.get('stats_per_domain_species',
                                                   list(self.sysmgr.domains_by_name.keys()))
        if self.stats_per_domain_species == 'hybridizable':
            self.stats_per_domain_species = list(self.sysmgr.domain_pairs.keys())


        self.stats_per_strand_species = config.get('stats_per_strand_species',
                                                   list(sysmgr.strands_by_name.keys()))
        self.stats_per_strand_fields = config.get('stats_per_strand_fields',
                                                  ['n_total', 'n_hybridized', 'n_fully_hybridized'])
        self.stats_per_strand_sysfields = config.get('stats_per_strand_sysfields',
                                                     ['tau', 'sim_system_time', 'temperature'])


        self.stats_total_file = config.get('stats_total_file')
        if self.stats_total_file:
            print("Writing stats_total_file to file:", self.stats_total_file)
            self.stats_total_file = open(self.stats_total_file, 'w')
            header = self.stats_field_sep.join(self.stats_total_fields)
            self.stats_total_file.write(header+"\n")

        self.stats_per_domain_file = config.get('stats_per_domain_file')
        if self.stats_per_domain_file:
            print("Writing stats_per_domain_file to file:", self.stats_per_domain_file)
            self.stats_per_domain_file = fh = open(self.stats_per_domain_file, 'w')
            self.open_files.append(fh)
            # Has two headers, the first is for the "outer" fields: [sys-specs, domA, domB, ...]
            # The second header line is for the "inner" fields: [tau, sim_system_time, ...], [n_hybridized, ...]
            header1 = self.stats_field_sep.join(["sysfields"] + self.stats_per_domain_species)
            header2 = self.stats_field_sep2.join(self.stats_per_domain_sysfields +
                                                 self.stats_per_domain_fields*len(self.stats_per_domain_species))
            self.stats_per_domain_file.write(header1+"\n"+header2+"\n")

        self.stats_per_strand_file = config.get('stats_per_strand_file')
        if self.stats_per_strand_file:
            print("Writing stats_per_strand_file to file:", self.stats_per_strand_file)
            self.stats_per_strand_file = fh = open(self.stats_per_strand_file, 'w')
            self.open_files.append(fh)
            # Has two headers, the first is for the "outer" fields: [sys-specs, domA, domB, ...]
            # The second header line is for the "inner" fields: [tau, sim_system_time, ...], [n_hybridized, ...]
            header1 = self.stats_field_sep.join(["sysfields"] + self.stats_per_strand_species)
            header2 = self.stats_field_sep2.join(self.stats_per_strand_sysfields +
                                                 self.stats_per_strand_fields*len(self.stats_per_strand_species))
            self.stats_per_strand_file.write(header1+"\n"+header2+"\n")

        self.stats_post_simulation_file = config.get('stats_post_simulation_file')

        self.complex_state_count_file = config.get('complex_state_count_file')
        if self.complex_state_count_file:
            self.complex_state_count_file = open(self.complex_state_count_file, 'wb') # msgpack - open in binary mode
            self.open_files.append(self.complex_state_count_file)


    def close_all(self):
        for fh in self.open_files:
            fh.close()


    def write_stats(self, tau, sys_time):
        """ Write all stats to all file. """
        if self.stats_total_file is not None:
            self.write_total_stats(tau, sys_time)
        if self.stats_per_domain_file is not None:
            self.write_per_domain_stats(tau, sys_time)
        if self.stats_per_strand_file is not None:
            self.write_per_strand_stats(tau, sys_time)
        # pdb.set_trace()


    def write_total_stats(self, tau, sys_time=None, fields=None, init_stats=None):
        """
        Write system stats to file.
        """
        if fields is None:
            fields = self.stats_total_fields
        sysmgr = self.sysmgr
        simulator = self.simulator
        if init_stats is None:
            stats_total = {}
        else:
            stats_total = init_stats.copy()

        ## Collect data to a single stats line dict:
        stats_total['tau'] = tau
        ## System/reaction/component manager stats:
        for attr in ('temperature', 'N_domains', 'N_strands'):
            stats_total[attr] = getattr(sysmgr, attr)
        for getter in ('n_hybridized_domains',
                       'n_hybridizable_domains',
                       'n_stacked_ends',
                       'n_partially_hybridized_strands',
                       'n_fully_hybridized_strands'):
            stats_total[getter] = getattr(sysmgr, getter)()
        stats_total['n_complexes'] = len(sysmgr.complexes)
        ## Simulator stats:
        if simulator is not None:
            for attr in ('sim_system_time', ):
                stats_total[attr] = getattr(simulator, attr)
        if sys_time is not None:
            stats_total['sim_system_time'] = sys_time

        ## Write data to file:
        line = self.stats_field_sep.join(str(stats_total[field]) for field in fields)
        self.stats_total_file.write(line + "\n")


    def write_per_domain_stats(self, tau, sys_time=None, sysfields=None, fields=None, species=None, init_stats=None):
        """
        Will first write constant fields, then domain species.
        sim_system_time, temperature, N_domains, n_hybridized_domains, N_strands, ...
        """
        if fields is None:
            sysfields = self.stats_per_domain_sysfields
        if fields is None:
            fields = self.stats_per_domain_fields
        if species is None:
            species = self.stats_per_domain_species
        sysmgr = self.sysmgr
        simulator = self.simulator
        if init_stats is None:
            stats_total = {}
        else:
            stats_total = init_stats.copy()

        ## Collect data to a single stats line dict:
        stats_total['tau'] = tau
        ## System/reaction/component manager stats:
        for attr in ('temperature', 'N_domains'):
            stats_total[attr] = getattr(sysmgr, attr)
        # for getter in ('n_hybridized_domains',
        #                'n_hybridizable_domains',
        #                'n_partially_hybridized_strands',
        #                'n_fully_hybridized_strands'):
        #     stats_total[getter] = sysmgr.getter()
        #stats_total['n_complexes'] = len(sysmgr.complexes)
        ## Simulator stats:
        if simulator is not None:
            for attr in ('sim_system_time', ):
                stats_total[attr] = getattr(simulator, attr)
        if sys_time is not None:
            stats_total['sim_system_time'] = sys_time

        ## Collect per-domain (specie) stats:
        domain_stats = {}
        for name, domains in sysmgr.domains_by_name.items():
            domain_stats[name] = {}
            domain_stats[name]['n_total'] = len(domains)
            domain_stats[name]['n_hybridized'] = sum(1 for domain in domains if domain.partner is not None)

        ## Write data to file:
        line = self.stats_field_sep.join(
            [self.stats_field_sep2.join([str(stats_total[field]) for field in sysfields])]+
            [self.stats_field_sep2.join([str(domain_stats[name][field]) for field in fields])
            for name in self.stats_per_domain_species])

        self.stats_per_domain_file.write(line + "\n")


    def write_per_strand_stats(self, tau, sys_time=None, sysfields=None, fields=None, species=None, init_stats=None):
        """
        Will first write constant fields, then domain species.
        sim_system_time, temperature, N_domains, n_hybridized_domains, N_strands, ...
        """
        if fields is None:
            sysfields = self.stats_per_domain_sysfields
        if fields is None:
            fields = self.stats_per_domain_fields
        if species is None:
            species = self.stats_per_strand_species
        sysmgr = self.sysmgr
        simulator = self.simulator
        if init_stats is None:
            stats_total = {}
        else:
            stats_total = init_stats.copy()

        ## Collect data to a single stats line dict:
        stats_total['tau'] = tau
        ## System/reaction/component manager stats:
        for attr in ('temperature', 'N_domains'):
            stats_total[attr] = getattr(sysmgr, attr)

        ## Simulator stats:
        if simulator is not None:
            for attr in ('sim_system_time', ):
                stats_total[attr] = getattr(simulator, attr)
        if sys_time is not None:
            stats_total['sim_system_time'] = sys_time

        ## Collect per-domain (specie) stats:
        strand_stats = {}
        for name, strands in sysmgr.strands_by_name.items():
            strand_stats[name] = {}
            strand_stats[name]['n_total'] = len(strands)
            strand_stats[name]['n_hybridized'] = sum(1 for strand in strands if strand.is_hybridized())
            strand_stats[name]['n_fully_hybridized'] = sum(1 for strand in strands if strand.is_fully_hybridized())

        ## Write data to file:
        # tau, sim_system_time, temperature \t
        line = self.stats_field_sep.join(
            [self.stats_field_sep2.join([str(stats_total[field]) for field in sysfields])]+
            [self.stats_field_sep2.join([str(strand_stats[name][field]) for field in fields])
             for name in self.stats_per_strand_species])

        self.stats_per_strand_file.write(line + "\n")
        #pdb.set_trace()


    def write_complex_state_count(self, ):
        """
        What to record?
        For each complex:
        - N_strands, N_domains, N_domains_hybridized
        - Complex state
        Globally:
        - Complex state encounters (only register changes).
        - Complex states at time t.
        """
        sysmgr = self.sysmgr
        complex_state_count = dict(Counter([c.state_fingerprint() for c in sysmgr.complexes]))
        msgpack.pack(self.complex_state_count_file, complex_state_count)


    def load_complex_state_count(self, fn):
        with open(fn) as fp:
            data = list(msgpack.Unpacker(fp, encoding='utf-8'))
        return data


    def write_post_simulation_stats(self):
        if self.stats_post_simulation_file is None:
            print("statsmgr.stats_post_simulation_file is None, cannot collect post simulation stats...")
            return
        sysmgr = self.sysmgr
        simulator = self.simulator
        stats = {}
        # remove defaultdict, etc:
        stats['reaction_throttle_cache'] = sysmgr.reaction_throttle_cache
        stats['sysmgr_cache'] = sysmgr.cache
        stats['reaction_attrs'] = sysmgr.reaction_attrs
        stats['reaction_invocation_count'] = sysmgr.reaction_invocation_count
        stats = simplify(stats)
        with open(self.stats_post_simulation_file, 'a') as fp:
            # dump yaml as list to make it easy to append.
            # otherwise, use msgpack format.
            before = fp.tell()
            yaml.dump([stats], fp)
            n_bytes = fp.tell() - before
        print("\nwrite_post_simulation_stats:", n_bytes, "bytes written to file", self.stats_post_simulation_file)



class StatsReader():

    def __init__(self, config):
        if config is None:
            config = {}
        self.config = config

        self.stats_field_sep = config.get('stats_field_sep', '\t')
        self.stats_field_sep2 = config.get('stats_field_sep2', ', ')


    def load_total_stats_file(self, fn):
        """
        Then the fields header:
            tau, sim_system_time, temperature, N_domains, n_hybridized_domains, ...
        Returns a list of stats-dicts

        """
        with open(fn) as fp:
            headers = next(fp).strip().split(self.stats_field_sep)
            rows = ([float(val) for val in line.strip().split(self.stats_field_sep)] for line in fp)
            # each row is like: [[0.2, 1.2, 330], [10, 5, 2], [...], ...]
            #sysstats, *domainstats =
            stats = [dict(zip(headers, row)) for row in rows]
        return stats


    def load_per_domain_stats_file(self, fn):
        """
        First line is a super-header spec:
            # sysspecs, domainA, domainB, ... domain names.
        Then the fields header:
            tau, sim_system_time, temperature \t n_total, n_hybridized, n_fully_hubridized \t (continue for next domain)
        The first three fields before the tab specifies the system state.
        The last three fields are "per-strand" stats field.
        Next line is values corresponding to the header.

        Returns a list of (system-dict, {domain-name: dom-stats-dict}) tuples.

        """
        with open(fn) as fp:
            first_row = next(fp).strip().split(self.stats_field_sep)
            names = first_row[1:]

            second_row = next(fp).strip().split(self.stats_field_sep)
            # [[tau, sim_system_time, temperature], [n_total, n_hybridized, n_fully_hubridized], [...], ...]
            headers = [substr.split(self.stats_field_sep2) for substr in second_row]

            rows = ([[float(val) for val in substr.split(self.stats_field_sep2)]
                     for substr in line.strip().split(self.stats_field_sep)]
                    for line in fp)
            # each row is like: [[0.2, 1.2, 330], [10, 5, 2], [...], ...]
            #sysstats, *domainstats =
            stats = [(dict(zip(headers[0], row[0])),
                      {name: dict(zip(headers[1:], domvalues)) for name, domvalues in zip(names, row[1:])})
                     for row in rows]
        return stats


    def load_per_strands_stats_file(self, fn):
        """
        First line is a super-header spec:
            # sysspecs, domainA, domainB, ... domain names.
        Then the fields header:
            tau, sim_system_time, temperature \t n_total, n_hybridized, n_fully_hubridized \t (continue for next domain)
        The first three fields before the tab specifies the system state.
        The last three fields are "per-strand" stats field.
        Next line is

        """
        with open(fn) as fp:
            first_row = next(fp).strip().split(self.stats_field_sep)
            names = first_row[1:]

            second_row = next(fp).strip().split(self.stats_field_sep)
            # [[tau, sim_system_time, temperature], [n_total, n_hybridized, n_fully_hubridized], [...], ...]
            headers = [substr.split(self.stats_field_sep2) for substr in second_row]

            rows = ([[float(val) for val in substr.split(self.stats_field_sep2)]
                     for substr in line.strip().split(self.stats_field_sep)]
                    for line in fp)
            # each row is like: [[0.2, 1.2, 330], [10, 5, 2], [...], ...]
            #sysstats, *domainstats =
            stats = [(dict(zip(headers[0], row[0])),
                      {name: dict(zip(headers[1:], domvalues)) for name, domvalues in zip(names, row[1:])})
                     for row in rows]
        return stats
