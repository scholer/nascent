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

from __future__ import absolute_import, print_function, division
import os
import pdb
from collections import Counter, defaultdict
from itertools import chain
import yaml
import networkx as nx
try:
    import msgpack
except ImportError as e:
    print("Error importing msgpack library: %s" % e)
    print(" - msgpack output not available to stats writer.")


# from nascent.graph_sim_nx.reactionmgr import ReactionAttrs
from .nx_utils import draw_graph_and_save, layout_graph


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


def load_complex_state_count(fn):
    with open(fn) as fp:
        data = list(msgpack.Unpacker(fp, encoding='utf-8'))
    return data


class StatsWriter():
    """
    TODO: Implement support for msgpack format.
    """

    def __init__(self, sysmgr, simulator, config):

        self.sysmgr = sysmgr
        self.simulator = simulator
        self.config = config
        self.open_files = []
        self.monitored_strands = []   # Add strands to monitor throughout the simulation.

        self.stats_field_sep = config.get('stats_field_sep', '\t')
        self.stats_field_sep2 = config.get('stats_field_sep2', ', ')

        # Use standard fields,
        default_totals_fields = [
            'tau', 'system_time', 'temperature',
            'N_domains', 'n_hybridizable_domains', 'n_hybridized_domains',
            'N_strands', 'n_partially_hybridized_strands', 'n_fully_hybridized_strands',
            'n_stacked_ends',
            'n_complexes'
        ]
        self.stats_total_fields = config.get('stats_total_fields', default_totals_fields)
        # Or use a custom line format:
        # self.stats_total_fmt = config.get('stats_total_fmt')

        self.stats_per_domain_sysfields = config.get('stats_per_domain_sysfields',
                                                     ['tau', 'system_time', 'temperature'])
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
                                                     ['tau', 'system_time', 'temperature'])


        ## General "total" stats (for each step, file is kept open)
        self.stats_total_file = config.get('stats_total_file')
        if self.stats_total_file:
            print("Writing stats_total_file to file:", self.stats_total_file)
            self.stats_total_file = open(self.stats_total_file, 'w')
            header = self.stats_field_sep.join(self.stats_total_fields)
            self.stats_total_file.write(header+"\n")

        ## Per-domain stats (for each step, file is kept open)
        self.stats_per_domain_file = config.get('stats_per_domain_file')
        if self.stats_per_domain_file:
            print("Writing stats_per_domain_file to file:", self.stats_per_domain_file)
            self.stats_per_domain_file = fh = open(self.stats_per_domain_file, 'w')
            self.open_files.append(fh)
            # Has two headers, the first is for the "outer" fields: [sys-specs, domA, domB, ...]
            # The second header line is for the "inner" fields: [tau, system_time, ...], [n_hybridized, ...]
            header1 = self.stats_field_sep.join(["sysfields"] + self.stats_per_domain_species)
            header2 = self.stats_field_sep2.join(self.stats_per_domain_sysfields +
                                                 self.stats_per_domain_fields*len(self.stats_per_domain_species))
            self.stats_per_domain_file.write(header1+"\n"+header2+"\n")

        ## Per-strand stats (for each step)
        self.stats_per_strand_file = config.get('stats_per_strand_file')
        if self.stats_per_strand_file:
            print("Writing stats_per_strand_file to file:", self.stats_per_strand_file)
            self.stats_per_strand_file = fh = open(self.stats_per_strand_file, 'w')
            self.open_files.append(fh)
            # Has two headers, the first is for the "outer" fields: [sys-specs, domA, domB, ...]
            # The second header line is for the "inner" fields: [tau, system_time, ...], [n_hybridized, ...]
            header1 = self.stats_field_sep.join(["sysfields"] + self.stats_per_strand_species)
            header2 = self.stats_field_sep2.join(self.stats_per_strand_sysfields +
                                                 self.stats_per_strand_fields*len(self.stats_per_strand_species))
            self.stats_per_strand_file.write(header1+"\n"+header2+"\n")


        ## Collect complex_state_count (for each step) - using msgpack format
        self.complex_state_count_file = config.get('stats_complex_state_count_file')
        if self.complex_state_count_file:
            self.complex_state_count_file = open(self.complex_state_count_file, 'wb') # msgpack - open in binary mode
            self.open_files.append(self.complex_state_count_file)


        ## File to collect stats for monitored strands:
        self.stats_monitored_strands_file = config.get('stats_monitored_strands_file')
        # grouping-fields, then index field, then weight (tau), then values or secondary group fields (?)
        monitored_strands_fields = ("strand_uid, strand_name, system_time, tau, "
                                    "complex_uid, complex_state, N_domains, N_hybridized_domains").split(", ")
        if self.stats_monitored_strands_file:
            print("Writing monitored strands stats file:", self.stats_monitored_strands_file)
            self.stats_monitored_strands_file = fh = open(self.stats_monitored_strands_file, 'w')
            self.open_files.append(fh)
            header = self.stats_field_sep.join(monitored_strands_fields)
            self.stats_monitored_strands_file.write(header+"\n")


        ## File to collect complex states:
        self.stats_complex_state_file = config.get('stats_complex_state_file')
        complex_state_fields = [
            "complex_uid", "system_time", "tau", "reaction_str",
            #"cached_state_hash",   # Is always None
            "state_hash",           # Obtained with cmplx.state_fingerprint()
            #"cached_state_hash2",  # Same as state_hash
            "N_strands", "N_domains", "N_hybridized_pairs", "N_stacked_pairs",
            'total_dH', 'total_dS', 'volume_dH', 'volume_dS', 'shape_dH', 'shape_dS',
            'hybridization_dH', 'hybridization_dS', 'stacking_dH', 'stacking_dS'
        ]
        if self.stats_complex_state_file:
            print("Writing stats to stats_complex_state_file:", self.stats_complex_state_file)
            self.stats_complex_state_file = fh = open(self.stats_complex_state_file, 'w')
            self.open_files.append(fh)
            header = self.stats_field_sep.join(complex_state_fields)
            self.stats_complex_state_file.write(header+"\n")


        #### Post-simulation stats: ####

        ## Collected to a single file using msgpack, each entry is a single end-of-simulation dict.
        self.stats_post_simulation_file = config.get('stats_post_simulation_file')

        ## Reaction graph
        self.reaction_graph_output_directory = config.get('reaction_graph_output_directory')
        self.reaction_graph_output_fnfmt = config.get('reaction_graph_output_fnfmt')
        self.reaction_graph_output_formats = config.get('reaction_graph_output_formats')
        if isinstance(self.reaction_graph_output_formats, str):
            # Ensure that is is a list/tuple:
            self.reaction_graph_output_formats = [self.reaction_graph_output_formats]


    def close_all(self):
        """ Explicitly close all open files. """
        for fh in self.open_files:
            fh.close()


    def write_stats(self, tau, reaction_attr=None, result=None):
        """ Write all stats to all file. """
        ## TODO: Consider including reaction_spec_pair, reaction_attr (for the *next* or the *previous* reaction?)
        if self.stats_total_file is not None:
            self.write_total_stats(tau)
        if self.stats_per_domain_file is not None:
            self.write_per_domain_stats(tau)
        if self.stats_per_strand_file is not None:
            self.write_per_strand_stats(tau)
        if self.stats_monitored_strands_file is not None and self.monitored_strands:
            self.write_monitored_strands_stats(tau)
        # write_complex_state_stats currently works by writing state changes...
        if self.stats_complex_state_file and result is not None:
            self.write_complex_state_stats(tau, result, reaction_attr)



    def write_total_stats(self, tau, fields=None, init_stats=None):
        """
        Write system stats to file.
        """
        if fields is None:
            fields = self.stats_total_fields
        sysmgr = self.sysmgr
        # simulator = self.simulator
        if init_stats is None:
            stats_total = {}
        else:
            stats_total = init_stats.copy()

        ## Collect data to a single stats line dict:
        stats_total['tau'] = tau
        ## System/reaction/component manager stats:
        for attr in ('system_time', 'temperature', 'N_domains', 'N_strands'):
            stats_total[attr] = getattr(sysmgr, attr)
        for getter in ('n_hybridized_domains',
                       'n_hybridizable_domains',
                       'n_stacked_ends',
                       'n_partially_hybridized_strands',
                       'n_fully_hybridized_strands'):
            stats_total[getter] = getattr(sysmgr, getter)()
        stats_total['n_complexes'] = len(sysmgr.complexes)

        ## Write data to file:
        line = self.stats_field_sep.join(str(stats_total[field]) for field in fields)
        self.stats_total_file.write(line + "\n")


    def write_per_domain_stats(self, tau, sysfields=None, fields=None, species=None, init_stats=None):
        """
        Will first write constant fields, then domain species.
        system_time, temperature, N_domains, n_hybridized_domains, N_strands, ...
        """
        if fields is None:
            sysfields = self.stats_per_domain_sysfields
        if fields is None:
            fields = self.stats_per_domain_fields
        if species is None:
            species = self.stats_per_domain_species
        sysmgr = self.sysmgr
        if init_stats is None:
            stats_total = {}
        else:
            stats_total = init_stats.copy()

        ## Collect data to a single stats line dict:
        stats_total['tau'] = tau
        ## System/reaction/component manager stats:
        for attr in ('system_time', 'temperature', 'N_domains'):
            stats_total[attr] = getattr(sysmgr, attr)

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


    def write_per_strand_stats(self, tau, sysfields=None, fields=None, species=None, init_stats=None):
        """
        Will first write constant fields, then domain species.
        system_time, temperature, N_domains, n_hybridized_domains, N_strands, ...
        """
        if fields is None:
            sysfields = self.stats_per_domain_sysfields
        if fields is None:
            fields = self.stats_per_domain_fields
        if species is None:
            species = self.stats_per_strand_species
        sysmgr = self.sysmgr
        if init_stats is None:
            stats_total = {}
        else:
            stats_total = init_stats.copy()

        ## Collect data to a single stats line dict:
        stats_total['tau'] = tau
        ## System/reaction/component manager stats:
        for attr in ('system_time', 'temperature', 'N_domains'):
            stats_total[attr] = getattr(sysmgr, attr)

        ## Collect per-domain (specie) stats:
        strand_stats = {}
        for name, strands in sysmgr.strands_by_name.items():
            strand_stats[name] = {}
            strand_stats[name]['n_total'] = len(strands)
            strand_stats[name]['n_hybridized'] = sum(1 for strand in strands if strand.is_hybridized())
            strand_stats[name]['n_fully_hybridized'] = sum(1 for strand in strands if strand.is_fully_hybridized())

        ## Write data to file:
        # tau, system_time, temperature \t
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



    def write_monitored_strands_stats(self, tau):
        """
        Args:
            :result:    The result dict with case and changed complexes from react_and_process method.

        Will, for each monitored strand, write a line with:
            strand_uid, strand_name, system_time, tau, complex_uid, complex_state, N_domains, N_hybridized_domains
        This should be easy. The tricky part is parsing the result...
        You will have to create a table *for each complex* before you start calculating the duration
        of each complex state.

        To read the data, simply use something like:
            pandas.read_table(filename)
        """
        system_time = self.sysmgr.system_time
        for strand in self.monitored_strands:
            cmplx = strand.complex
            cuid, cstate = (-1, 0) if cmplx is None else (cmplx.cuid, cmplx.state_fingerprint())
            line = self.stats_field_sep.join((#"%s" % val for val in (
                "%s" % strand.suid, strand.name,
                "%0.05f" % system_time, "%0.04e" % tau, #reaction_str,
                "%s" % cuid, "%s" % cstate,
                "%s" % len(strand.domains), "%s" % sum(1 for d in strand.domains if d.partner is not None)
            ))
            self.stats_monitored_strands_file.write(line+"\n")


    def write_complex_state_stats(self, tau, result, reaction_attr):
        """
        Args:
            :result:    The result dict with case and changed complexes from react_and_process method.

        Will, for each complex, write a line with:
            <complex uuid>, system_time, N_strands, N_hybridized_domains, N_stacked_domains,
                and all fields from complex.energies_dHdS
        This should be easy. The tricky part is parsing the result...
        You will have to create a table *for each complex* before you start calculating the duration
        of each complex state.
        """
        system_time = self.sysmgr.system_time
        if result['changed_complexes'] is not None and result['new_complexes'] is not None:
            changed_complexes = result['changed_complexes'] + result['new_complexes']
        elif result['changed_complexes'] is not None:
            changed_complexes = result['changed_complexes']
        elif result['new_complexes'] is not None:
            changed_complexes = result['new_complexes']
        else:
            return
        reaction_attr_str = (reaction_attr.reaction_type + ("+" if reaction_attr.is_forming else "-")
                             + ("*" if reaction_attr.is_intra else " "))

        for cmplx in changed_complexes:
            ## in ReactionMgr we track N_domains_hybridized explicitly, but not for complexes.
            ## Using %s string interpolation seems to be the fastest way to produce strings:
            line = self.stats_field_sep.join((
                "%s" % cmplx.cuid, "%0.04e" % system_time, "%0.04e" % tau,
                reaction_attr_str,
                # cmplx._state_fingerprint,
                "%s" % cmplx.state_fingerprint(),
                # cmplx._state_fingerprint,
                "%s" % len(cmplx.strands), "%s" % len(list(cmplx.domains())),
                "%s" % len(cmplx.hybridized_pairs), "%s" % len(list(cmplx.stacked_pairs)),
                "%0.04f" % cmplx.energy_total_dHdS[0], "%0.04f" % cmplx.energy_total_dHdS[1],
                "%0.02f" % cmplx.energies_dHdS['volume'][0],        "%0.03f" % cmplx.energies_dHdS['volume'][1],
                "%0.02f" % cmplx.energies_dHdS['shape'][0],         "%0.03f" % cmplx.energies_dHdS['shape'][1],
                "%0.02f" % cmplx.energies_dHdS['hybridization'][0], "%0.03f" % cmplx.energies_dHdS['hybridization'][1],
                "%0.02f" % cmplx.energies_dHdS['stacking'][0],      "%0.03f" % cmplx.energies_dHdS['stacking'][1],
            ))
            self.stats_complex_state_file.write(line+"\n")


    def write_post_simulation_stats(self, fnpostfix=""):
        """
        Append a dict with stats:
            sysmgr_cache, reaction_attrs, reaction_throttle_cache, reaction_invocation_count
        to post_simulation_stats.yaml file.
        Note that I write the stats by appending the stats dict as:
            yaml.dump([stats])
        Doing this, I can append [stats] multiple times to the same file
        and still get a readable yaml-formatted list.
        (I typically append [stats] before and after a simulation.)
        """
        if self.stats_post_simulation_file is None:
            print("statsmgr.stats_post_simulation_file is None, cannot collect post simulation stats...")
            return
        sysmgr = self.sysmgr
        # systime = sysmgr.system_time
        stats = {}
        # remove defaultdict, etc:
        stats['reaction_throttle_cache'] = sysmgr.reaction_throttle_cache
        # ReactionMgr.cache has: domain_hybridization_energy, intracomplex_activity, stochastic_rate_constant,
        stats['sysmgr_cache'] = sysmgr.cache
        stats['reaction_attrs'] = sysmgr.reaction_attrs
        stats['reaction_spec_pairs'] = sysmgr.reaction_spec_pairs
        stats['reaction_invocation_count'] = sysmgr.reaction_invocation_count
        stats['possible_hybridization_reactions'] = sysmgr.possible_hybridization_reactions
        stats['possible_stacking_reactions'] = sysmgr.possible_stacking_reactions
        stats = simplify(stats)
        fn = self.stats_post_simulation_file.format(
            fnpostfix=fnpostfix,
            system_time=sysmgr.system_time, T=sysmgr.temperature)
        with open(fn, 'a') as fp:
            # dump yaml as list to make it easy to append.
            # otherwise, use msgpack format.
            before = fp.tell()
            yaml.dump([stats], fp)
            n_bytes = fp.tell() - before
        print("\nwrite_post_simulation_stats:", n_bytes, "bytes written to file", self.stats_post_simulation_file)


    def save_reaction_graph(self, **kwargs):
        """
        Save sysmgr.reaction_graph to file.
        See also ReactionMgr.save_reaction_graph (although this one actually saves all system graphs,
        not just the reaction graph.)
        And:
            This method saves in more formats..
            This method has different filename formatting, uses ReactionMgr.reaction_graph_output_fnfmt

        """

        # self.reaction_graph_output_directory = config.get('reaction_graph_output_directory')
        # self.reaction_graph_output_fnfmt = config.get('reaction_graph_output_fnfmt')
        # self.reaction_graph_output_formats = config.get('reaction_graph_output_formats')
        g = self.sysmgr.reaction_graph
        systime = self.sysmgr.system_time
        output_funcs = {method: getattr(nx, "write_"+method)
                        for method in ("yaml", "edgelist", "adjlist", "multiline_adjlist", "gexf", "pajek")}
        # output_funcs['png'] = draw_graph_and_save  # save reaction graph to png using e.g. graphviz
        if not os.path.exists(self.reaction_graph_output_directory):
            os.makedirs(self.reaction_graph_output_directory)
        if not os.path.isdir(self.reaction_graph_output_directory):
            print("Warning: output dir %s is not a director!" % self.reaction_graph_output_directory)
        for ext in self.reaction_graph_output_formats:
            path = os.path.join(self.reaction_graph_output_directory,
                                self.reaction_graph_output_fnfmt.format(ext=ext, systime=systime, **kwargs))
            # e.g. nx.write_gexf(g, path) or draw_graph_and_save(g, path)
            try:
                output_funcs[ext](g, path)
            except Exception as e:
                print("\nStatsManager: Error saving reaction graph using %r(%r, %r)" % (output_funcs[ext], g, path))
                print(" - exception type and msg:", type(e), e)




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
            tau, system_time, temperature, N_domains, n_hybridized_domains, ...
        Returns a list of stats-dicts

        """
        with open(fn) as fp:
            headers = next(fp).strip().split(self.stats_field_sep)
            rows = ([float(val) for val in line.strip().split(self.stats_field_sep)] for line in fp)
            stats = [dict(zip(headers, row)) for row in rows]
        return stats


    def load_per_domain_stats_file(self, fn):
        """
        First line is a super-header spec:
            # sysspecs, domainA, domainB, ... domain names.
        Then the fields header:
            tau, system_time, temperature \t n_total, n_hybridized, n_fully_hubridized \t (continue for next domain)
        The first three fields before the tab specifies the system state.
        The last three fields are "per-strand" stats field.
        Next line is values corresponding to the header.

        Returns a list of (system-dict, {domain-name: dom-stats-dict}) tuples.

        """
        with open(fn) as fp:
            first_row = next(fp).strip().split(self.stats_field_sep)
            names = first_row[1:]

            second_row = next(fp).strip().split(self.stats_field_sep)
            # [[tau, system_time, temperature], [n_total, n_hybridized, n_fully_hubridized], [...], ...]
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
            tau, system_time, temperature \t n_total, n_hybridized, n_fully_hubridized \t (continue for next domain)
        The first three fields before the tab specifies the system state.
        The last three fields are "per-strand" stats field.
        Next line is

        """
        with open(fn) as fp:
            first_row = next(fp).strip().split(self.stats_field_sep)
            names = first_row[1:]

            second_row = next(fp).strip().split(self.stats_field_sep)
            # [[tau, system_time, temperature], [n_total, n_hybridized, n_fully_hubridized], [...], ...]
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


    def load_complex_state_stats(self, fn):
        """
        returns
            stats, stats_by_cuid
        Where stats is a list of dicts with all entries

        and stats_by_cuid is the entries grouped by ComplexID.
        (This grouping could probably also be done by Pandas and a filter...)
        """
        with open(fn) as fp:
            headers = next(fp).strip().split(self.stats_field_sep)
            rows = ([float(val) for val in line.strip().split(self.stats_field_sep)] for line in fp)
            stats = [dict(zip(headers, row)) for row in rows]
            stats_by_cuid = defaultdict(list)
            # This could probably also be done by Pandas and a filter...
            for stat in stats:
                stats_by_cuid[stat['ComplexID']].append(stat)
        return stats, stats_by_cuid
