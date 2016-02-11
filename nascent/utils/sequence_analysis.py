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


Possible refs:
* https://github.com/joaks1/partition_summary
* https://github.com/htailor/simulation_dna, https://github.com/htailor/new_simulation_dna

Consider using Directed Acyclic Word Graphs (DAWG):
* https://github.com/kmike/DAWG

"""

from __future__ import absolute_import, print_function
import os
import sys
#import difflib
from difflib import SequenceMatcher
import networkx as nx
from networkx.drawing import nx_agraph
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot
import webbrowser
from pprint import pprint
from collections import defaultdict, deque
import string
import math
import pdb
import subprocess
import itertools



def domains_from_sequences(group1, group2=None, minoverlap=8, maxoverlap=32, match_self=False, promptlimit=3):
    """
    :group1: list of dicts with keys {'name', 'seq', 'rcompl'}.
    :promptlimit: When splitting up primary fragments, prompt if the overlap is equal-to-or-less than this limit.
    Maybe just blasting them is easier?

    So many alternative strategies:
    * Binary masks
    * Graph - binary heap
    * Clasps
    * All-fragments map
    * Largest-fragment map
    * Find-all-first vs find-longest-then-check-reduce

    How about this:
    For each strand we make a structure where the top node is the strand it self.
    Child nodes are fragments that match against other strands.         . AGT
    We do not check all fragment permutation, which would be:    AGTC -<
                                                                        ` GTC
    Instead, we only add fragments that have *additional matches* compared to
    the parent, or if we have *overlapping matches* that share parent as base,
    e.g. in the tree above, the fragments "AGT" and "GTC" share "AGTC" as parent.
    Note that while this is still a heap tree, it is not binary.
    (I note this because "binary" is sometimes assumed when talking about heaps.)
    What about:               .- ABCDEFG
                ABCDEFGHIJK -<
                              `- FGHIJK
    And then we have another match, "DEFGH"? Just add it to the parent and analyse
    afterwards. To resolve these shared overlaps, we might have to create bipartite
    graphs (considering only the overlapping nodes).
    Edit: Maybe have just a single graph?
                                     .- XYWV
                       XYWVABCDEFG -<             .- ABCDE
                                     `- ABCDEFG -<
     NABCDEFGHIJKLM -- ABCDEFGHIJK -<             `- FG
                                     `- FGHIJK  -<
                       OPQRSFGHIJK -<             `- HIJK
                                     `- OPQRS
    Each node has name=seq and attributes:
    * length,
    Each edge has attributes:
    * Offset, rcompl (true/false),
    Maybe first build a 1-level "parent, overlaps" bipartite tree,
    then process this, doing sub-division and making "agg" nodes?
    """
    answer = ""
    # Basic, 2-level bipartite directed graphs:
    strands = group1 + (group2 if group2 is not None else [])

    ptoc, ctop, rcgraph, overlaps = overlap_graphs(group1, group2=group2, ptoc=None, ctop=None, rcgraph=None,
                                                   minoverlap=minoverlap, maxoverlap=maxoverlap,
                                                   match_self=match_self)
    ptoc.graph['graph'] = {'rankdir': 'LR'}
    fn = os.path.join(os.path.abspath('.'), "seq_analysis_1_bipartite.dot")
    # make_graphviz_graph(fn, ptoc, rcgraph, add_same_rank=True, add_labels=True)

    # nx.write_dot(ptoc, fn)
    # fn = os.path.join(os.path.abspath('.'), "seq_analysis_bipartite_1.agraph.dot")
    # nx_agraph.write_dot(ptoc, fn)
    # fn = os.path.join(os.path.abspath('.'), "seq_analysis_bipartite_1.png")
    # pos = nx.pydot_layout(ptoc, prog='dot')
    # nx.draw_networkx(ptoc, pos=pos, with_labels=True)
    # pyplot.savefig(fn)
    # webbrowser.open(fn)

    ## Annotate fragment competitors:
    #answer = input("\nPress enter to annotate fragment base competitors...")
    if answer and answer[0].lower() == 'q': return
    fcompg = annotate_fragment_competitors(ptoc, ctop, rcgraph, strands, fcompg=None)
    fn = os.path.join(os.path.abspath('.'), "seq_analysis_2-3_competitors.dot")
    # make_graphviz_graph(fn, fcompg, ptoc, rcgraph, add_same_rank=True, add_labels=True)


    ### Re-partition fragments overlapping the same strand with an amount less than a specified limit:
    # answer = input("\nPress enter to re-partition fragments with minor overlaps or type 'q' to abort...")
    if answer and answer[0].lower() == 'q': return
    find_compl(ptoc, rcgraph) # find complementary nodes in ptoc and add edges between these in rcgraph.
    repartition_overlapping_fragments(ptoc, ctop, rcgraph, strands, overlap_limit=4)
    make_graphviz_labels(ptoc)
    fn = os.path.join(os.path.abspath('.'), "seq_analysis_2_repartitioned.dot")
    make_graphviz_graph(fn, ptoc, rcgraph, add_same_rank=True, add_labels=True)


    ## Remove short fragments (caused by re-partitioning):
    answer = input("\nPress enter to remove short fragments reduced by re-partitioning...")
    if answer and answer[0].lower() == 'q': return
    remove_short_fragments(ptoc, ctop, rcgraph, limit=5)
    fn = os.path.join(os.path.abspath('.'), "seq_analysis_2-2_short_removed.dot")
    make_graphviz_graph(fn, ptoc, rcgraph, add_same_rank=True, add_labels=True)


    ## Annotate fragment competitors:
    answer = input("\nPress enter to annotate fragment base competitors...")
    if answer and answer[0].lower() == 'q': return
    fcompg = annotate_fragment_competitors(ptoc, ctop, rcgraph, strands, fcompg=None)
    fn = os.path.join(os.path.abspath('.'), "seq_analysis_2-3_competitors.dot")
    make_graphviz_graph(fn, fcompg, ptoc, rcgraph, add_same_rank=True, add_labels=True)

    ## Find overlapping fragments and split the overlaps to separate domains:
    answer = input("\nPress enter to split remaining overlapping fragments...")
    if answer and answer[0].lower() == 'q': return
    split_overlapping_fragments(ptoc, ctop, rcgraph, strands)
    fn = os.path.join(os.path.abspath('.'), "seq_analysis_3_subfragmented.dot")
    make_graphviz_graph(fn, ptoc, rcgraph, add_same_rank=True, add_labels=True)
    # nx.write_dot(ptoc, fn)

    ## Go over all newly-created sub-fragments and either revert the split to one side or another,
    ## or propagate the split to complement sequences.
    # for sf_seq in subfragments:
    #     #if len(ctop[sf_seq]) >
    #     if len(sf_seq) <= promptlimit:
    #         #print("Two or more fragments are overlapping on strand %s:\n  %s" % (strand['Name'], strand['seq']))
    #         print("The sub-fragment %s is shorter than the promptlimit %s." % (sf_seq, promptlimit))
    #         answer = input(" Do you want to try to reverse this? ")
    #         if answer and answer[0].lower() == 'y':
    #             print(" Which fragments do you want to split?")
    #             sense_fragments = list(ctop[sf_seq].keys())
    #             # sense_fragments = [(fseq, eattrs) for fseq, eattrs in ptoc[pseq].items() if eattrs['rcompl'] is False]
    #             for i, fseq in enumerate(sense_fragments, 1):
    #                 print(i, fseq)
    #             answer = input(" ")
    #             answer = [int(sel) for sel in answer.split(",")]
    #             remove_edges = [sense_fragments[sel] for sel in answer]
    #             for target in remove_edges:


    ## Add single/un-paired sequence nodes to graphs:
    answer = input("\nPress enter to add unpaired strand domains...")
    if answer and answer[0].lower() == 'q': return
    for strand in itertools.chain(group1, group2 if group2 is not None else []):
        pseq = strand['seq']
        # Only have to consider sense, not rcompl:
        fragments = sorted([(eattr['offset'], eattr['end']) for fseq, eattr in ptoc[pseq].items()
                            if not eattr['rcompl']])
        # if seq is rcompl, then offset < end! (because reversed)
        coverage = unify_all_regions(fragments)
        noncovered = invert_regions(coverage, [0, len(pseq)])
        for start, end in noncovered:
            nc_seq = pseq[start:end]
            ptoc.add_node(nc_seq)
            ctop.add_node(nc_seq)
            ptoc.add_edge(pseq, nc_seq, offset=start, end=end, rcompl=False)
            ctop.add_edge(nc_seq, pseq, offset=start, end=end, rcompl=False)


    # make_graphviz_labels(ptoc)
    # find_compl(ptoc, rcgraph) # find complementary nodes in ptoc and add edges between these in rcgraph.
    fn = os.path.join(os.path.abspath('.'), "seq_analysis_3_unpaired_domains.dot")
    print("\nWriting dot file to", fn)
    # joined = nx.DiGraph(ptoc)
    # joined.add_nodes_from(rcgraph.nodes())
    # joined.add_edges_from(rcgraph.edges())
    # nx.write_dot(joined, fn)
    # fn = os.path.join(os.path.abspath('.'), "seq_analysis_3_rcgraph.dot")
    # make_graphviz_graph(fn, rcgraph)
    for strands in [group for group in (group1, group2) if group is not None]:
        extract_domains(ptoc, strands, rcgraph=rcgraph, namegen=None)
    print("\nStrands after extracting domains:")
    pprint(strands)
    # make_graphviz_labels(ptoc)
    make_graphviz_graph(fn, ptoc, rcgraph, add_same_rank=True, add_labels=True)
    # pdb.set_trace()
    return ptoc, ctop, rcgraph



def overlap_graphs(group1, group2=None, ptoc=None, ctop=None, rcgraph=None,
                   minoverlap=8, maxoverlap=32, match_self=False, negative_idx=True):
    """
    For sequences in group1 vs group2, generate basic, 2-level bipartite directed graphs.
    :ptoc: DiGraph going from parents (full length strands) to child fragments.
    :ctop: DiGraph from chlild overlap fragments to full-length parent sequences.
    :overlaps: Set of overlap sequences.
    """

    if group2 is None:
        group2 = group1

    # overlap_trees = {} # Indexed by strand name => tree
    if ptoc is None:
        ptoc = nx.DiGraph() # directed graph with edges from parent to children.
    if ctop is None:
        ctop = nx.DiGraph() # DiGraph with edges from children to parents.
    if rcgraph is None:
        rcgraph = nx.DiGraph()
    fragments = set()
    revc = rcgraph

    for strand in group1:
        ptoc.add_node(strand['seq'], attr_dict=strand, rank=0)
        ctop.add_node(strand['seq'], attr_dict=strand, rank=0)
        if 'rcompl' not in strand:
            strand['rcompl'] = "".join(WC[b] for b in reversed(strand['seq']))
        for strand2 in group2:
            if strand2 == strand:
                continue
            if 'rcompl' not in strand2:
                strand2['rcompl'] = "".join(WC[b] for b in reversed(strand2['seq']))
            ptoc.add_node(strand2['seq'], attr_dict=strand2, rank=0)
            ctop.add_node(strand2['seq'], attr_dict=strand2, rank=0)
            for rcompl in (True, ): #(False, True):
                overlaps = find_overlaps(strand['seq'], strand2['rcompl'] if rcompl else strand2['seq'])
                # each match is a named tuple, Match
                for overlap in overlaps:
                    if overlap.size < minoverlap:
                        continue
                    fseq = strand['seq'][overlap.a:overlap.a+overlap.size]
                    if rcompl:
                        seqb = strand2['seq'][-1-overlap.b:-1-overlap.b-overlap.size:-1]
                        rseq = strand2['rcompl'][overlap.b:overlap.b+overlap.size]
                        assert fseq == rseq
                    else:
                        seqb = strand2['seq'][overlap.b:overlap.b+overlap.size]
                        assert fseq == seqb
                    fragments.add(fseq)
                    if fseq not in ptoc:
                        assert fseq not in ctop
                        ptoc.add_node(fseq, length=overlap.size, rank=1)
                        ctop.add_node(fseq, length=overlap.size, rank=1)
                        revc.add_node(fseq, length=overlap.size, rank=1)
                    ptoc.add_edge(strand['seq'], fseq, offset=overlap.a, end=overlap.a+overlap.size, rcompl=False)
                    ctop.add_edge(fseq, strand['seq'], offset=overlap.a, end=overlap.a+overlap.size, rcompl=False)
                    if negative_idx and rcompl:
                        offset, end = [-1 * (overlap.b + i -1) for i in (0, overlap.size)]
                    else:
                        offset, end = overlap.b, overlap.b+overlap.size
                    ptoc.add_edge(strand2['seq'], fseq, offset=offset, end=end, rcompl=rcompl)
                    ctop.add_edge(fseq, strand2['seq'], offset=offset, end=end, rcompl=rcompl)

    find_compl(ptoc, rcgraph)
    #find_compl(ctop, rcgraph)
    return ptoc, ctop, rcgraph, fragments



def repartition_overlapping_fragments(ptoc, ctop, rcgraph, strands, overlap_limit=4):

    # if strands is None:
    #     pseqs = [strand['seq'] for strand in strands]
    # else:
    #     pseqs = [node for node in ptoc.nodes() if ptoc.in_degree(node) == 0]
    # for pseq in pseqs:
    strands_by_seq = {s['seq']: s for s in strands}
    revc = rcgraph
    for strand in strands:
        pseq = strand['seq']
        overlapping_edges = has_overlap_shorter_than(ptoc[pseq].items(), overlap_limit=overlap_limit)
        while overlapping_edges:
            e1, e2, overlap_size = overlapping_edges
            assert edge_overlap_size(e1, e2)
            assert e1[1]['offset'] < e2[1]['offset']
            (f1seq, f1attr, e1attr), (f2seq, f2attr, e2attr) = [(e[0], ptoc.node[e[0]], e[1]) for e in (e1, e2)]
            e1_competition_score = 0
            e2_competition_score = 0
            print("\nTwo or more fragments are overlapping on strand %s:\n   %s" % (strand['Name'], strand['seq']))
            if e1[0] in rcgraph.node and len(rcgraph[e1[0]]) > 0:
                assert len(rcgraph[e1[0]]) == 1
                rseq = next(iter(rcgraph[e1[0]].keys()))
                rattr = ptoc.node[rseq]
                if 'base_competitors_count' in rattr:
                    e1_competition_score += sum(rattr['base_competitors_count'])/2/len(rseq)
                    print("   " + " "*e1[1]['offset'] + "".join(str(i) for i in rattr['base_competitors_count'][::-1]))
            if 'base_competitors_count' in f1attr:
                print("   " + " "*e1[1]['offset'] + "".join(str(i) for i in f1attr['base_competitors_count']))
                e1_competition_score += sum(f1attr['base_competitors_count'])/len(f1seq)
            # print("  " + "-"*e2[1]['offset'] + "".join('\\' for i in range(0, overlap_size+1)))
            for i, (fseq, eattrs) in enumerate((e1, e2), 1):
                print(("e%s " % i) + " "*eattrs['offset'] + fseq)
            if 'base_competitors_count' in f2attr:
                print("   " + "-"*e2[1]['offset'] + "".join(str(i) for i in f2attr['base_competitors_count']))
                e2_competition_score += sum(f2attr['base_competitors_count'])/len(f2seq)
            if e2[0] in rcgraph.node and len(rcgraph[e2[0]]) > 0:
                assert len(rcgraph[e2[0]]) == 1
                rseq = next(iter(rcgraph[e2[0]].keys()))
                rattr = ptoc.node[rseq]
                if 'base_competitors_count' in rattr:
                    e2_competition_score += sum(rattr['base_competitors_count'])/2/len(rseq)
                    print("   " + " "*e2[1]['offset'] + "".join(str(i) for i in rattr['base_competitors_count'][::-1]))
            print("   " + "-"*(e2[1]['offset']-1) + "".join('^' for i in range(0, overlap_size+1)))
            print("   " + "-"*(e2[1]['offset']-1) + "".join(str(i) for i in range(0, overlap_size+1)))
            print("The overlap is less than the specified overlap_limit %s." % overlap_limit)
            print("What fragment do you want to partition this fragment to?",
                  "(Type e1 or e2 or a number to split the overlap in the middle.",
                  "Type s to skip this for now, type 'b' to break to next strand.) ")
            if (e1_competition_score > 0 or e2_competition_score > 0) and e1_competition_score != e2_competition_score:
                answer = min(((e1_competition_score, 'e1'), (e2_competition_score, 'e2')))[1]
            else:
                answer = ""
            # pdb.set_trace()
            answer = input("Input any of 'e1' / 'e2' / 0...N / 's' / 'b' " +
                           ("[default: %s]" % answer if answer else "") + " : ") or answer
            while not answer:
                answer = input()
            if answer.lower()[0] == 'b':
                break
            elif answer.lower()[0] == 's':
                continue
            elif answer.lower()[0] == 'q':
                sys.exit(1)
            save, remove = [], []

            #e1_size = 0 if answer == 'e1' else (overlap_size if answer == 'e2' else int(answer))
            if answer == 'e1': # save e1
                e1_size = overlap_size
            elif answer == 'e2':
                e1_size = 0
            else:
                e1_size = int(answer)

            if e1_size == 0:
                save.append((e2, None))
                # If overlap_size is 5, cut e1 seq at [:-5]
                remove.append((e1, (None, -overlap_size), (e1[1]['offset'], e1[1]['end']-overlap_size)))
                # Recording redundant (edge, fseq_slice, pseq_slice) -- to assert that fseq[fslice] == pseq[pslice]
            elif e1_size == overlap_size:
                save.append((e1, None))
                # If overlap_size is 5, cut e2 seq at [5:]
                remove.append((e2, (overlap_size, None), (e2[1]['offset']+overlap_size, e2[1]['end'])))
            else:
                assert 0 < e1_size < overlap_size
                cutting_point = overlap_size - e1_size
                remove.append((e1, (None, -cutting_point), (e1[1]['offset'], e1[1]['end']-cutting_point)))
                remove.append((e2, (cutting_point, None), (e2[1]['offset']+cutting_point, e2[1]['end'])))

            for (fseq, eattr), fslice, pslice in remove:
                # pstart, pend = pslice
                fattr = ptoc.node[fseq]
                fslice, pslice = slice(*fslice), slice(*pslice)
                assert fseq[fslice] == pseq[pslice]
                fseq_new = fseq[fslice]
                ## Add new node and edge from pseq:
                nattrs = {'length': len(fseq_new),
                          'rank': 1}
                if 'base_competitors' in fattr:
                    nattrs['base_competitors'] = fattr['base_competitors'][fslice]
                if 'base_competitors_count' in fattr:
                    nattrs['base_competitors_count'] = fattr['base_competitors_count'][fslice]
                ## TODO: Instead of removing old node and adding new, just update the node+edge attributes.
                ## - use nx.relabel_nodes(G, mapping, copy=False)
                ptoc.add_node(fseq_new, attr_dict=nattrs)
                ctop.add_node(fseq_new, attr_dict=nattrs)
                revc.add_node(fseq_new, attr_dict=nattrs)
                # Add edge from parent strand to new/shortened fragment seq node:
                ptoc.add_edge(pseq, fseq_new, offset=pslice.start, end=pslice.stop, rcompl=eattr['rcompl'])
                ctop.add_edge(fseq_new, pseq, offset=pslice.start, end=pslice.stop, rcompl=eattr['rcompl'])
                ## Make sure to propagate the change to all originating parent strands:
                if len(ctop[fseq]) > 0:
                    for other_pseq, oeattr in ctop[fseq].items():
                        if oeattr['rcompl']:
                            # offset and end may be negative for rcomplement edges
                            offset, end = oeattr['offset']-(fslice.stop or 0), oeattr['end']-(fslice.start or 0)
                        else:
                            offset, end = oeattr['offset']+(fslice.start or 0), oeattr['end']+(fslice.stop or 0)
                        if other_pseq == pseq:
                            assert offset == pslice.start and end == pslice.stop
                        else:
                            ptoc.add_edge(other_pseq, fseq_new, offset=offset, end=end, rcompl=oeattr['rcompl'])
                            ctop.add_edge(fseq_new, other_pseq, offset=offset, end=end, rcompl=oeattr['rcompl'])
                ## If fseq has rcompl fragment, check if we should also propagate this:
                if len(rcgraph[fseq]) > 0:
                    assert len(rcgraph[fseq]) == 1
                    rseq = next(iter(rcgraph[fseq].keys()))
                    rattr = ptoc.node[rseq]
                    rseq_parents = list(ctop[rseq].items())
                    rfslice = slice((-fslice.stop if fslice.stop else None), (-fslice.start if fslice.start else None))
                    rseq_new = rseq[rfslice]
                    assert rseq_new == "".join(WC[b] for b in pseq[pslice][::-1])
                    print("Sequence %s, which is being re-partitioned, has complement %s from parents %s." %
                          (fseq, rseq, [rp[0]+(" (rev)" if rp[1]['rcompl'] else "") for rp in rseq_parents]))
                    propagate = input("Alter this complement sequence accordingly? [Y/n] ")
                    if not (propagate and propagate[0] == 'n'):  # default = "yes, alter complement accordingly."
                        # Removing node will also remove edge
                        ptoc.remove_node(rseq)
                        ctop.remove_node(rseq)
                        revc.remove_node(rseq)
                    ## Add reverse complement node
                    nattrs = {'length': len(rseq_new),
                              'rank': 1}
                    if 'base_competitors' in rattr:
                        nattrs['base_competitors'] = rattr['base_competitors'][rfslice]
                    if 'base_competitors_count' in rattr:
                        nattrs['base_competitors_count'] = rattr['base_competitors_count'][rfslice]
                    ptoc.add_node(rseq_new, attr_dict=nattrs)
                    ctop.add_node(rseq_new, attr_dict=nattrs)
                    revc.add_node(rseq_new, attr_dict=nattrs)
                    ## Connect new fseq and rseq in rcgraph:
                    revc.add_edge(fseq_new, rseq_new, rcompl=True)
                    revc.add_edge(rseq_new, fseq_new, rcompl=True)
                    ## Re-create rcomplement's edges to parent strands:
                    for rpseq, reattr in rseq_parents:
                        if reattr['rcompl']:
                            # offset and end may be negative for rcomplement edges
                            offset, end = reattr['offset']-(rfslice.stop or 0), reattr['end']-(rfslice.start or 0)
                        else:
                            offset, end = reattr['offset']+(rfslice.start or 0), reattr['end']+(rfslice.stop or 0)
                        ptoc.add_edge(rpseq, rseq_new, offset=offset, end=end, rcompl=reattr['rcompl'])
                        print("ptoc.add_edge(%s, %s, offset=%s, end=%s, rcompl=%s)" %
                              (rpseq, rseq_new, offset, end, reattr['rcompl']))
                        ctop.add_edge(rseq_new, rpseq, offset=offset, end=end, rcompl=reattr['rcompl'])
                # Remove old/original fragment seq node (and all connecting edges):
                ptoc.remove_node(fseq)
                ctop.remove_node(fseq)
                revc.remove_node(fseq)
            # Re-check for overlapping edges
            overlapping_edges = has_overlap_shorter_than(ptoc[pseq].items(), overlap_limit=overlap_limit)
            # (and go back to while-loop check)
        # end while-loop: no more overlapping strand edges (or user-aborted through 'b' input).
    print("Strand fragments overlap re-partition complete.")


def remove_short_fragments(ptoc, ctop, rcgraph, limit=5):
    """
    Remove fragments shorter than :limit: from the graphs.
    """
    revc = rcgraph
    for short_leaf in short_leaf_node_gen(ptoc, limit):
        print("Removing short fragment %s" % (short_leaf, ))
        ptoc.remove_node(short_leaf)
        ctop.remove_node(short_leaf)
        if short_leaf in rcgraph.node:
            if len(rcgraph[short_leaf]) > 0:
                assert len(rcgraph[short_leaf]) == 1
                rseq = list(rcgraph[short_leaf])[0]
                print(" - also removing short complementary fragment %s" % (rseq, ))
                ptoc.remove_node(rseq)
                ctop.remove_node(rseq)
                revc.remove_node(rseq)
            revc.remove_node(short_leaf)
        # short_leaf = get_short_leaf_node(ptoc, limit)


def leaf_node_gen(ptoc, include_unbound=False):
    nodes = ptoc.nodes()
    for node in nodes:
        if ptoc.out_degree(node) == 0 and (ptoc.in_degree(node) > 0 or include_unbound):
            yield node



def short_leaf_node_gen(ptoc, limit=5):
    """ Can be optimized as a single-run generator. """
    nodes = ptoc.nodes_iter()
    yield next(seq for seq in nodes
               if len(seq) < limit and ptoc.out_degree(seq) == 0 and ptoc.in_degree(seq) > 0)



def has_overlap_shorter_than(edges, overlap_limit=4):
    edges = [(fseq, eattr) for fseq, eattr in edges if eattr['rcompl'] is False]
    edges = sorted(edges, key=lambda e: e[1]['offset'])
    for e1, e2 in itertools.combinations(edges, 2):
        overlap_size = edge_overlap_size(e1, e2)
        if overlap_size and overlap_size < overlap_limit:
            return (e1, e2, overlap_size)
    return False


def edge_overlap_size(e1, e2):
    (start1, end1), (start2, end2) = [(eattr['offset'], eattr['end']) for fseq, eattr in (e1, e2)]
    if start2 < start1:
        (start1, end1), (start2, end2) = (start2, end2), (start1, end1)
    if end1 <= start2:
        return None
    return end1 - start2


def annotate_fragment_competitors(ptoc, ctop, rcgraph, strands, fcompg=None):
    """
    For the strands that a fragment binds to (i.e. ctop edges with rcompl=True)
    find other fragments that also binds to the same parts of those strands.
    """
    if fcompg is None:
        #fcompg = nx.MultiDiGraph()
        fcompg = nx.DiGraph()
        fcompg.graph['graph'] = {'rankdir': 'LR'}
    for fseq in leaf_node_gen(ptoc):
        fattr = ptoc.node[fseq]
        if 'base_competitors_count' not in fattr:
            fattr['base_competitors_count'] = [0]*len(fseq)
        if 'base_competitors' not in fattr:
            fattr['base_competitors'] = [set() for _ in range(len(fseq))]
        fbcc = fattr['base_competitors_count']
        fbcomp = fattr['base_competitors']
        for pseq in ctop[fseq]:
            eattr = ptoc[pseq][fseq]
            if eattr.get('rcompl') is True:
                # We only want to consider sense fragments:
                continue
            fpoffset, fpend = eattr['offset'], eattr['end']
            for fseq2, eattr2 in ptoc[pseq].items():
                if fseq2 == fseq or ptoc[pseq][fseq2].get('rcompl') is True:
                    # We only want to consider sense fragments:
                    continue
                poffset2, pend2 = eattr2['offset'], eattr2['end']
                overlap_size = edge_overlap_size((fseq, eattr), (fseq2, eattr2))
                if overlap_size: # it overlaps
                    # fragment index where the overlap starts. starts at base fseq[0].
                    # fpoffset and poffset2 idxs are all relative to parent seq.
                    # By subtracting fpoffset we get idx relative to fragment fseq.
                    overlap_fstart = max((0, poffset2-fpoffset)) # 0 if poffset2 <= fpoffset else poffset2 - fpoffset
                    # offset and end can be negative for rcompl edges:
                    assert pend2 > 0
                    #overlap_fend = overlap_size if pend2 >= fpend else overlap_size - (fpend - pend2)
                    overlap_fend = min((pend2-fpoffset, fpend-fpoffset))
                    for bidx in range(overlap_fstart, overlap_fend):
                        fbcc[bidx] += 1
                        fbcomp[bidx].add(fseq2)
                        #fcompg.add_edge(fseq, fseq2, bidx=bidx)
                    fcompg.add_edge(fseq, fseq2, overlap_fstart=overlap_fstart, overlap_fend=overlap_fend)
    return fcompg


def clear_fragment_competition_annotations(ptoc):
    for node, fattr in ptoc.nodes(data=True):
        if 'base_competitors_count' in fattr:
            del fattr['base_competitors_count']
        if 'base_competitors' in fattr:
            del fattr['base_competitors']



def split_overlapping_fragments(ptoc, ctop, rcgraph, strands):
    needs_splitting_up = defaultdict(set)
    # aborted_splits = set()
    for strand in strands:
        # Make base index => fragments map
        base_to_fragments = defaultdict(set)
        pseq = strand['seq']
        sense_fragments = [(fseq, eattrs) for fseq, eattrs in ptoc[pseq].items() if eattrs['rcompl'] is False]
        for fseq, eattrs in ptoc[pseq].items(): # graph[node] => dict with {target: edge_attrs}
            # graph.node[node] => node_attrs
            # Probably don't have to include fragments that are purely overlap and has no rcompl?
            # In fact, only break down sense fragments.
            if eattrs['rcompl']:
                continue
            offset, end, length = eattrs['offset'], eattrs['end'], ptoc.node[fseq]['length']
            assert len(fseq) == length
            for i in range(offset, offset+length):
                base_to_fragments[i].add((fseq, pseq, offset, end, True))
        overlap_length = sum(len(fset)-1 for fset in base_to_fragments.values())
        # dont_split_fseq = None
        # if overlap_length == 0: # No overlaps
        #     print("No overlapping fragments on strand %s" % (strand['Name'], ))
        #     continue
        # elif overlap_length <= promptlimit:
        #     print("Two or more fragments are overlapping on strand %s:\n  %s" % (strand['Name'], strand['seq']))
        #     for i, (fseq, eattrs) in enumerate(sense_fragments, 1):
        #         print(i, " "*eattrs['offset'], fseq)
        #     print("However, the overlap is less than the specified promptlimit %s." % promptlimit)
        #     answer = ""
        #     while not answer:
        #         answer = input("Do you still want to split this strand up? "
        #                        "(y/n or type the number for the strand you want to keep intact) ")
        #     if answer.lower()[0] == 'n':
        #         continue
        #     elif answer.lower()[0] == 'y':
        #         pass
        #     else:
        #         dont_split_fseq = sense_fragments[int(answer)-1][0]
        # else:
        #     print("Overlap %s > promptlimit %s" % (overlap_length, promptlimit))
        # Go over every index of pseq and determine if that base has multiple fragments.
        for idx, fset in base_to_fragments.items():
            # fset is a set of tuples, one tuple for each fragment present at this base idx.
            if len(fset) > 1:
                # first check whether any tups are unnecessary:
                # Actually, might as well do this in a separate function before doing this...
                # obsolete = set()
                # for tup in fset:
                #     other_tups = [otup for otup in fset if otup != fset]
                #     for otup in other_tups:
                #         # offset and end:
                #         if otup[2] < tup[2] and otup[3] > tup[3]:
                #             obsolete.add(tup)
                # if len(obsolete) > 0:
                #     fset -= obsolete
                for tup in fset: # tup = (fseq, pseq, offset, end)
                    # if dont_split_fseq:
                    #     if tup[0] == dont_split_fseq:
                    #         print("(Not splitting fragment %s)" % dont_split_fseq)
                    #         continue
                    #     # Modify fset to set breaking=False
                    #     fset = {tup[0:4] + (False,) for tup in fset}
                    needs_splitting_up[tup] |= fset
    # How to know what parts are covered by the "non-broken" fragment??
    # We now have a map:
    # fragment => {set of other fragments with shared sub-fragments}
    subfragments = set()
    for (fseq, pseq, fstart, fend, breaking), other_frags in needs_splitting_up.items():
        break_idx = set((fstart, fend))
        for other_frag in other_frags:
            # if other_frag[4] == False:
            #     # this overlapping strand does not actually break up the
            #     continue
            # other_frag = (fseq, pseq, fstart, fend) tuple for other fragment
            for idx in other_frag[2:4]: #start, end = other_frag[2:3]
                if fstart < idx < fend:
                    break_idx.add(idx)
        break_idx = sorted(break_idx)
        # aborted_breaks = {other_fset for other_fset in other_frags
        #                   if other_fset[0] != fseq and }
        if len(break_idx) <= 2:
            # Can happen if we have one domain that is fully covered by another domain.
            # This can happen after re-partitioning fragments.
            # There might be some bugs here; I'll weed them out as/when/if-needed.
            # Edit: Doing this in a separate "cleaning" function...
            # for oseq, opseq, ofstart, ofend, obreaking in other_frags:
            #     if ofstart < fstart and ofend > fend:
            #         # other_fragment fully covers frag on pseq:
            #         assert pseq == opseq
            #         ptoc.remove_edge(pseq, fseq)
            #         ctop.remove_edge(fseq, pseq)
            pass
        else:
            for new_start, new_end in zip(break_idx, break_idx[1:]):
                sf_seq = fseq[new_start-fstart:new_end-fstart]
                subfragments.add(sf_seq)
                assert sf_seq == pseq[new_start:new_end]
                print("Creating sub-fragment: %s from fragment %s (pseq=%s)" % (sf_seq, fseq, pseq[0:6]+"..."))
                ptoc.add_node(sf_seq)
                ctop.add_node(sf_seq)
                ptoc.add_edge(fseq, sf_seq, rcompl=False, undo=not breaking,
                              offset=new_start, end=new_end, foffset=new_start-fstart, fend=new_end-fstart)
                ctop.add_edge(sf_seq, fseq, rcompl=False, undo=not breaking,
                              offset=new_start, end=new_end, foffset=new_start-fstart, fend=new_end-fstart)



def extract_domains(ptoc, strands, rcgraph=None, namegen=None):
    """
    """
    if namegen is None:
        nleafs = len(leaf_nodes(ptoc))
        uppers = string.ascii_uppercase
        # combinations = 26^r >= nleafs  <==> r = ln(nleafs)/ln(26) = log26(nleafs)
        r = math.log(nleafs, len(uppers))
        import itertools
        if r > 1:
            namegen = ("".join(tup) for tup in itertools.combinations(uppers, r))
        else:
            namegen = iter(uppers)

    for strand in strands:
        root = strand['seq']
        strand['domain_seqs'] = get_strand_domains(ptoc, root)
        if namegen:
            name_domains(ptoc, rcgraph, domainnodes=strand['domain_seqs'], namegen=namegen)
        strand['domain_names'] = [ptoc.node[dseq]['Name'] for dseq in strand['domain_seqs']]


def gen_strand_domains(ptoc, root, offset=0):
    """
    Generate domain names.
    Assumes that nodes have been added for single/unpaired sequences.
    Arguments:
    :ptoc: parent-to-child graph - must be acyclic (tree/trie) graph.
    :root: node to start from in the tree.
    """
    if len(ptoc[root]) == 0:
        yield (offset, root)
    else:
        # Not sorted - do this afterwards...
        nodes = sorted([(eattr['offset'], fseq) for fseq, eattr in ptoc[root].items() if eattr['rcompl'] is False])
        #return [node for offset, fseq in nodes for node in get_strand_domains(ptoc, fseq)]
        for offset, fseq in nodes:
            # or python 3.3+: yield from gen_strand_domains(ptoc, fseq, offset)
            for offset, leaf in gen_strand_domains(ptoc, fseq, offset):
                yield (offset, leaf)



def get_strand_domains(ptoc, root):
    """
    Assumes that nodes have been added for single/unpaired sequences.
    ptoc must be acyclic (tree/trie) graph.
    """
    # if len(ptoc[root]) == 0:
    #     return [(None, root)]
    # nodes = sorted([(eattr['offset'], fseq) for fseq, eattr in ptoc[root].items() if eattr['rcompl'] is False])
    # return [(offset, leaf) for offset, fseq in nodes for offset2, leaf in get_strand_domains(ptoc, fseq)]
    return [fseq for offset, fseq in sorted(set(gen_strand_domains(ptoc, root)))]


def name_domains(ptoc, rcgraph, domainnodes=None, namegen=None,
                 rcnamefunc=lambda name: name.lower() if name.upper() == name else name.upper()):
    """
    :ptoc: graph going from parent (larger) strands/fragments to child fragment/domain.
    :rcgraph: graph where complementary domains/fragments are 2-way connected.
    """
    if domainnodes is None:
        domainnodes = [node for node in ptoc.nodes() if len(ptoc[node] == 0)]
    for node in domainnodes:
        if ptoc.node[node].get('Name') is None:
            name = next(namegen)
            print("Name: %s -> %s" % (name, node))
            ptoc.node[node]['Name'] = name
            if node in rcgraph:
                rcname = rcnamefunc(name)
                for rseq, reattr in rcgraph[node].items():
                    print(" - name rcompl: %s -> %s" % (rseq, rcname))
                    ptoc.node[rseq]['Name'] = rcname
        else:
            print("Domain node %s already has name: %s" % (node, ptoc.node[node]['Name']))
            # name_domain(ptoc, rcgraph, node, name) # Takes care to also name complement nodes


def name_domain(ptoc, rcgraph, domain, name, rcname=None,
                rcnamefunc=lambda name: name.lower() if name.upper() == name else name.upper()):
    """
    Name a single domain and propagate that change to reversed complement domains in rcgraph.
    """
    ptoc.node[domain]['Name'] = name
    if domain in rcgraph.node:
        unnamed = [(node, eattr['rcompl']) for node, eattr in rcgraph[domain].items()
                   if ptoc.node[domain].get('Name') is None]
        for node, rcompl in unnamed:
            kwargs = {'name': (rcname if rcname is not None else rcnamefunc(name)) if rcompl else name,
                      'rcname': name if rcompl else rcname,
                      'rcnamefunc': rcnamefunc}
            name_domain(ptoc, rcgraph, node, **kwargs)



def invert_regions(regions, base=None):
    """ This would have been easier using a bitarray... """
    if base is None:
        base = [0, max(region[1] for region in regions)]
    else:
        base = list(base)
    regions.sort()
    regions = deque(regions)
    inverted = [base] # regions not covered by :regions:
    # for region in regions:
    #     code
    while regions:
        first = regions.popleft()
        if first[0] <= inverted[-1][0]:
            inverted[-1][0] = first[1]
        elif first[0] <= inverted[-1][1]:
            new_inv = [first[1], inverted[-1][1]]
            inverted[-1][1] = first[0]
            inverted.append(new_inv)
    if inverted[-1][0] == inverted[-1][1]:
        # Remove empty last region:
        inverted.pop()
    return inverted



def unify_all_regions(regions):
    regions.sort()
    regions = deque(regions)
    unified = []
    first = regions.popleft()
    while regions:
        second = regions.popleft()
        first, second = unify_regions(first[0], first[1], second[0], second[1])
        if second is not None:
            unified.append(first)
            first = second
    unified.append(first)
    return unified


def unify_regions(start1, end1, start2, end2):
    assert start1 < end1 and start2 < end2
    if start2 < start1:
        if end2 < end1:
            start1, end1, start2, end2 = start2, end2, start1, end1
        else:
            # start2 < start1 < end1 < end2:
            return (start2, end2), None
    if start1 <= start2 and end2 <= end1:
        return (start1, end1), None
    elif start1 < start2 < end1 <= end2:
        return (start1, end2), None
    return (start1, end1), (start2, end2)



def make_graphviz_graph(filename, G, *more, prog="dot", outputtype="png", open_png=True,
                        add_same_rank=False, add_labels=False, add_style=True, add_rank_dir=None):
    """
    Make a graphviz dot file and process it with :prog: to produce render the graph as an image.
    :add_labels: can be True or a non-empty dict. If it is a non-empty dict,
                 **add_labels is passed to make_graphviz_labels.
    """
    if add_rank_dir:
        # should be one of "TB", "LR", "BT", "RL", c.f. http://www.graphviz.org/doc/info/attrs.html
        G.graph['graph'] = {'rankdir': add_rank_dir}
    if more:
        print("Merging %s with %s before parsing with graphviz..." % (G, more))
        G = G.copy()
        for g in more:
            if G.is_multigraph() and not g.is_multigraph():
                g = nx.MultiDiGraph(g) if G.is_directed() else nx.MultiGraph(g)
            G.add_nodes_from(g.nodes(data=True))
            G.add_edges_from(g.edges(data=True))
    if add_labels:
        if add_labels is True:
            add_labels = {}
        make_graphviz_labels(G, **add_labels)
    if add_style:
        apply_style(G, add_style)
    err = None
    if add_same_rank:
        rank_groups = defaultdict(list)
        for node, attrs in G.nodes(data=True):
            if 'rank' in attrs:
                rank_groups[attrs['rank']].append(node)
        if add_same_rank == 'pydot' or G.is_multigraph():
            # Use the pydot interface to graphviz:
            P = nx.to_pydot(G)
            # from networkx.drawing.nx_pydot import load_pydot
            # pydot = load_pydot()
            # DOESN'T WORK
            # for rank, seqs in rank_groups.items():
            #     S = pydot.Subgraph(rank='same')
            #     for seq in seqs:
            #         nodes = P.get_node(seq) # returns a list of nodes with the given name.
            #         if len(nodes) != 1:
            #             pdb.set_trace()
            #         node = nodes[0]
            #         # We get multiple nodes with the same pydot 'sequence' identifier
            #         #node = pydot.Node(seq) # This makes dot stop working during write_dot...
            #         S.add_node(node)
            #     P.add_subgraph(S)
            P.write_dot(filename) # same as P.write(filename, format='dot', prog=None). Default prog is 'dot'.
            is_ok = P.write_png(filename + "." + outputtype, prog='dot') # returns True on success.
        else:
            # Gohlke's pygraphviz for Windows doesn't seem to support strict=False.
            # Use the pygraphviz Agraph interface to graphviz:
            A = nx.to_agraph(G)
            for nodes in rank_groups.values():
                if len(nodes) > 1:
                    A.add_subgraph(nodes, rank="same")
            A.write(filename)
            # Use pygraphviz to invoke graphviz dot (simply does a call with subprocess module...)
            # returns None if written to file; returns data if path=None; raises exception if error.
            try:
                A.draw(path=filename + "." + outputtype, format='png', prog='dot')
            except OSError as e:
                err, is_ok = e, False
            else:
                is_ok = True
    else:
        nx.write_dot(G, filename) # uses pydot by default
        cmd = [prog, "-T%s" % outputtype, "-O", filename]
        err = subprocess.call(cmd) # use check_output if you need the output
        is_ok = not err
    if is_ok:
        if open_png:
            webbrowser.open(filename + "." + outputtype)
    else:
        print("Graphviz graph drawing not successful:", is_ok, err)


def apply_style(G, style=None):
    """
    Apply style to graph :G:
    :G: is a networkx graph.
    :style: is a dict with keys 'graph', 'nodes', and 'edges'.
            G.graph is updated with elements from style['graph'] (if given),
            all node attrs is updated with style['nodes'] (if given),
            all edge attrs is updated with style['edges'] (if given).
    """
    if style is None or style is True:
        style = {'nodes': {'shape': 'box'}}
    if 'graph' in style:
        G.graph.update(style['graph'])
    if 'nodes' in style:
        for n, attr in G.nodes(data=True):
            attr.update(style['nodes'])
    if 'edges' in style:
        for n, attr in G.edges(data=True):
            attr.update(style['edges'])




def make_graphviz_labels(G, edge_labels=True, xlabel_key='Name'):
    """ Add labels (node and edge label/xlabel attrs) to graph. """
    if edge_labels is True:
        for s, t, eattr in G.edges(data=True):
            eattr['label'] = "".join(
                str(v) for v in
                [eattr['offset'] if 'offset' in eattr else "",
                 (", %s" % eattr['end']) if 'end' in eattr else "",
                 " (rev)" if 'rcompl' in eattr and eattr['rcompl'] else "",
                 (" at %s" % eattr['bidx']) if 'bidx' in eattr else "",
                 (" %s-%s" % (eattr['overlap_fstart'], eattr['overlap_fend'])) if 'overlap_fstart' in eattr else ""])
    if xlabel_key:
        for n in G.nodes():
            if xlabel_key in G.node[n]:
                G.node[n]['xlabel'] = G.node[n][xlabel_key]


def leaf_nodes(G):
    return [n for n in G.nodes() if G.out_degree(n) == 0]


WC = dict(zip("ATGC", "TACG"))
def dna_rcompl(seq):
    """ Reversed complement of DNA sequence. """
    return "".join(WC[b] for b in reversed(seq))


def find_compl(G, rcg=None):
    if rcg is None:
        rcg = G
    for n1 in G.nodes():
        for n2 in G.nodes():
            if n2 == "".join(WC[b] for b in reversed(n1)):
                rcg.add_edge(n1, n2, rcompl=True)


def find_overlaps(seq1, seq2):
    """
    https://docs.python.org/3/library/difflib.html#sequencematcher-objects
    https://pypi.python.org/pypi/pydna/0.9.9
    https://pypi.python.org/pypi/biopython
    http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/
    http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/diffseq.html
    http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/seqmatchall.html
    http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/wordcount.html

    """
    sm = SequenceMatcher(a=seq1, b=seq2)
    blocks = sm.get_matching_blocks()
    return blocks


def export_strands(strands, path, ptoc=None, fields=('Name', 'domain_names', 'domain_seqs'),
                   header=None, sep="\t", sep2=", ",
                   list_fields={'domain_names': 'Domains', 'domain_seqs': 'Sequence'}):
    """
    :strands:     List of dicts, one for each strand, with fields describing what to be exported.
    :path:        Write to this file.
    :ptoc:        If field names are not found in the strand dict, fall back to this graph for determining field values.
    :fields:      List of fields to export.
    :header:      If given, use this as header for each field column instead of the names given in fields/list_fields.
    :sep:, :sep2: Separators used as sep.join(fields) or sep2.join(list) for each strand in strands.
    :list_fields: dict with key: value pairs.
                  Keys are fields where strand[field] is a list to be converted to a string via sep2.join.
                  Values are the used as header name in the file (Can be None to just use the original field name).
    """
    if header is None:
        header = [list_fields.get(field) or field for field in fields]
    with open(path, 'w') as fp:
        fp.write(sep.join(header) + "\n")
        for strand in strands:
            for key in list_fields:
                if key not in strand:
                    strand[key] = ptoc.node[strand['seq']][key]
            fp.write(sep.join(sep2.join(strand[field]) if field in list_fields else strand[field]
                              for field in fields) + "\n")



def main():
    import re
    modre = re.compile(r"(?P<mod5p>\/\w+\/)?(?P<seq>[ATGC]*)(?P<mod3p>\/\w+\/)?")
    """
    """
    filename = "oligos.txt"
    basename, ext = os.path.splitext(filename)
    lines = """
Name	Sequence
AJB023	/5DiGN/TTTTTTTTTTCGTGTAGCCAATTAGACTGA
AJB024	/5Phos/TCTAACTTACAGAGCATGGCTTTTTTTTTT/3Bio/
BN045	/5Phos/TG CAC AAC TGT AAG GTC CTT CCG CCG GGC GTG AGT AGG GTT GAC TAA GAG TTT TCT CTT AGT CAA CCC TAC TCA CGC CCG GCG GTT GTC ACG CAG ACG ACG GCC CAG GGA GAC CAC TGA C
BN046	/5Phos/GA GTA TCA AAC CCG TAC ATT AAT TAA TAT AAG CAT TTT CGG AGG TTC TCT TTT TAG AGA ACC TCC GAA AAT GCT TAT ATT AAT TTT AAG GCA ACA TCC TCG GCA TCA GGG TTA CTT TTT G
BN047	AAA AGC CAT GCT CTG TAA GTT AGA CAA AAA GTA ACC CTG ATG CC
BN048	GAG GAT GTT GCC TTT TTG TAC GGG TTT GAT ACT CGT CAG TGG TCT CCC TGG GCC
BN049	GTC GTC TGC GTG ACT TGG ACC TTA CAG TTG TGC ACT CAG TCT AAT TGG CTA CAC GAA AA
""".strip().split("\n")
    rows = (line.strip().split("\t") for line in lines)
    header = next(rows)
    strands = [dict(zip(header, row)) for row in rows]
    for sd in strands:
        match = modre.match(sd['Sequence'].replace(" ", ""))
        assert match is not None
        sd.update(match.groupdict())
        sd.pop('Sequence')

    ptoc, ctop, rcgraph = domains_from_sequences(strands)

    path = basename + "_stdoms" + ext
    export_strands(strands, path, ptoc=ptoc)





if __name__ == '__main__':
    main()
