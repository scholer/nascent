# -*- coding: utf-8 -*-
##    Copyright 2016 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
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

# pylint: disable=C0103,W0142,W0212,R0902

"""

OLD looptracking module, for regular InterfaceGraph, not InterfaceMultiGraph.


Provides classes and functions for dealing with Complex loops to calculate loop energies.

Contains LoopGraph (former InterfaceGraph), IfNodes classes.


Note from commit ec5fc4b:
    Next up:
    - Make InterfaceGraph a MultiGraph where domain edges are un-keyed (causing overlap when merging) while backbone
      edges between domains are keyed by (hXend3p, hXend5p). That should help prevent loop path degeneracies.
    - Move loop tracking to InterfaceGraph. I'm also considering composing InterfaceGraph on a per-complex basis
      rather than on a system level. Usually when i'm using the interface graph, it is in the context of a Complex.
      (I'm currently using in Graphmanager to shortest-paths for new loops, and I'm modifying interface_graph in
      ComponentMgr.(de)hybridize/(un)stack, but that could arguably be better done in the complex-specific
      join/break_complex_at methods anyway.) This would allow me to move all the scattered Ifnode-delegation and
      loop tracking to a single, self-contained class. Maybe call this class "MultiLoopGraph" or "InterfaceLoopGraph"
      or just "LoopTracker" (composing the interface_graph as an attribute)?
    - Note: I'm really not that concerned with the "cost of Complex object instantiation", since I expect to move to
      a "1-permanent-complex-per-stand with delegation" scheme. But this is currently a pre-mature optimization.



## InterfaceGraph: System-level vs complex-level, pros and cons: ##

Complex-level pros:
- We often use InterfaceGraphs very explicitly in the context of complexes, e.g. in Cmplx.effectutate_loop_changes()

System-level pros:
- We can add strands to the system-level graph just once and then we no longer have to add/remove any nodes,
    instead we are just adding/removing edges.
- In fact, since both hybridization and stacking interactions are represented by ifnodes merging,
    we don't even have to create or delete any edges (in terms of attr dict objects),
    we just re-arrange the edge's source or target node.

Maybe as a "compromise": Don't use classes/methods at all for loop tracking, but have
*everything* loop-related exposed as functions operating operating on arguments (complex, interface_graph).

Or, if we need configurational parameters, have a LoopTracker class. But, no separate LoopGraph class.


"""



from __future__ import absolute_import, print_function, division
from collections import defaultdict, Counter, deque
import itertools
from itertools import groupby
from operator import itemgetter
import networkx as nx
# from networkx.algorithms.shortest_paths import shortest_path
from networkx import single_source_dijkstra, dijkstra_path
# This is actually the shortest path function that we are always using
# shortest_path = dijkstra_path # defaults to 'weight' as edge length
import pdb
import math
from math import log as ln
from pprint import pprint

# Local imports:
from .domain import DomainEnd
from .system_graphs import InterfaceMultiGraph
from .domain import DomainEnd
from .utils import (sequential_number_generator, tupleify)
from .utils import prod
from .utils import unique_gen
# from .constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION)
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        AVOGADRO_VOLUME_NM3,
                        #HELIX_XOVER_DIST, HELIX_WIDTH, HELIX_STACKING_DIST,
                        chain_model_gamma_exponent as gamma)


## Module constants:
SEGMENT_STIFFNESS, SEGMENT_LENGTHS, SEGMENT_CONTOUR_LENGTH, SEGMENT_LENGTH_SQUARED = 0, 1, 0, 1

## Module singletons:
loopid_sequential_number_generator = sequential_number_generator()


def shortest_path(interface_graph, source, target, weight='len_contour'):
    return dijkstra_path(interface_graph, source, target, weight)
    # dijkstra_path calculates length to all nodes until target is found.
    # return single_source_dijkstra(interface_graph, source, target, 'len_contour')
    # single_source_dijkstra returns a length and paths dict with lengths from source to [node]
    # perhaps check single_source_dijkstra and _dijkstra to see if anyting can be optimized


def shortest_multigraph_path_edges(interface_graph, source, target, weight='len_contour'):
    """
    Return a tuple of (path, edges) where path is the shortest path between source and target,
    as a list of N nodes: [source, node1, node2, node3, ..., target]
    and edges is a list of N-1 edges connecting each pair of nodes in path, such that:
        [list of (src, tgt, key) = zip(path, path[1:], edges)]
    """
    path = dijkstra_path(interface_graph, source, target, weight)
    edges = [min([(eattr[weight], key) for key, eattr in interface_graph[n1][n2]])
             for n1, n2 in zip(path, path[1:])]
    return path, edges



# TODO: Remove reacted_ifnodes - only used for debug assertions.
def effectuate_loop_changes(cmplx, loop_effects, is_forming, reacted_ifnodes=None, is_copy=False):
    """
    loop_effects is a dict describing what the effects of
    a new loop being formed or an existing loop being broken.
    it should contain a 'changed_loops_by_hash' entry which is a dict with:
    loop_hash:

    Note: The loop hashes applies to the state *before* the loop is formed, so in order to re-create
    loops, this method should be called BEFORE asserting state change (but AFTER the reaction has
    actually been performed, i.e. the interaction has been broken).

    The timing of this needs to be exactly right. Make sure you keep track of when and where
    hashes (state specs) are re-calculated for ifnodes and loops:
    ifnode state fingerprint is reset and re-calculated:
        IfNode.state_fingerprint() is called by:
            Complex.calculate_loop_hash(loop_path) similarly calls ifnode.state_fingerprint() for all src,tgt.
            Complex.rebuild_loopid_by_hash_index() calls ifnode.state_fingerprint(), which simply calls
        ifnode.state_fingerprint() simply calls:
            ifnode.domain_end.state_fingerprint() for ifnode in IfNode.delegated_edges.keys()
        DomainEnd.state_fingerprint() calls:
            (DomainEnd.domain.state_fingerprint(), DomainEnd.end, DomainEnd.stack_partner is not None)
        Domain.state_fingerprint() uses:
            Domain._specie_state_fingerprint - which is recalculated if None as:
            Domain._specie_state_fingerprint = (dspecie, self.partner is not None, c_state, in_complex_identifier)
        Domain._specie_state_fingerprint reset by (and only by):
            Domain.state_change_reset()
        Domain.state_change_reset is currently called from multiple places:
            Complex.reset_state_fingerprint()
            ReactionMgr.hybridize_and_process (but not stack_and_process?) - debugging...
            ComponentMgr.join/break_complex_at() - if Domain.strand.complex is None
            ReactionMgr.post_reaction_processing - but only for free strands (!)
        Complex.reset_state_fingerprint is called from:
            ComponentMgr.join/break_complex_at()
            Complex.assert_state_change()
            # Complex.get_ifnode_by_hash -- but just for debugging
            # ReactionMgr.hybridize_and_process (but not stack_and_process?) - but only temporarily, for debugging...

    Now, we need to be able to update Complex loops based on the abstract hash-based
    loop_hash and path_spec values loop_effects.
    These hashes are all based on the Complex state *before* performing the reaction.
    Thus, we need to either effectuate_loop_changes(loop_effects) or at least update loop_effects
    before calling Complex.reset_state_fingerprint(), i.e.:
        in ComponentMgr.join/break_complex_at, or
        in ReactionMgr.post_reaction_processing() before calling Complex.assert_state_change.

    I think I will try to effectuate_loop_changes() it in ComponentMgr.join/break_complex_at
        Edit: Is a little tedious since we do not have reaction_spec_pair
        ComponentMgr generally doesn't know much about reactions...
        We DO calculate dHdS for stacking and hybridization in ComponentMgr.hybridize/stack,
        but they are state independent and does not rely on reaction_spec_pair.

    Moved effectuate_loop_changes invocation to ReactionMgr.post_reaction_processing

    Note: Tools to generate call graphs:
        pycallgraph:
        pyan: https://github.com/dafyddcrosby/pyan, https://github.com/davidfraser/pyan
        yapppi: https://bitbucket.org/sumerc/yappi/
    Profiling blog posts:
        TalkPython podcast #28, http://blog.thehumangeo.com/2015/07/28/profiling-in-python/
        PyCharm professional has a very nice profiler with call-graph:
        http://blog.jetbrains.com/pycharm/2015/05/pycharm-4-5-eap-build-141-988-introducing-python-profiler/

    """
    if not isinstance(loop_effects, dict):
        loop_effects = dict(loop_effects)
    if not is_copy:
        # loop_effects = deepcopy(loop_effects) # Don't, we have ifnodes in here which have recursive things.
        loop_effects = dict(tupleify(loop_effects))

    # pdb.set_trace()
    # TODO: WAAIT a moment, are you modifying the loop_effects directive??
    # TODO: Separate mutable vs immutable (and cachable) parts of the loop info!! Only add immutables to cache!

    if is_forming:
        # For forming reactions, loop_effects are calculated on the *previous* state.
        assert loop_effects['complex_state_fingerprint'] == cmplx._historic_fingerprints[-1]
        assert loop_effects['complex_loop_ensemble_fingerprint'] == cmplx._historic_loop_ensemble_fingerprints[-1]
    else:
        # For breaking reactions, loop_effects are calculated on the *current* state.
        assert loop_effects['complex_state_fingerprint'] == cmplx.state_fingerprint()
        assert loop_effects['complex_loop_ensemble_fingerprint'] == cmplx.loop_ensemble_fingerprint

        # It seems when we are breaking we should update the index:
        cmplx.rebuild_ifnode_by_hash_index()

    if reacted_ifnodes is None and is_forming:
        reacted_ifnodes_given = reacted_ifnodes # for debug..
        reacted_ifnodes = [cmplx.ifnode_by_hash[ifnode_hash] for ifnode_hash in loop_effects['reactant_node_specs']]

    # Add/remove loop:
    if is_forming:
        # 1a. Add new loop:
        new_loop_values = dict(loop_effects['new_loop'])
        # new_loop_values is the "new loop template" from the cached loop_effects directive.
        # We create a new loop (instance) and add the cached result values to it
        new_loop = {
            'path': None,
            'loop_hash': None,
            'activity': None,
            'loop_entropy': None, # TODO: Rename to "loop_entropy",
            # TODO: Remove loop debug_info when done debugging
            'debug_info': {
                'update_history': []
            }
        }
        # new_loop['debug_info'] = dict(new_loop.get('debug_info', []))
        # new_loop.update(new_loop_values)
        # only update specific keys:
        new_loop['loop_hash'] = new_loop_hash = new_loop_values['loop_hash']
        new_loop['activity'] = new_loop_values['activity']
        new_loop['loop_entropy'] = new_loop_values['loop_entropy']
        new_loop['debug_info']['update_history'].append(new_loop_values)
        # if 'debug_info' in new_loop:
        new_loop_id = next(loopid_sequential_number_generator)
        cmplx.loops[new_loop_id] = new_loop
        cmplx.loopid_by_hash[new_loop_hash] = new_loop_id

        # Make sure ifnode_by_hash is up-to-date:
        # TODO: Check if this is needed or even a good idea. Shouldn't the index already be up-to-date?
        # Yes, it is definitely a bad idea to re-build ifnode_by_hash index.
        # Even if most ifnode state hashes haven't changed, the ifnode delegation may have changed,
        # and we use THAT when re-building the index.
        # cmplx.rebuild_ifnode_by_hash_index()
        # Edit: We need an up to date index when re-building after loop breakage (=ifnode splitting)
        # To know the hash for the recently-undelegated ifnode.
        # Then again, maybe for loop-breakage, we should first assert the complex state change and THEN
        # effectuate the loop_breakage_effects directive. In this also makes sense considering the order
        # of when the directives are calculated:
        #   loop_formation_effects directives are calculated BEFORE the any actual reaction is done. (predictively)
        #   loop_breakage_effects directives are calculated AFTER a reaction has been done.
        # TODO: Split this method into two separate methods to highlight difference in when the methods are used.


        # Update loop path to use this Complex's ifnodes (instances):
        # new_loop_template_path = new_loop['path']
        new_loop['path'] = cmplx.recreate_loop_path_from_spec(new_loop_values['path_spec'])
        new_loop['debug_info']['path_before_ifnode_processing'] = tuple(new_loop['path'])

        # Testing for changes in ifnode delegation:
        assert new_loop['path'][0].top_delegate() == new_loop['path'][-1].top_delegate()
        new_loop['path'][0] = new_loop['path'][0].top_delegate()
        new_loop['path'].pop() # Remove the last ifnode (which is also the first after reaction)
        # Make sure path consists of top_delegates only:
        assert [ifnode.top_delegate() for ifnode in new_loop['path']] == new_loop['path']
        # new_loop['path'] = tuple(new_loop['path'])
        new_loop['debug_info']['path_after_ifnode_processing'] = tuple(new_loop['path'])
        if max(Counter(ifnode.top_delegate() for ifnode in new_loop['path']).values()) > 1:
            print("\n\nWARNING: new_loop['path'] has duplicate top_delegate ifnodes:")
            print(new_loop['path'])
            print([ifnode.top_delegate() for ifnode in new_loop['path']])
            print(Counter(ifnode.top_delegate() for ifnode in new_loop['path']).most_common(3))
            pdb.set_trace()
            # Isn't it natural that delegation will have changed? The loop effects were calculated before
            # the reaction was done, but delegation has now been updated.
            # Either (a) Take changes in delegation into account when calculating the path,
            # or (b) update the path to account for delegation changes.
            # (a) is probably hard since we can't know in advance which node is selected as top_delegate.
            # We don't consider the edge between reactants when calculating activity, so we don't need to modify.
            # In fact, we could just say that the loop path is just always


        # Make sure the loop hash haven't changed:
        # new_loop_hash_from_path = cmplx.calculate_loop_hash(new_loop['path'])
        # new_loop_hash_from_path_spec = cmplx.calculate_loop_hash(None, new_loop['path_spec'])
        # assert new_loop['loop_hash'] == cmplx.calculate_loop_hash(new_loop['path'])
        # Edit: you cannot re-calculate the loop hash from the path ifnodes and expect to get the same hash
        # if you have performed a reaction because then the ifnode's delegation scheme will have changed
        # and that is used when calling src/tgt.state_fingerprint().
        # If you want to rely on a cached hash, you have to use cmplx.loops[loopid]['loop_hash']
        assert new_loop['loop_hash'] == cmplx.calculate_loop_hash(new_loop_values['path_spec'])
        for ifnode in new_loop['path']:
            cmplx.ifnode_loopids_index[ifnode].add(new_loop_id)
    else:
        # 1b: Remove existing loop
        del_loop_hash = loop_effects['del_loop_hash']
        del_loop_id = cmplx.loopid_by_hash[del_loop_hash]
        if del_loop_id == loop_effects['del_loop_id']:
            print("del_loop_id == loop_effects['del_loop_id']: %s == %s" %
                  (del_loop_id, loop_effects['del_loop_id']))
        del_loop = cmplx.loops[del_loop_id]
        print("%s.ifnode_loopids_index before deleting ifnodes in del_loop['path'] = %s: %s" % (
            cmplx, del_loop['path'], cmplx.ifnode_loopids_index))
        # There should be at least ONE loop to delete:
        assert len(set.union(*[cmplx.ifnode_loopids_index[ifnode] for ifnode in del_loop['path']])) > 0
        for ifnode in del_loop['path']:
            try:
                cmplx.ifnode_loopids_index[ifnode].remove(del_loop_id)
            except KeyError:
                if reacted_ifnodes:
                    assert ifnode in reacted_ifnodes
                # I'm not sure whether the ifnode delegation is invariant. I think it is random ATM, but
                # ifnode delegation should probably be made deterministic based on state fingerprints..
                # Although that might be hard because delegation can take different paths.
                assert ifnode.state_fingerprint() in (
                    loop_effects['previous_top_delegates_specs'] | loop_effects['newly_undelegated_ifnodes_specs'])
                # pdb.set_trace()
                # print("This was unexpected! - no, not really.")
        assert not any(del_loop_id in loopids_set for loopids_set in cmplx.ifnode_loopids_index.values())
        del cmplx.loops[del_loop_id]
        del cmplx.loopid_by_hash[del_loop_hash]

        # Update ifnodes hashes: - No. See above reason for why not.
        # cmplx.rebuild_ifnode_by_hash_index()

    ifnode_loopids_before = dict(cmplx.ifnode_loopids_index)
    # 2. Update changed existing loops:
    # TODO: (a) Is this needed when breaking loops? Wouldn't this already have been updated by loop_breakage_effects?
    # TODO: (b) Split this out to separate method?
    # TODO: (c) We really shouldn't be looking just at loop_effects; it would be more reliable to go over all loops
    # TODO      currently touching the reacted ifnddes.
    # When doing ifnodes-are-merged loop path updates, it is better to go over all loops in the union of
    # set.union(cmplx.loops_by_ifnode[ifnode] for ifnode in reactant_node_specs)
    # set.union(cmplx.loops_by_ifnode[.ifnode_by_hash[ifnode_hash]] for ifnode_hash in reactant_node_specs)

    for loop_hash, loop_updated_values in loop_effects['changed_loops_by_hash']: #.items():  # edit, using immutable tuples
        new_loop_values = dict(loop_updated_values)

        # Find the loop instance from hash
        loop_id = cmplx.loopid_by_hash[loop_hash]
        loop = cmplx.loop_by_loopid[loop_id]
        old_loop_nodes = set(loop['path'])  # Make sure this is BEFORE you re-set loop path to new

        # Update the loop instance with new values
        # loop.update(loop_info) # Edit: Explicitly updating entries:
        loop['loop_hash'] = new_loop_values['loop_hash']
        loop['activity'] = new_loop_values['activity']
        loop['loop_entropy'] = new_loop_values['loop_entropy']
        loop['debug_info']['update_history'].append(new_loop_values)

        ## Update path (instances)
        loop['debug_info']['path_before_overwriting_to_changed_path'] = tuple(loop['path'])
        # Uh, there is a problem re-creating path ifnodes after a breaking reaction:
        # The newly_undelegated_ifnode isn't present in the ifnode_by_hash index.
        # Maybe make sure both top-delegates *and* constituent (electorate/delegator) ifnodes?
        # Edit: the problem is more that the ifnode_by_hash index wasn't up to date:
        loop['path'] = cmplx.recreate_loop_path_from_spec(new_loop_values['path_spec'])

        # Update loopids_by_ifnode index:
        new_loop_nodes = set(loop['path'])
        # Update loopids_by_ifnode index:
        # Remove old ifnode entries from cmplx.ifnode_loopids_index and add new:
        print("\nUpdating cmplx.ifnode_loopids_index =", cmplx.ifnode_loopids_index)
        print("old_loop_nodes:", old_loop_nodes)
        print("new_loop_nodes:", new_loop_nodes)
        print("old_loop_nodes - new_loop_nodes:", old_loop_nodes - new_loop_nodes)
        print("new_loop_nodes - old_loop_nodes:", new_loop_nodes - old_loop_nodes)
        for ifnode in old_loop_nodes - new_loop_nodes:
            print("- cmplx.ifnode_loopids_index[%s].discard(%s)" % (ifnode, loop_id))
            cmplx.ifnode_loopids_index[ifnode].discard(loop_id)
        for ifnode in new_loop_nodes: # - old_loop_nodes: # Edit: Just marking loop for all ifnodes on path
            print("- cmplx.ifnode_loopids_index[%s].add(%s)" % (ifnode, loop_id))
            cmplx.ifnode_loopids_index[ifnode].add(loop_id)
        print("cmplx.ifnode_loopids_index after update: ", cmplx.ifnode_loopids_index, '\n')

        # Note: Above only touches ifnodes on this loop's path, not other

        # At this point we haven't updated the loop indexes, so new loop ifnodes will not be present yet.
        # New ifnodes will typically be present for breaking reactions (is_forming=False)
        # TODO: You should probably make a broader more general check that the ifnode_loopids_index index is up-to-date
        # if is_forming:
        # path_ifnodes_not_in_index = [ifnode for ifnode in loop['path'] if ifnode not in cmplx.ifnode_loopids_index]
        # ifnode_loopids_index is a defaultdict: ifnode => set(of all loops going through ifnode)
        # ifnodes_not_registered_for_path = [ifnode for ifnode in loop['path'] if loop_id not in cmplx.ifnode_loopids_index[ifnode]]
        # Error encountered when unstacking somewhere else, then performing a stacking loop formation
        try:
            assert all(ifnode in cmplx.ifnode_loopids_index and loop_id in cmplx.ifnode_loopids_index[ifnode]
                       for ifnode in loop['path'])
        except AssertionError:
            print("\n\nWARNING (effectuate_loop_changes):\nProblem with Complex.ifnode_loopids_index index; rebuilding and asserting again...\n\n")
            print("loop['path']:", loop['path'])
            print("old_loop_nodes:", old_loop_nodes)
            print("new_loop_nodes:", new_loop_nodes)
            print("old_loop_nodes - new_loop_nodes:", old_loop_nodes - new_loop_nodes)
            print("new_loop_nodes - old_loop_nodes:", new_loop_nodes - old_loop_nodes)
            print("ifnode_loopids_index:", cmplx.ifnode_loopids_index)
            print("ifnode_loopids_before:", ifnode_loopids_before)
            pdb.set_trace()
            cmplx.rebuild_ifnode_loopids_index()
            assert all(ifnode in cmplx.ifnode_loopids_index and loop_id in cmplx.ifnode_loopids_index[ifnode]
                       for ifnode in loop['path'])
    # end updating changed loops (loops with modified paths)

    ## Update loops with changed ifnode delegations:
    loopids_with_ifnode_delegation_changes = set.union(
        *[cmplx.ifnode_loopids_index[ifnode] for ifnode in reacted_ifnodes])


    if is_forming is False:
        # TODO: Isn't this currently covered by GraphManager.loop_breakage_effects??
        # TODO: - A good reason to split this functionality out to separate methods!!
        # # This is probably a bit harder. Except if the new loop path starts and ends with the node that
        # # is being split/expanded. Then it is just about the same:
        # # Note: This really should not be needed for loop breakage directives
        # # where the path_spec was created *after* breaking the connection,
        # # and so the path_spec already reflects the *current* state.
        # # But: There may be loops touching splitted ifnodes that are not changed.
        # updated_loop_path = list(unique_gen(ifnode.top_delegate() for ifnode in loop['path']))
        # # Assert that ifnode delegatio doesn't change path:
        # assert loop['path'] == updated_loop_path
        # loop['debug_info']['path_after_ifnode_processing'] = tuple(loop['path'])
        # print("Loop %s path changes after forming a new intra-complex connection:" % loop_id)
        # print("Old      loop path:", loop['debug_info']['path_before_update'])
        # print("Changed  loop path:", loop['debug_info']['path_after_update_before_ifnode_processing'])
        # print("Delegate loop path:", loop['debug_info']['path_after_ifnode_processing'])
        pass
    else:
        # Remove all loops from ifnodes who's representation has been delegated to another ifnode:
        for ifnode in reacted_ifnodes:
            if ifnode.delegatee is not None: # i.e. it has been delegated to another
                print("\n\nTransferring all loops (%s) associated with ifnode %s to new top_delegatee to %s" % (
                    cmplx.ifnode_loopids_index[ifnode], ifnode, ifnode.delegatee))
                cmplx.ifnode_loopids_index[ifnode.delegatee].update(cmplx.ifnode_loopids_index[ifnode]) # NB: update is in-place, union returns a new set...
                cmplx.ifnode_loopids_index[ifnode].clear()
        for loop_id in loopids_with_ifnode_delegation_changes:
            # NOTE: THIS DOES NOT replace looping over loop_effects['changed_loops_by_hash'] - which is updating loops
            # with changed/shortened(formation)/extended(breaking) paths.
            # However, the functionality does overlap ATM
            loop = cmplx.loop_by_loopid[loop_id]

            # if is_forming:

            loop['debug_info']['path_after_update_before_ifnode_processing'] = tuple(loop['path'])
            # As with the new loop, we should update the path to reflect current ifnode delegation:
            # I think it is guaranteed that the first and last ifnodes in the new_loop path above
            # (before removing the last ifnode of the path) will be merged into a single ifnode.
            # There might be more nodes joined, e.g. for hybridization.
            # if is_forming:
            # Should be the easy part, just check for merged ifnodes by raising to t
            # How is the updated path? e3a+e2+e3b or e2+e3b+e3a or ..?
            updated_loop_path = list(unique_gen((ifnode.top_delegate() for ifnode in loop['path'])))
            loop['path'] = updated_loop_path
            loop['debug_info']['path_after_ifnode_processing'] = tuple(loop['path'])
            print("Loop %s path changes after forming a new intra-complex connection:" % loop_id)
            print("Old      loop path:", loop['debug_info'].get('path_before_overwriting_to_changed_path'))
            print("Changed  loop path:", loop['debug_info']['path_after_update_before_ifnode_processing'])
            print("Delegate loop path:", loop['debug_info']['path_after_ifnode_processing'])
            # else:
            # cmplx.ifnode_loopids_index[ifnode].discard(loop_id)


    ## Check that ALL loops (a) uses top_delegate ifnodes, and (2) does not contain duplicate ifnodes
    ## - DEBUG - TODO: Remove excessive assertion checks (or move to separate method).
    ## TODO: It IS possible to have loops that have the same "path" but who are still affected because of ifnode changes.
    for loop in cmplx.loops.values():
        top_delegate_path = [ifnode.top_delegate() for ifnode in loop['path']]
        assert top_delegate_path == loop['path']
        assert max(Counter(top_delegate_path).values()) == 1

    # Update loop hashes (state specs):
    # Note: We need the old loop state specs (fingerprints, hashes) when updating existing loops above
    # cmplx.rebuild_loopid_by_hash_index(update_ifnodes=False)

    ifnode_loopids_after_manual_update = {k: v for k, v in cmplx.ifnode_loopids_index.items() if len(v) > 0}
    cmplx.rebuild_ifnode_loopids_index() # Does not have any empty sets
    ifnode_loopids_after_rebuild = dict(cmplx.ifnode_loopids_index)
    assert ifnode_loopids_after_manual_update == ifnode_loopids_after_rebuild



    # TODO: Do we need to update loops when ifnode delegation changes?
    # E.g. if stacking, then four ifnodes are represented as one ifnode.
    # All four ifnodes will have the same hash, because ifnode.state_fingerprint uses all delegates.
    # So if a path includes the edge {src, tgt} and the edge is collapsed by e.g. stacking to a single ifnode,
    # then that is just {ifnode}, which should be OK.
    # What if we have a path that includes ifnode2 and ifnode2 is expanded to {ifnode2, ifnode3}?
    # The system-level InterfaceGraph only includes edges between the top delegates.
    # InterfaceGraph.delegate() will move the delegator's edges to the top delegatee.
    # So new paths DOES NOT include
    # But either way, the all affected paths should be updated to reflect the new ifnode top_delegate path, right?

# end effectuate_loop_changes()








def interfaces_shortest_path(interface_graph, ifnode1, ifnode2):
    """
    TODO: This should certainly be cached.
    """
    if isinstance(ifnode1, DomainEnd):
        ifnode1, ifnode2 = ifnode1.ifnode.top_delegate(), ifnode2.ifnode.top_delegate()
    return shortest_path(interface_graph, ifnode1, ifnode2)



def group_interfaces_path_by_stiffness(interface_graph, path):
    """
    Returns an iterator of structural elements based on a interface-level path (list of InterfaceGraph nodes),
    where neighboring edges with the same stiffness has been grouped:
    path_edges = [(stiffness, [(length, length_sq, stiffness, source, target), ...]),
                  (stiffness, [(length, length_sq, stiffness, source, target), ...]),
                  (stiffness, [(length, length_sq, stiffness, source, target), ...]),
                  ...]
    That is, if the path is:
        4 nm ss-backbone + 3 nm ss-backbone + 5 nm ds+backbone + 1 nm ss-backbone + 2 nm ds-backbone + 8 nm ds-backbone,
    The result is:
        [(0, [(4 nm, 16 nm2, source, target), (3 nm, 9 nm2, source, target)]),
         (1, [(5 nm, 25 nm2, source, target)]),
         (0, [(1 nm, 1 nm2, source, target)]),
         (1, [(2 nm, 4 nm2, source, target), (8 nm, 64 nm2, source, target)])]

    Where stiffness=0 indicates a list of single-stranded edges,
    stiffness=1 indicates a hybridized duplex edge, and stiffness=3 is used for helix-bundles.
    """
    ## First create a generator of (length, length_sq, stiffness, source, target)
    ## Then use itertools.groupby to group elements by stiffness
    ## https://docs.python.org/3/library/itertools.html#itertools.groupby
    path_source_target_eattr = ((source, target, interface_graph[source][target])
                                for source, target in zip(path, path[1:]))
    # We really should have all three lengths: len_contour, dist_ee_nm, and dist_ee_sq.
    path_tuples = ((edge_attrs['len_contour'], edge_attrs['dist_ee_sq'], edge_attrs['stiffness'], source, target)
                   for source, target, edge_attrs in path_source_target_eattr)
    stiffness_key = itemgetter(2) # Because stiffness = path_tup[2] for path_tup i path_tuples
    return itertools.groupby(path_tuples, key=stiffness_key)


def interfaces_path_segments(interface_graph, path):
    """
    Arguments:
        :path:  A list of ifnodes describing the path.
    Return a list of (stiffness, [length, sum_length_squared]) tuples
    where each tuple represents a path segment with the same stiffness.
    This basically reduces the path to its constituent parts, which are usually all
    that is needed to calculate end-to-end probability distribution functions.
    """
    grouped_path_edges = group_interfaces_path_by_stiffness(interface_graph, path)
    # Trim away stiffness, source, target from (length, length_sq, stiffness, source, target) tuples:
    # zip(*[(1,1), (2, 4), (3, 9), (4, 16)]) = [(1, 2, 3, 4), (1, 4, 9, 16)]
    return [(stiffness, [sum(lengths) for lengths in zip(*[tup[:2] for tup in elem_tuples])])
            for stiffness, elem_tuples in grouped_path_edges]

def grouped_path_edges_to_segments(grouped_path_edges):
    """
    Reference function, should be inlined.
    """
    return [(stiffness, [sum(lengths) for lengths in zip(*[tup[:2] for tup in elem_tuples])])
            for stiffness, elem_tuples in grouped_path_edges]


def join_two_edge_groups(interface_graph, part_a, part_b, simulate_reaction, reaction_elems=(0, None)):
    """
    Join two edges groups.
        Assumes that the reaction is between the start of part_a and the end of part_b, i.e.:
            [part_a_last_segment, ..., part_a_first_segment] <::> [part_b_last_segment, ..., part_b_first_segment]

    Will join the first segment/group of part_a a with the last segment/group of part b:
    The join is done in-place by updating part_b. Make sure to make a copy if you want to retain the original!
    For stacking reactions, this is quite simple.
    For hybridization reactions, this requires a bit more thought..

    :part_a: and :part_b: are lists of segment edges groups with the form of:
        [(segment_stiffness, [segment_edge_attr_tuples]), ...],
    e.g.
        [(1, [(e1_length, e1_len_sq, e1_stiffness, e1_source, e1_target), (e2 attr tuple), ...])
         (0, [edge tuples for next segment with stiffness=0]),
         ...]
    """
    ## Adjust path segments as they would be *after* reaction.
    part_a_first_group = part_a[0]
    part_b_last_group = part_b[-1]
    if simulate_reaction is STACKING_INTERACTION:
        # Here we can assume that stacking *goes through* the duplexes:
        # TODO: Not sure this assertion is correct!
        # One case where it is not the case is the backbone edge connecting a stacked ifnode (self-loop)
        assert part_a[0][SEGMENT_STIFFNESS] > 0
        assert part_b[-1][SEGMENT_STIFFNESS] > 0
        # First join the first segment of a to the last segment of b: (using extend in-place)
        part_b[-1][1].extend(part_b[0][1])
        # Then append the rest of the segments from part a to part b:
        part_b.extend(part_a[1:])
    ## TODO: Determine this for hybridization reactions.
    elif simulate_reaction is HYBRIDIZATION_INTERACTION:
        assert part_a[0][SEGMENT_STIFFNESS] == 0
        assert part_b[-1][SEGMENT_STIFFNESS] == 0
        # Remove the domain edge from part_a and part_b:
        domain_edge_a = part_a[0][1].pop(0)
        domain_edge_b = part_b[-1][1].pop()
        ## TODO: MAKE SURE TO GET EVERYTHING RIGHT:
        ## part_a    d1end2    d1end1           . part_a    d1end2    d1end1           .
        ##  <------o-o---------o                .  <------o-o---------o                .
        ##            \         \       part_b  .            \         \       part_b  .
        ##             o---------o-------<-<-<--.             o---------o-------<-<-<--.
        ##             d2end2    d2end1         .             d2end2    d2end1         .
        # length, len_sq, stiffness, source, target
        d1end1_ifnode = domain_edge_a[3]
        d1end2_ifnode = domain_edge_a[4]
        d2end1_ifnode = domain_edge_b[3]
        d2end2_ifnode = domain_edge_b[4]
        d1 = d1end1_ifnode.domain_end.domain
        d2 = d2end1_ifnode.domain_end.domain
        assert d1end1_ifnode != d1end2_ifnode != d2end1_ifnode != d2end2_ifnode
        assert d1end2_ifnode.domain_end.domain == d1
        assert d2end2_ifnode.domain_end.domain == d2
        assert d1end1_ifnode.domain_end.end == d2end2_ifnode.domain_end.end
        assert d2end1_ifnode.domain_end.end == d1end2_ifnode.domain_end.end
        # Use the ifnode guaranteed to be connected in interface_graph:
        duplex_tup = (d1.ds_len_contour, d1.ds_dist_ee_sq, 1, d2end1_ifnode, d1end2_ifnode)
        # Add new 1-element segment group to part_b:
        part_b.append((1, [duplex_tup]))
        # Append modified part_a to part_b:
        part_b.extend(part_a[1:])
        # NB: the newly hybridized duplex is still connected by flexible phosphate backbone links,
        # since it has not had any oppertunity to stack, so the surrounding segments are still flexible,
        # and I don't have to worry about joining the stiff duplex edge with neighboring segments.
    else:
        raise NotImplementedError("Only STACKING_INTERACTION and HYBRIDIZATION_INTERACTION currently supported.")
    return part_b



def join_edges_groups(interface_graph, edge_groups, simulate_reaction=None, reaction_elems=(0, None)):
    """
    Returns
        closed_loop_groups:  A list of segments making up the loop after the loop has been closed by the given reaction.
                             Each segment is a list of edges with the same stiffness.
    Note that edges_groups are slightly different from path segments:
        segments are tuples of    (stiffness, [segment_contour_length, segment_dist_ee_sq])
        path groups are tuples of (stiffness, [(length, length_sq, stiffness, source, target), ...])
        i.e. each path group has a full list of path edges for each stiffness group,
        while each segment only has information about the segment's contour length and sum-of-squared-links.
    Raises:
        IndexError if the path is already only a single, stiff duplex.
    Arguments:
    :loop_path: List of ifnodes describing the full loop path.
        (InterfaceGraph is not a multigraph so adj[source][target] is unambiguous.)
    :simulate_reaction: Assume that the first and last path elements have undergone this type of reaction.
        Only valid argument is None or STACKING_INTERACTION.
    :reaction_elems: If the elements undergoing the reaction described by simulate_reaction,
        add them here. (NOT CURRENTLY SUPPORTED)
    Alternative method names:
        join_path_groups, join_path_edge(s)_groups
        join_edge_segment_groups, join_segment_edge_groups
        close_loop_edges_groups
    """
    # edge_groups = group_interfaces_path_by_stiffness(path)
    # Grouped path edges is a list of:
    # path_edges_groups = [(stiffness, [(length, length_sq, stiffness, source, target), ...]), ...]

    ## Adjust path segments as they would be *after* reaction.
    if simulate_reaction is STACKING_INTERACTION:
        # A duplex end can either be connected through its own duplex (the domain backbone),
        # or throught the backbone connection to the next domain.
        # We have both cases for both duplex ends.
        first_group = edge_groups[0]
        last_group = edge_groups[-1]
        if last_group[SEGMENT_STIFFNESS] > 0 and first_group[SEGMENT_STIFFNESS] > 0:
            if len(edge_groups) == 1:
                # The path may or may not go through the duplexed double-helix for the ends being stacked:
                if first_group[SEGMENT_STIFFNESS] > 0:
                    # Trying to stack opposite ends of duplexed double-helix:
                    raise IndexError("Path consists of a single, stiff element, cannot circularize")
            else: # if len(edge_groups) > 1:
                ## Both duplex ends are connected via their own (stiff) duplex;
                ## we should put their segments "together" as they would be when stacked:
                # Remove the last segment and apply it onto the first.
                # Even though the group is a tuple, the edge attrs in group_tup[1] is still a mutable list, i.e:
                #       first_group = (first_group_stiffness, [list of edge tuples]
                first_group[1] += last_group[1]   # Extend first group edges with edges from last group
                edge_groups.pop()  # Remove the last segment and append it to the first
                # Edit: You should not sum the square, but re-calculate the square based on actual length:
                # first_group[1][1] = first_group[1][0]**2
                # We re-calculate sum-of-squares for stiff groups/segments later, not here.
        # else:
        # Ends are connected by one or more flexible segments, either a single phosphate-backbone link,
        # or through one or more unhybridized domains.
    ## TODO: Consider hybridization reactions, c.f. join_two_edge_groups(...) method.
    else:
        raise NotImplementedError("This needs to be implemented!")
    return edge_groups


def closed_loop_segments(interface_graph, loop_path, simulate_reaction=None, reaction_elems=(0, None)):
    """
    Like closed_loop_path_groups, but returns segments rather than path groups.
    It does this a little different, though:
        First get segments (not edge groups), then:
        * if STACKING reaction and first and last segments are both stiff, merge first and last segment.
    This is the order used by the original calculate_loop_activity.
    """
    # path_segments is a list of (stiffness, [length, sum_length_squared]) tuples
    # where stiffness is (0 for ssDNA, 1 for dsDNA and 2 for helix-bundles).
    path_segments = interfaces_path_segments(interface_graph, loop_path)

    ## Adjust path segments as they would be *after* reaction.
    if simulate_reaction is STACKING_INTERACTION:
        # A duplex end can either be connected through its own duplex (the domain backbone),
        # or throught the backbone connection to the next domain.
        # We have both cases for both duplex ends.
        first_segment = path_segments[0]
        if len(path_segments) > 1:
            # Remove the last segment and apply it onto the first.
            last_segment = path_segments[-1]
            # SEGMENT_STIFFNESS, SEGMENT_LENGTHS, SEGMENT_CONTOUR_LENGTH, SEGMENT_LENGTH_SQUARED = 0, 1, 0, 1
            if last_segment[SEGMENT_STIFFNESS] > 0 and first_segment[SEGMENT_STIFFNESS] > 0:
                ## Both duplex ends are connected via their own (stiff) duplex;
                ## we should put their segments "together" as they would be when stacked:
                # stiffness = first_segment[0] if first_segment[0] >= last_segment[0] else last_segment[0]
                # Even though segment is a tuple, the segment lengths in segment[1] is still a mutable list:
                #path_segments[0] = (first_segment[0], []
                path_segments.pop()  # Remove the last segment and append it to the first
                first_segment[1][0] += last_segment[1][0]   # Add segment contour lengths
                # first_segment[1][1] += last_segment[1][1]   # Add segment squared lengths
                # Edit: You should not sum the square, but re-calculate the square based on actual length:
                first_segment[1][1] = first_segment[1][0]**2
            # else:
                # If the duplexes are not connected via their own stiff/duplexed domains, then downstream
                # calculation using LRE/SRE should produce correct result...
        else:
            # Only a single path segment:
            if first_segment[SEGMENT_STIFFNESS] > 0:
                # We have a single, fully-stacked segment; this cannot ever stack back upon it interface_graph.
                # return 0
                raise IndexError("Path consists of a single, stiff element, cannot circularize")
            else:
                # A single, flexible connection...
                # This must be just a single phosphate connection between upstream 3' and downstream 5' ends, right?
                # Nah, it's the same for duplex hybridization of duplex domains separated by T-loop.
                # Thus, we also need to check for loops, etc... Let's just let the general routine take its course.
                pass
                # mean_sq_ee_dist = first_segment[1][1]
                # effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2)
                # activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
                # return 2.0
    return path_segments



def calculate_loop_activity(interface_graph, loop_path, simulate_reaction=None, reaction_elems=(0, None), circularize_path=False):
    """
    Returns
        activity - the "probability" that the loop ends are within a critical radius (relative to a 1 M solution).
    Arguments:
    :loop_path: List of ifnodes describing the full loop path.
        (InterfaceGraph is not a multigraph so adj[source][target] is unambiguous.)
    NOTE: loop_path should describe the FULL path, i.e. it should start AND end on the same ifnode.
    (simulate_reaction can be used to further TRIM AWAY overlapping path, but it will not extend the path).
    :simulate_reaction: Assume that the first and last path elements have undergone this type of reaction.
        Only valid argument is None or STACKING_INTERACTION.
    :reaction_elems: If the elements undergoing the reaction described by simulate_reaction,
        add them here. (NOT CURRENTLY SUPPORTED)
    :circularize_path: If True, transform path so it starts and ends on the same element: path = path + path[0:1]
    """
    # cyclize vs circularize?
    # path_segments is a list of (stiffness, [length, sum_length_squared]) tuples
    # where stiffness is (0 for ssDNA, 1 for dsDNA and 2 for helix-bundles).
    print("\ncalculate_loop_activity called with %s, %s, %s" % (loop_path, simulate_reaction, reaction_elems))
    # pdb.set_trace()
    if circularize_path:
        loop_path = loop_path + loop_path[0:1]
    try:
        path_segments = closed_loop_segments(interface_graph, loop_path, simulate_reaction, reaction_elems)
        print(" - path_segments:", path_segments)
    except IndexError:
        return 0

    ## Find LRE and SRE parts
    ## TODO: Replace "element" by "segment", which better describes that a segment is composed of one or more
    ## edges with same stiffness.
    try:
        # LRE: "Longest Rigid Element"; "SRE": "Sum of Remaining Elements".
        _, LRE_len, LRE_len_sq, LRE_idx = max((stiffness, elem_length, elem_len_sq, i)
                                              for i, (stiffness, (elem_length, elem_len_sq))
                                              in enumerate(path_segments)
                                              if stiffness > 0)
        # No need to check if LRE_len <= HELIX_WIDTH - InterfaceGraph only includes backbone connections.
    except ValueError:
        # No stiff elements: (summarizing segment_len_sq is OK for non-stiff, single-stranded segments)
        LRE_len, LRE_len_sq = 0, 0
        SRE_len = sum(elem_length for (stiffness, (elem_length, elem_len_sq)) in path_segments if stiffness == 0)
        SRE_len_sq = sum(elem_len_sq for (stiffness, (elem_length, elem_len_sq)) in path_segments if stiffness == 0)
    else:
        # Exclude LRE when calculating SRE length:
        # Don't sum segment_len_sq for stiff segments; You must, instead, square the total segment_length.
        if len(path_segments) > 1:
            SRE_lengths_and_sq = [(elem_length, elem_length**2 if stiffness > 0 else elem_len_sq)
                                  for sub_path in (path_segments[:LRE_idx], path_segments[LRE_idx+1:])
                                  for stiffness, (elem_length, elem_len_sq) in sub_path]
            SRE_len, SRE_len_sq = [sum(vals) for vals in zip(*SRE_lengths_and_sq)]
        else:
            # SRE_lengths, SRE_sq_lengths = [], []
            SRE_len, SRE_len_sq = 0, 0

    ## Minor corrections to SRE_len and SRE_len_sq:
    # SRE_len -= 0.25
    # SRE_len_sq -= 0.35
    print(" - LRE_len, SRE_len:", LRE_len, SRE_len)

    ## Check if the contour-length is long enough for the elements to "touch":
    if LRE_len > SRE_len:
        # The domains cannot reach each other. (Unless we allow helical bending which is not yet implemented)
        return 0

    # If LRE is significant, then we cannot use the simple end-end-distance < rc approximation. Instead, we must
    # consider the end-to-end PDF of the SRE part at radius r=LRE_len and within critical volume v_c.
    # Can this be approximated by a simple LRE factor?
    # LRE_factor = math.exp(-3*LRE_len_sq / (2*SRE_len_sq)) if LRE_len_sq > 0 else 1  - Uh, yes.
    # But only relevant for LRE_len_sq > SRE_len_sq/4:

    if LRE_len_sq > 0 and LRE_len_sq*4 > SRE_len_sq:
        # Domains can reach, but requires the SRE to extend beyond the mean squared end-to-end distance.
        # Use end-to-end distance PDF for the of the SRE segments using LRE as radius
        activity = (3/(2*math.pi*SRE_len_sq))**gamma * math.exp(-3*LRE_len_sq/(2*SRE_len_sq)) * AVOGADRO_VOLUME_NM3
    else:
        mean_sq_ee_dist = LRE_len_sq + SRE_len_sq       # unit of nm
        activity = (3/(2*math.pi*mean_sq_ee_dist))**gamma * AVOGADRO_VOLUME_NM3
    print(" - activity =", activity)

    ## Calculate mean end-to-end squared distance between the two domains, aka E_r_sq:
    # We already have the squared length, Nᵢbᵢ², so we just need to sum LRE and SRE:
    ## EDIT, no, wait - if we have stacked, stiff double-helices/helix-bundles, then we CANNOT
    ## just add the squared lengths together. We have to square the entire element length.
    ## We can do it for ss-helices, since the "squared length" is already just N_nt*ss_rise_per_nt*ss_kuhn_length,
    ## but that isn't valid for fully-stacked continuous segments of ds-helices.
    ## Done: Re-calculate SRE_len_sq based on the summarized segment lengths, not the path edges.

    # There shouldn't be any need to test for if mean_sq_ee_dist <= HELIX_XOVER_DIST**2 in the InterfaceGraph

    ## TODO: (Optimization) Factor out and pre-calculate (3/(2*math.pi))**gamma * AVOGADRO_VOLUME_NM3
    ##       such that activity = PREFACTOR * mean_sq_ee_dist^(-gamma)

    ## Note: "effective_volume" is just an informal term for P_loop(rc) / P_v0(rc) x AVOGADRO_VOLUME_NM3
    # effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**gamma
    # activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
    # if gamma_corr > 1:
    #     activity = activity**gamma_corr

    return activity




def loop_formation_activity(interface_graph, elem1, elem2, reaction_type):
    r"""
    Returns
        :intracomplex_activity:
    between domain1 and domain2, so that
        c_j = k_j * intracomplex_activity
    The intracomplex activity is basically just:
        activity = 1 / (N_A * effective_volume) = N_A⁻¹ * Ω⁻¹    [unit: M = mol/L]
    where NA is Avogadro's constant, 6.022e23/mol.

    """
    activity, loop_info = loop_formation_effects(interface_graph, elem1, elem2, reaction_type)
    return activity


def loop_breakage_effects(interface_graph, elem1, elem2, reaction_type):
    r"""
    Returns tuple with:
        loop_info
    Where
      * loop_info is a dict describing the loop(s) that will be changed (and, perhaps, split)
        by deleting the connection between elem1 and elem2.
    Arguments:
        :elem1:, :elem2:, :reaction_type: are the same as for intercomplex_activity().

    """
    if reaction_type is HYBRIDIZATION_INTERACTION:
        domain1, domain2 = elem1, elem2
        cmplx = domain1.strand.complex
        d1end5p, d2end5p = domain1.end5p, domain2.end5p
        d1end3p, d2end3p = domain1.end3p, domain2.end3p
        # Domains are currently NOT hybridized and thus also NOT stacked. There should be NO delegation!
        assert d1end5p.ifnode.top_delegate() == d1end5p.ifnode
        assert d1end3p.ifnode.top_delegate() == d1end3p.ifnode
        assert d2end5p.ifnode.top_delegate() == d2end5p.ifnode
        assert d2end3p.ifnode.top_delegate() == d2end3p.ifnode
        reactant_nodes = (d1end5p.ifnode, d1end3p.ifnode, d2end5p.ifnode, d2end3p.ifnode)
        ## TODO: (OPTIMIZATION) Instead of finding the shortest path, try to
        ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP.
        ## If they are, then the existing loop should give the shortest path.
        ## (For now, just assert that IF the two loops are part of an existing, non-obsolete loop,
        ## then the shortest path is in that loop.)
        elem1_ifnode, elem2_ifnode = d1end5p.ifnode, d2end5p.ifnode
        # new_loop_path = shortest_path(interface_graph, d1end5p.ifnode, d2end5p.ifnode)
        # slice_start = 1 if new_loop_path[1] == d1end3p.ifnode else 0
        # slice_end = -1 if new_loop_path[-2] == d2end3p.ifnode else None
        # pdb.set_trace()
        # if slice_start == 1 or slice_end is not None:
        #     new_loop_path = new_loop_path[slice_start:slice_end]
    elif reaction_type is STACKING_INTERACTION:
        (h1end3p, h2end5p), (h2end3p, h1end5p) = elem1, elem2
        cmplx = h1end5p.domain.strand.complex
        #  DUPLEX/DOUBLE-HELIX 1 (d: domain, h: Single helix, dh: double-helix / duplex-ends)
        #             dh1, d1         dh2, d4
        #             h1end3p         h1end5p
        # Helix 1   ----------3' : 5'----------
        # Helix 2   ----------5' : 3'----------
        #             h2end5p         h2end3p
        #             dh1, d2         dh2, d3
        elem1_ifnode, elem2_ifnode = h1end3p.ifnode.top_delegate(), h2end3p.ifnode.top_delegate()
        # new_loop_path = shortest_path(interface_graph, dh1_delegate, dh2_delegate)
    else:
        assert isinstance(elem1, DomainEnd) and isinstance(elem2, DomainEnd)
        d1end5p, d2end5p = elem1, elem2
        cmplx = d1end5p.domain.strand.complex
        elem1_ifnode, elem2_ifnode = elem1.ifnode.top_delegate(), elem2.ifnode.top_delegate()

    # We assume elem1 and elem2 has been split apart, so they should no longer be represented by a single ifnode.
    assert elem1_ifnode != elem2_ifnode

    # We assume connection is already broken between elem1 and elem2.
    if len(cmplx.loops) == 0:
        return None

    # TODO, FIX: While stacking/unstacking always creates/deletes a loop, that is not necessarily the case
    # TODO, FIX: for HYBRIDIZATION reactions.
    # TODO, FIX: Actually that is not the problem. The problem is that loops can be affected in different ways.
    # TODO, FIX: We assumed that the "shared path" was broken. But that may not be true for duplex dehybridization.
    # TODO, FIX: This should be fixed when looping over affected loops.

    # 1. Find all affected loops:
    # cmplx.ifnode_loopids_index[elem1_ifnode] => {set of loopids for loops touching ifnode}
    # Wait: If elem1 and elem2 are currently hybridized or stacked, wouldn't they have the
    # same top_delegate and thus share the exact same loops?
    # No: At this point, the breaking reaction has occured, so the two elem no longer has the same ifnode.
    # However, since we haven't updated cmplx.ifnode_loopids_index index yet, one of them should be empty
    # and the other should be full.
    if reaction_type is HYBRIDIZATION_INTERACTION:
        # A hybridized duplex (two hybridized domains) will have been represented by TWO ifnode delegates, one for each end.
        d1end5p_loops = cmplx.ifnode_loopids_index[d1end5p.ifnode]
        d1end3p_loops = cmplx.ifnode_loopids_index[d1end3p.ifnode]
        d2end5p_loops = cmplx.ifnode_loopids_index[d2end5p.ifnode]
        d2end3p_loops = cmplx.ifnode_loopids_index[d2end3p.ifnode]
        if len(d1end5p_loops) > 0:
            assert len(d2end3p_loops) == 0
            end1_loops = d1end5p_loops
            end1_former_delegate = d1end5p.ifnode
            end1_new_delegate = d2end3p.ifnode
        else:
            # assert len(d2end3p_loops) > 0  # There is no guarantee we have any loops
            assert len(d1end5p_loops) == 0
            end1_loops = d2end3p_loops
            end1_former_delegate = d2end3p.ifnode
            end1_new_delegate = d1end5p.ifnode
        if len(d2end5p_loops) > 0:
            assert len(d1end3p_loops) == 0
            end2_loops = d2end5p_loops
            end2_former_delegate = d2end5p.ifnode
            end2_new_delegate = d1end3p.ifnode
        else:
            # assert len(d1end3p_loops) > 0   # There is no guarantee we have any loops
            assert len(d2end5p_loops) == 0
            end2_loops = d1end3p_loops
            end2_former_delegate = d1end3p.ifnode
            end2_new_delegate = d2end5p.ifnode
        previous_top_delegates = ends_former_top_delegates = (end1_former_delegate, end2_former_delegate)
        reactant_node_specs = {end1_former_delegate.state_fingerprint(), end2_former_delegate.state_fingerprint()}

        newly_undelegated_ifnodes = ends_new_delegates = (end1_new_delegate, end2_new_delegate)
        affected_loopids = set.union(*[cmplx.ifnode_loopids_index[ifnode] for ifnode in reactant_nodes])
    elif reaction_type is STACKING_INTERACTION:
        # Two stacked duplex ends are represented by just ONE ifnode delegate.
        # So, for STACKING (and probably backbone ligation/nicking) we can do it simpler:
        elem1_loops = cmplx.ifnode_loopids_index[elem1_ifnode]
        elem2_loops = cmplx.ifnode_loopids_index[elem2_ifnode]
        if len(elem1_loops) > 0 and len(elem2_loops) > 0:
            # have we considered this? This probably shouldn't happen when unstacking: The two ifnodes should
            # have been represented by just ONE ifnode delegatee.
            assert False
        elif len(elem1_loops) > 0:
            previous_top_delegate = elem1_ifnode    # Or maybe there just wasn't any loops?
            previous_top_delegates = [elem1_ifnode]    # Or maybe there just wasn't any loops?
            newly_undelegated_ifnode = elem2_ifnode
            newly_undelegated_ifnodes = [elem2_ifnode]
            assert len(elem2_loops) == 0
            affected_loopids = elem1_loops
        else:
            previous_top_delegate = elem2_ifnode
            newly_undelegated_ifnode = elem1_ifnode
            previous_top_delegates = [elem2_ifnode]
            newly_undelegated_ifnodes = [elem1_ifnode]
            assert len(elem1_loops) == 0
            affected_loopids = elem2_loops
        reactant_node_specs = {previous_top_delegate.state_fingerprint()}
    else:
        # For backbone ligation, there is NO change in ifnode delegation. We do NOT merge backbone-connected ifnodes.
        raise ValueError("Only %s and %s interactions are supported!" % (HYBRIDIZATION_INTERACTION, STACKING_INTERACTION))
    if len(affected_loopids) == 0:
        return None
    # Make sure not to modify affected_loopids - we use that later! (This is just for debugging)
    # unprocesses_affected_loopids = affected_loopids.copy()
    intact_path_loops = set()
    broken_path_loops = set()

    # 2. Update all loops to include the ifnode split.
    # Edit: TODO: Should this be described genericly? Otherwise, maybe move this function to Complex
    #       to indicate that this is actually modifying the loops, not just "calculating loop effects" !
    #       OTOH, I'm making use of interface_graph which is not available within Complex.
    # TODO: Is it possible to do this as an OPTIONAL operation and perhaps split it to a separate method?
    # - Pro: It would be convenient if this method was idempotent and only had side-effects if explicitly requested.
    # - Con: Having paths reflecting current ifnode delegation makes it easier to query path edges' length.
    # - Pro: We could maybe simulate these changes like we do in loop_formation_effects where we simulate the formation reaction.

    # We need to have the updated ifnodes in the paths because we query the edge length in
    # find_alternative_shortest_path using interface_graph[source][target]['len_contour']
    affected_loop_tuples = sorted([(cmplx.loops[loopid]['activity'], cmplx.loops[loopid]['loop_hash'],
                                    cmplx.loops[loopid]['path_spec'], cmplx.loops[loopid]['path'], loopid)
                                   for loopid in affected_loopids], reverse=True)
    # affected_loops = [cmplx.loops[loopid] for loopid in changed_loopids_lst]
    # for loopid, affected_loop in zip(changed_loopids_lst, affected_loops):
    print("Pre-processing affected loop paths:")
    print("\n".join(str(tup) for tup in affected_loop_tuples))
    for activity, loop_old_hash, path_spec, path, loopid in affected_loop_tuples:
        # path = affected_loop['path']
        cmplx.loops[loopid]['debug_info']['path_before_breakage_preprocessing'] = tuple(path)  # TODO: Remove debug
        path_node_index = {ifnode: idx for idx, ifnode in enumerate(path)}
        # TODO, FIX: For duplex dehybridization, we should perhaps pre-process the path such that if
        # the previous_top_delegate ifnode is no longer connected to the next ifnode in the path,
        # then check if the newly_delegated_ifnode is and if so, swap them out.
        # Then, if both duplex ends are still connected to the rest of the path, then
        # no further path changes are necessary. You can add effects caused by changes in contour length,
        # so perhaps re-calculate loop activity, but you don't need to change the path ifnodes.
        if reaction_type is HYBRIDIZATION_INTERACTION:
            ends_path_idxs = [path_node_index.get(ifnode) for ifnode in ends_former_top_delegates]
            end_on_path = [idx is not None for idx in ends_path_idxs]
            assert any(end_on_path) # path should go through duplex somehow.
            if all(end_on_path):
                # Both duplex ends are part of the path. Figure out which is up/down stream in the path's direction
                assert abs(ends_path_idxs[0] - ends_path_idxs[1]) == 1
                # path_is_parallel = ends_path_idxs[0] < ends_path_idxs[1]
                end_downstream_idx, end_upstream_idx = ((ends_path_idxs[0], ends_path_idxs[1])
                                                        if ends_path_idxs[0] < ends_path_idxs[1] else
                                                        (ends_path_idxs[1], ends_path_idxs[0]))
                neighbor_idxs = (end_downstream_idx-1, # OK, even if idx=0
                                 (end_upstream_idx+1) % len(path))  # modulu list len to wrap around
                neighbors = [path[idx] for idx in neighbor_idxs]
                for (neighbor, top_delegate, new_delegate, idx) in zip(
                        neighbors, ends_former_top_delegates, ends_new_delegates, ends_path_idxs):
                    if top_delegate not in interface_graph.adj[neighbor]:
                        # I can't think of any situations where neither delegates are connected in the path:
                        assert new_delegate in interface_graph.adj[neighbor]
                        # Swap ifnodes:
                        print("Swapping old top delegate %s with new delegate %s on loop %s path %s" %
                              (top_delegate, new_delegate, loopid, path))
                        path[idx] = new_delegate
                    else:
                        pdb.set_trace()
                # Check if the two path-connected delegates are connected to each other.
                # TODO: This check belongs in the "find-new-loop-path-maybe" loop below!
                # Path is updated to use the best-connected ifnode.
                path_connected_ifnodes = [path[idx] for idx in ends_path_idxs]
                if path_connected_ifnodes[0] in interface_graph.adj[path_connected_ifnodes[1]]:
                    print("Loop %s: Path-connected Ifnodes %s are still connected to each other" % (
                        loopid, path_connected_ifnodes))
                    intact_path_loops.add(loopid)
                else:
                    broken_path_loops.add(loopid)
            else:
                # If not all(ends_on_path), then I don't think we need to do any pre-processing.
                # The loop path must surely change, it is just a question of how.
                # (It could only happen if the end is stacked, but then it couldn't dehybridize..
                # Also, we don't allow duplexes to be directly
                ifnode_idx, former_top_delegate = next((idx, ifnode) for idx, ifnode in zip(
                    ends_path_idxs, ends_former_top_delegates) if idx is not None)
                neighbor_before, neighbor_after = path[ifnode_idx-1], path[(ifnode_idx+1) % len(path)]
                if former_top_delegate in interface_graph.adj[neighbor_before]:
                    assert neighbor_before in interface_graph.adj[former_top_delegate]
                    assert neighbor_after not in interface_graph.adj[former_top_delegate]
                    assert former_top_delegate not in interface_graph.adj[neighbor_after]
                else:
                    assert neighbor_after in interface_graph.adj[former_top_delegate]
                    assert former_top_delegate in interface_graph.adj[neighbor_after]
                    assert neighbor_before not in interface_graph.adj[former_top_delegate]
                    assert former_top_delegate not in interface_graph.adj[neighbor_before]
                broken_path_loops.add(loopid)


        else:
            # STACKING_INTERACTION
            # TODO: If unstacking backbone-connected domain-ends, then after ifnode-split path updating below,
            # the path is essentially intact. True, we will also find the path with find_alternative_shortest_path,
            # but that search isn't needed. Instead, just check if after processing the ifnode split, the
            # two ifnodes are (a) still connected to each other, and (b) still connected to the path.
            # if so, add loop to intact_path_loops.
            if len(path) == 1:
                # Special case, backbone-connected stacked duplex ends;
                path.append(newly_undelegated_ifnode)
                print(("Inserted newly_undelegated_ifnode %s AFTER previous_top_delegate %s in loop %s path") % (
                    newly_undelegated_ifnode, previous_top_delegate, loopid))
                broken_path_loops.add(loopid)
                print(" - loop %s = %s marked as broken (length 1)." % (loopid, path))
                continue
            ifnode_idx = path.index(previous_top_delegate)
            ifnode_before = path[ifnode_idx-1] # = (ifnode_idx % len(path))-1
            ifnode_after = path[(ifnode_idx+1) % len(path)]
            if len(path) == 2:
                # For path of length 2, we should still be connected..  A--X --> A-B--X <=> A--X--B
                # We can still short-circuit the "where to insert newly_undelegated_ifnode (B)" check, but
                # we must still do the following "are A and B still connected or is loop broken" test below.
                path.append(newly_undelegated_ifnode)
                print(("Appending newly_undelegated_ifnode %s AFTER previous_top_delegate %s at end of loop %s."
                       " path: %s ") % (newly_undelegated_ifnode, previous_top_delegate, loopid, path))
            elif previous_top_delegate in interface_graph.adj[ifnode_before]:
                # former_top_delegate is still connected to ifnode_before
                # newly_undelegated_ifnode should be inserted *after* former_top_delegate in the path.
                # Not sure about these assertions, verify if you encounter AssertionErrors:
                # One reason why it may be correct is that ifnode_idx is not at the start or end of the path,
                # so the path must be 3 or more elements long.
                assert ifnode_before in interface_graph.adj[previous_top_delegate]
                assert previous_top_delegate not in interface_graph.adj[ifnode_after]
                assert ifnode_after not in interface_graph.adj[previous_top_delegate]
                assert newly_undelegated_ifnode in interface_graph.adj[ifnode_after]
                assert newly_undelegated_ifnode not in interface_graph.adj[ifnode_before]
                # TODO: Remove the excessive assertions above.
                # Insert newly_undelegated_ifnode AFTER previous_top_delegate in the path at position idx+1:
                if ifnode_idx+1 == len(path): # former_top_delegate at end of list; we can just append.
                    path.append(newly_undelegated_ifnode)
                    print(("Appending newly_undelegated_ifnode %s AFTER previous_top_delegate %s at end of loop %s."
                           " path: %s ") % (newly_undelegated_ifnode, previous_top_delegate, loopid, path))
                else:
                    path.insert(ifnode_idx+1, newly_undelegated_ifnode)
                    print(("Inserted newly_undelegated_ifnode %s AFTER previous_top_delegate %s in loop %s path at "\
                           "position %s+1") % (newly_undelegated_ifnode, previous_top_delegate, loopid, ifnode_idx))
                # Edit: Shouldn't we ALWAYS make sure the two ifnodes are at opposite ends of the path, just to make sure??
                # Well: A loop is a loops and the path is cyclic, so it should never matter where the path stats and ends.
            else:
                # former_top_delegate is NO longer connected to ifnode_before (but should be to ifnode_after)
                # newly_undelegated_ifnode should be connected to ifnode_before
                # newly_undelegated_ifnode should be inserted *before* former_top_delegate in the path.
                assert ifnode_before not in interface_graph.adj[previous_top_delegate]
                assert previous_top_delegate in interface_graph.adj[ifnode_after]
                assert ifnode_after in interface_graph.adj[previous_top_delegate]
                assert newly_undelegated_ifnode not in interface_graph.adj[ifnode_after]
                assert newly_undelegated_ifnode in interface_graph.adj[ifnode_before]
                # TODO: Remove the excessive assertions above.
                # Insert newly_undelegated_ifnode *before* previous_top_delegate in the path:
                if ifnode_idx == 0:
                    # previous_top_delegate is at start of path; newly_undelegated_ifnode before = at end.
                    path.append(newly_undelegated_ifnode)
                    print(("Appending newly_undelegated_ifnode %s BEFORE previous_top_delegate %s at end of loop %s."
                           " path: %s ") % (newly_undelegated_ifnode, previous_top_delegate, loopid, path))
                else:
                    path.insert(ifnode_idx, newly_undelegated_ifnode)
                    print(("Inserted newly_undelegated_ifnode %s BEFORE previous_top_delegate %s in loop %s path at "\
                           "position %s yielding the path: %s") % (
                        newly_undelegated_ifnode, previous_top_delegate, loopid, ifnode_idx, path))
            # Check if loop path is broken:
            if newly_undelegated_ifnode in interface_graph.adj[previous_top_delegate]:
                # Loop path still valid, re-use:
                assert previous_top_delegate in interface_graph.adj[newly_undelegated_ifnode]
                intact_path_loops.add(loopid)
                print(" - loop %s = %s marked as INTACT." % (loopid, path))
            else:
                broken_path_loops.add(loopid)
                print(" - loop %s = %s marked as BROKEN." % (loopid, path))




    # 3. Find the loop with the highest activity. This is the one we want to remove.
    # Concern: Do we need the path ifnode-split processing of this loop? Well, yeah, because we use it
    # as alternative path when piecing together new paths for broken loops.
    # Concern: Do we really need a deque? Or can we just use an iterable? Do we need to append new loops to the
    # list of affected loops? I don't think so.
    # del_loop_id = changed_loopids_deque.popleft()
    # del_loop_hash = cmplx.loops[del_loop_id]['loop_hash']
    affected_loop_tups_iter = iter(affected_loop_tuples)
    activity, del_loop_hash, del_loop_path_spec, del_loop_path, del_loop_id = next(affected_loop_tups_iter)

    # changed_loops = {del_loop_id: None}
    changed_loops_by_hash = {}
    changed_loops_hash_by_id = {}
    alternative_loop_paths = [del_loop_path]
    # reaction loop_effects directive = the collection of loop changes;
    # loop_update info = how to update each affected loop
    loop_effects = {
        'is_forming': False,
        'del_loop_id': del_loop_id,
        'del_loop_hash': del_loop_hash,
        'del_loop_path_spec': del_loop_path_spec, # for debugging
        'del_loop_activity': activity,
        # TODO: Combine all 'del_loop_xxx' into a single entry, similar to 'new_loop' in loop_formation_effects()
        'del_loop': {
            'loop_id': del_loop_id,
            'loop_hash': del_loop_hash,
            'path_spec': del_loop_path_spec, # for debugging
            'activity': activity,
            'del_loop_debug_info': {
                'del_path_tuple': tuple(del_loop_path), # Use immutable tuples to avoid mutation bugs
                'description': "loop_breakage_effects delete this loop"
            }
        },
        # 'changed_loops': changed_loops,
        'changed_loops_by_hash': changed_loops_by_hash,
        'reactant_node_specs': reactant_node_specs,
        'newly_undelegated_ifnodes_specs': {ifnode.state_fingerprint() for ifnode in newly_undelegated_ifnodes},
        'previous_top_delegates_specs': {ifnode.state_fingerprint() for ifnode in previous_top_delegates},
        # It is nice to permanently save complex state fingerprints to ensure that the complex
        # that tries to apply this loop_effects change directive has the same state
        'complex_loop_ensemble_fingerprint': cmplx.loop_ensemble_fingerprint,
        'complex_state_fingerprint': cmplx.state_fingerprint(),
        'directive_debug_info': {
            # Instance-specific debug info and other unused values
            'changed_loops_hash_by_id': changed_loops_hash_by_id,
            # 'previous_top_delegate': previous_top_delegate,         # Do not use, for debugging only
            # 'newly_undelegated_ifnode': newly_undelegated_ifnode,   # Do not use, for debugging only
            # 'previous_top_delegate_spec': previous_top_delegate.state_fingerprint(),
        },
        'description': "loop_breakage_effects directive",
    }

    # unprocesses_affected_loopids.remove(loop1_id) # changed_loopids

    # System:   Λ₀           .------.
    #            .-------,--´------. \
    #           /  Λa   /          \  \ c
    #         a \   eᵤᵥ/:/   Λb    /   \  Λc
    #            \      /       b /     \
    #             `-.--´---------´      |
    #                `-----------------´
    # We have three affected loops: Λa (aᵢ+e₁ᵢ), Λb (bᵢ+e₂ᵢ), and Λc (cᵢ+e₃ᵢ).
    # Note that the edge designation is a little weird, since it is basically a matrix:
    # Λa vs Λb:  Λa=a₂+e₂₁, Λb=b₁+e₁₂,  e₁₂==e₂₁
    # Λa vs Λc:  Λa=a₃+e₃₁, Λc=c₁+e₁₃,  e₁₃==e₃₁
    # Λb vs Λc:  Λb=b₃+e₃₂, Λc=c₂+e₂₃,  e₂₃==e₃₂

    # That is, the "shared subpath" depends on what loops we are comparing.
    # If we assume that two loops will only ever have one shared sub-path,
    # then any shared path should be eᵤᵥ, the path that is no longer usable.

    # All loops goes through eᵤᵥ which is broken apart, so all loops are affected.
    # We must: (1) Delete one loop, and (2) find new paths for the remaining loops.
    # When finding new loops, care must be taken to make the new loops (a) unique and
    # (b) use the shortest-path available.
    # In the case above:
    # affected_loopids = {Λa, Λb, Λc}
    # del_loop = Λa
    # changed_loopids_deque = (Λb, Λc)
    # alternative_loop_paths = [Λa]
    # Expected (if both Λb and Λc opts to use sub-path a for the new shortest path)
    #   Λb:  b₁+e₁₂ --> b₁+a₂
    #   Λc:  c₁+e₁₃ --> a₃+c₁
    # while loop 1:
    #   # loop1 is the loop before
    #   loop1 = Λb
    #   path2_


    # Propagate the change iteratively same way as you determine loop formation efffects:
    # while changed_loopids_deque:
    #     # For the first round, the loop considered is Λ₁ with loop_id None.
    #     # It is not appropriate to use "loop0" to describe this; you can use "loop1" or "new_loop"
    #     loop1_id = changed_loopids_deque.popleft()
    for old_activity, loop_old_hash, loop_old_path_spec, path, loopid in affected_loop_tups_iter:
        # looping over iter because we don't want to include the first to-be-deleted loop.
        if loopid in intact_path_loops:
            # Just need to update loop's activity using the updated loop path:
            path_is_broken = False
            new_loop_path = path # We can use old path (but this should've been updated with new ifnodes)
            # Note: calculate_loop_activity expects the path to START AND END on the same ifnode!
            # I.e. the edges considered are zip(path, path[1:]), it does NOT automatically include the
            # edge from path[-1] to path[0], unless circularize_path=True
            a2 = calculate_loop_activity(interface_graph, new_loop_path, simulate_reaction=None, circularize_path=True)
            if a2 == old_activity:
                print("\n\nUNEXPECTED: a2 == old_activity: %s vs %s" % (a2, old_activity))
                pdb.set_trace()
        else:
            # loop1_info = cmplx.loops[loopid]
            # loop_old_hash = loop1_info['loop_hash']
            # path = loop1_info['path']
            # loop0_path is the alternative_loop_path used to create the new_loop_path.
            assert loopid in broken_path_loops
            path_is_broken = True
            new_path_length, new_loop_path, e1, e2, e3, loop0_path = find_alternative_shortest_path(
                interface_graph, path, alternative_loop_paths)
            a2 = calculate_loop_activity(interface_graph, new_loop_path, simulate_reaction=None, circularize_path=True)
        # [ifnode.state_fingerprint() for ifnode in new_loop_path]
        new_loop_path_spec = cmplx.calculate_loop_path_spec(new_loop_path)
        new_loop_hash = cmplx.calculate_loop_hash(new_loop_path_spec)
        # new_loop_path_spec = tuple(new_loop_path_spec) # TODO: Remove cast
        activity_change_ratio = a2/old_activity
        loop_update_info = {
            'loop_hash': new_loop_hash,
            'path_spec': tuple(new_loop_path_spec), # We need the path-spec to re-create re-built loop
            'activity': a2,
            'loop_entropy': ln(a2), # loop_entropy in units of R
            'activity_change_ratio': activity_change_ratio, # = a2/a0
            'ddS': ln(activity_change_ratio),  # change in loop entropy
            # Values for the original loop (mostly for debugging)
            # 'a0': old_activity,
            # debug info for the last loop change:
            'loop_update_debug_info': {
                'path_is_broken': path_is_broken,
                'description': 'loop_breakage_effects() %s loop' % ("rebuilt" if path_is_broken else "modified"),
                'new_path_tuple': tuple(new_loop_path), # Only for debugging; the result is cached and uses between complexes
                'new_path_spec': tuple(new_loop_path_spec),
                'old_activity': old_activity,
                'old_loop_hash': loop_old_hash, # so we can relate this to the original/parent loop.
                'old_path_spec': loop_old_path_spec,
                'old_path_tuple': tuple(path),
                'source_loop_id': loopid,
            },
        }
        # changed_loops[loopid] = loop_update_info
        changed_loops_by_hash[loop_old_hash] = loop_update_info
        changed_loops_hash_by_id[loopid] = loop_old_hash
        # unprocesses_affected_loopids.remove(loop1_id)
    # assert len(unprocesses_affected_loopids) == 0  # TODO: Remove debug check
    loop_effects['total_loop_activity_change'] = (
        prod(loop['activity_change_ratio'] for loop in loop_effects['changed_loops_by_hash'].values())
        / loop_effects['del_loop']['activity']) # divide not multiply since we are DELETING not forming the loop.
    loop_effects['total_loop_entropy_change'] = ln(loop_effects['total_loop_activity_change'])
    return loop_effects


def find_alternative_shortest_path(interface_graph, loop1_path, alternative_loop_paths, avoid=None):
    r"""
    Again,
        loop1 is the "old" loop that no longer works because it is broken.
        alternative_loops are other loop paths that are also broken, but should be usable...
    The hypothesis is currently that we should be able to piece together a new shortest path
    using the previous (broken-up and processed) paths.
    ## Consider the system:   .------.
    ##            .-------,--´------. \
    ##           /  Λa   /          \  \ c = e1
    ##         a \   eᵤᵥ/:/   Λb    /   \  Λc = loop1
    ##            \      /       b /     \
    ##             `-.--´---------´      |
    ##                `-----------------´
    We are considering loop1_path = Λc = c+eᵤᵥ
    against alternative_loop_paths = [Λa, Λb]     "loop0"
    for each loop0 in alternative_loop_paths
    we find the shared part eᵤᵥ aka e3 and the not-shared-but-possibly-new-shorter sub-path e2.
    Thus, for loop1_path=Λc and alternative_loop_path=Λa, we have the following loop aliases:
        e3 = eᵤᵥ = e₁₃ = e₃₁
        e2 = a
        e1 = c
        loop2_path = e1 + e2 = c + a
    # (Yes, there was lots of changes to the nomenclature, but I think it is fixed now:)
    # EDIT EDIT EDIT: e3 is the shared part, e1 is the "unique to loop1" part, e2 is the "unique to loop0" part.

    # Note: When forming, e3 is the "possibly shorter" path and e1 the shared path part;
    # when breaking, e3 is the "broken part" (shared) and e2 is the "possibly shorter" part.
    Previous considerations:
        # Or should we make it the same: e3 is the "possibly shorter" path, e1 is the shared (but broken) path.
        # e3 is "possibly shorter", e1 is the shared part.
        # edit edit: e1 is the shared part, e3 is the "unique to loop1" part, e2 is the "unique to loop0" part.
        # edit edit edit: e3 is the shared part, e1 is the "unique to loop1" part, e2 is the "unique to loop0" part.

    Discussion: Where do the "intersection nodes" go?
    - in find_maybe_shorter_e1_e3_subpaths() we added the intersection nodes to the two candidate sub-paths
        because we wanted to compare the two candidate sub-paths directly
    We don't have that concern here: we are comparing the full loop path length, since the different
    alternatives have different shared/overlapping sub-paths.


    """
    # TODO, WARNING: THERE IS A RISK THAT LOOPS PATHS WILL BE DEGENERATE, I.E. TWO LOOPS WILL EFFECTIVELY SHARE
    # THE SAME IFNODES, BUT THE LOOPS WILL BE OF DIFFERENT ORIGIN, SOMETHING LIKE INNER VS OUTER.
    # This is e.g. the case with cross-stacked fourway junction.
    alternatives = []
    loop1_nodes = set(loop1_path)
    for loop0_path in alternative_loop_paths:
        loop0_nodes = set(loop0_path)
        e2_nodes = loop0_nodes - loop1_nodes # set(e3)
        e3_nodes = loop0_nodes & loop1_nodes # shared nodes
        # This is probably an easier and more reliable way:
        loop1_e3 = [ifnode for ifnode in loop1_path if ifnode in e3_nodes]
        loop0_e3 = [ifnode for ifnode in loop0_path if ifnode in e3_nodes]
        loop1_e1 = [ifnode for ifnode in loop1_path if ifnode not in e3_nodes]
        loop0_e2 = [ifnode for ifnode in loop0_path if ifnode not in e3_nodes]
        if loop1_nodes == loop0_nodes:
            # just return whatever is intact? Reasoning: e1 and e2 will both be empty.
            pdb.set_trace()
            loop2_path = loop1_e3[0:1]+loop1_e3[-1:]
            print("\n\nSpecial: find_alternative_shortest_path - using loop2_path", loop2_path)
            pprint(locals())
        else:
            if loop1_nodes <= loop0_nodes:
                # loop1 is covered by loop0..
                alternatives.append((path2_length, loop2_path, loop1_e1, loop0_e2, loop1_e3, loop0_path))
            if loop1_e3 == loop0_e3:
                e3_is_parallel = True
            elif loop1_e3 == loop0_e3[::-1]:
                e3_is_parallel = False
            else:
                # Have to do a more detailed check, one of the paths may start on the middle of e3.
                # This is why we usually do "e1_or_e3_subpaths" grouping and re-assembly.
                # There may also be special cases for short paths of 1-3 nodes where loops
                pdb.set_trace()
                e3_is_parallel = False
                # assert loop1_e3 == loop0_e3[::-1]
            # intersection_nodes = (loop1_e3[0], loop1_e3[-1])
            assert len(loop1_e3) >= 2
            assert loop1_e3[0] != loop1_e3[-1]
            loop1_e3_idxs = (loop1_path.index(loop1_e3[0]), loop1_path.index(loop1_e3[-1]))
            loop1_e3_idxs = (loop1_path.index(loop0_e3[0]), loop1_path.index(loop0_e3[-1]))
            if e3_is_parallel:
                # We need to reverse loop0_e2 to make it fit.
                loop2_path = loop1_e1 + loop1_e3[0:1] + loop0_e2[::-1] + loop1_e3[-1:]
                # It is generally faster to extend+append instead of concatenate, but no difference for small lists.
            else:
                loop2_path = loop1_e1 + loop1_e3[0:1] + loop0_e2 + loop1_e3[-1:]
        path2_length = sum(interface_graph[source][target]['len_contour']
                           for source, target in zip(loop2_path, loop2_path[1:]))
        alternatives.append((path2_length, loop2_path, loop1_e1, loop0_e2, loop1_e3, loop0_path))


        #
        # ## Old groupby approach:
        # # e3 is the shared part, e1 is the "unique to loop1" part, e2 is the "unique to loop0" part.
        # e1_or_e3_sub_paths = [(is_shared, list(group_iter)) for is_shared, group_iter in
        #                       groupby(loop1_path, lambda ifnode: ifnode not in loop0_nodes)] # pylint: disable=W0640
        # assert len(e1_or_e3_sub_paths) <= 3
        # if len(e1_or_e3_sub_paths) == 1:
        #     # We only have one sub-path: make sure it is e1
        #     assert e1_or_e3_sub_paths[0][0] is True  # verify that e₁ is *not* on Λ₀
        #     e1 = e1_or_e3_sub_paths[0][1]
        #     e3a, e3b = [e1[0]], [e1[-1]]
        # elif len(e1_or_e3_sub_paths) == 2:
        #     if e1_or_e3_sub_paths[0][0] is True:
        #         # We have e3b > 0; e3a starts at e1 (has zero length)
        #         e1 = e1_or_e3_sub_paths[0][1]
        #         e3a, e3b = [e1[0]], [e1[-1]]
        #         e3b.extend(e1_or_e3_sub_paths[1][1])
        #     else:
        #         # We have e3a > 0; e3b starts and ends on e1 (has zero length)
        #         e3a, e1 = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1]
        #         e3a, e3b = [e1[0]], [e1[-1]]
        #         e3a.append(e1[0])
        # else:
        #     # We cannot know whether loop1_path starts on e1 or e3, but we can easily check:
        #     if e1_or_e3_sub_paths[0][0] is True: # True means "on e1, not on e3".
        #         # loop1_path starts on e1
        #         assert tuple(tup[0] for tup in e1_or_e3_sub_paths) == (True, False, True)
        #         e1a, e3, e1b = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1], e1_or_e3_sub_paths[2][1]
        #         e1 = e1b + e1a
        #         e3 = [e1[-1]] + e3 + [e1[0]] # e3 starts with the last node in e1 and ends with the first e1 node
        #     else:
        #         # loop1_path starts on e3:
        #         assert tuple(tup[0] for tup in e1_or_e3_sub_paths) == (False, True, False)
        #         e3a, e1, e3b = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1], e1_or_e3_sub_paths[2][1]
        #         e3 = [e1[-1]] + e3b + e3a + [e1[0]] # e3 starts with last node in e1 and end in the first e1 node
        #
        # ## Unlike when forming new loops, here we are not very interested in e1.
        # ## However, we *are* interested in e2 and the total loop path length.
        #
        # # Find e2:
        # # loop1 and loop0 can potentially go in opposite directions. We need to be aware of that:
        # e3_is_parallel = loop0_path.index(e3[0]) < loop0_path.index(e3[-1])
        # e2_or_e3 = [(is_shared, list(group_iter)) for is_shared, group_iter
        #             in groupby(loop0_path, lambda node: node in e2_nodes)]  # pylint: disable=W0640
        # # e2_subpaths = [sub_path_group[1] for sub_path_group in e2_or_e3 if sub_path_group[0]]
        # assert 1 <= len(e2_or_e3) <= 3
        # if len(e2_or_e3) == 3:
        #     if e2_or_e3[0][0] is True: # we start in e2
        #         e2a, e3s, e2b = e2_or_e3[0][1] + e2_or_e3[1][1] + e2_or_e3[2][1]
        #         e2 = e2b + e2a
        #     else:
        #         e3a_, e2, e3b_ = e2_or_e3[0][1] + e2_or_e3[1][1] + e2_or_e3[2][1]
        #         e3_ = e3b_ + e3a_  # Using underscore to distinguish from e3 calculated using loop1_path
        # elif len(e2_or_e3) == 2:
        #     if e2_or_e3[0][0] is True: # we start in e2
        #         e2, e3_ = e2_or_e3[0][1] + e2_or_e3[1][1]
        #     else:
        #         e3_, e2 = e2_or_e3[0][1] + e2_or_e3[1][1]
        # else:
        #     # I don't think this should happen since e3 includes the "intersection nodes" (as does e1).
        #     assert len(e3) == 0
        #     assert e2_or_e3[0][0] is True # we start in e2
        #     e2 = e2_or_e3[0][1]
        # # alternatively, you could have looked at whether loop0_path[0] and loop0_path[-1]
        # # are both in e2_nodes
        #
        # ## Determine how e2 and e1 should be combined to form a new shortest candidate:
        # # Note: e1 and e3 e3a/e3b share the "intersection nodes", but e2 does not.
        # # Again: e3 is the shared/broken sub-path, e1 is "unique to loop1" sub-path, e2 is "unique to loop0".
        # if e3_is_parallel:
        #     # we must reverse e2 when concatenating with e1:
        #     assert e2[-1] in interface_graph.adj[e1[-1]] # the last nodes are next to each other
        #     assert e2[0] in interface_graph.adj[e1[0]]   # and the first nodes are next to each other
        #     loop2_path = e1 + e2[::-1]                        # so definitely need to reverse
        # else:
        #     assert e2[0] in interface_graph.adj[e1[-1]]
        #     assert e2[-1] in interface_graph.adj[e1[0]]
        #     loop2_path = e1 + e2
        # # loop2_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
        # #                 in group_interfaces_path_by_stiffness(loop2_path)]
        # # We only evaluate shortest paths based on length, we do not calculate activity or anything like that:
        # # path_source_target_eattr = ((source, target, interface_graph[source][target])
        # #                             for source, target in zip(loop2_path, loop2_path[1:]))
        # # # We really should have all three lengths: len_contour, dist_ee_nm, and dist_ee_sq.
        # # path_tuples = ((edge_attrs['len_contour'], edge_attrs['dist_ee_sq'], edge_attrs['stiffness'], source, target)
        # #                for source, target, edge_attrs in path_source_target_eattr)
        # path2_length = sum(interface_graph[source][target]['len_contour']
        #                    for source, target in zip(loop2_path, loop2_path[1:]))
        # alternatives.append((path2_length, loop2_path, e1, e2, e3, loop0_path))
    # end for loop0_path in alternative_loop_paths
    return min(alternatives)



def loop_formation_effects(interface_graph, elem1, elem2, reaction_type):
    r"""
    Returns tuple with:
        activity, loop_info
    Where
      * Activity is a numeric value describing the probability that the two elements are
        in close enough proximity to react, relative to the probability when the two elements
        are free in solution at a concentration of 1 M, and
      * loop_info is a dict describing the loop(s) that will be created (and, perhaps, split)
        by connecting elem1 and elem2.
    Arguments:
        :elem1:, :elem2:, :reaction_type: are the same as for intercomplex_activity().

    Making a new connection has two effects:
        1: Creating/splitting loops:
            1a: Creating a new loop.
            1b: Splitting an existing loop in two.
            1c: A split that is then merged into a single Ifnode (e.g. stacking of neighboring duplexes).
        2: Altering existing loops:
            e.g. stacking causes two helices to fuse.
    We must consider the actiity by evaluating the effect *after* performing the reaction.
    We should keep track of the effects:
        Create a new/primary loop.
        Split an existing loop in two.
        Ifnode merge.
        Joined helices.
        Other affected loops (side-effects).

    We typically want to get activity for the final state, i.e. for hybridization we want to calculate
    what the activity is with the duplex formed.
    Which strategy is better:
        (1) Perform the hybridization/stacking, calculate activity, then revert when done.
            Actually connecting the graph could make it hard to find shortest path (they are directly connected),
            but you could make other changes, e.g. adjust edge length/rigidity to what it would be when
            stacked/hybridized.
        (2) Perform calculations on unhybridized/unstacked state, but make adjustments as-needed
            to give results that are equivalent to the formed state.
            I'm currently doing this (2).

    """
    # TODO: Better handling of degenerate loops. See the "crossed fourway junction" test for an example of this.
    # TODO: For now, just forbid reactions leading to degenerate loops.

    # pdb.set_trace()

    print("\nCalculating loop_formation_effects for '%s' reaction between %s and %s" % (reaction_type, elem1, elem2))

    if reaction_type is HYBRIDIZATION_INTERACTION:
        domain1, domain2 = elem1, elem2
        cmplx = domain1.strand.complex
        d1end5p, d2end5p = domain1.end5p, domain2.end5p
        d1end3p, d2end3p = domain1.end3p, domain2.end3p
        # Domains are currently NOT hybridized and thus also NOT stacked. There should be NO delegation!
        assert d1end3p.ifnode.top_delegate() == d1end3p.ifnode
        assert d2end3p.ifnode.top_delegate() == d2end3p.ifnode
        assert d1end5p.ifnode.top_delegate() == d1end5p.ifnode
        assert d2end5p.ifnode.top_delegate() == d2end5p.ifnode
        reactant_nodes = {d1end3p.ifnode, d2end3p.ifnode, d1end5p.ifnode, d2end5p.ifnode}
        ## TODO: (OPTIMIZATION) Instead of finding the shortest path, try to
        ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP.
        ## If they are, then the existing loop should give the shortest path.
        ## (For now, just assert that IF the two loops are part of an existing, non-obsolete loop,
        ## then the shortest path is in that loop.)
        new_loop_path = shortest_path(interface_graph, d1end5p.ifnode, d2end5p.ifnode)
        # TODO: Is this really always optimial? Shouldn't we more likely EXTEND the path rather than shrink it?
        # Remember, "shortest path" is NOT automatically circularized, and we may expect
        # path[0] == path[-1] after merging ifnodes.
        #                                A-a
        #  t       5'___A________T___
        #  ,.~~---------a--- 5'      \
        #  `~------c-------- 3'      |
        #      ______________________/
        #   3'      C          B

        slice_start = 1 if new_loop_path[1] == d1end3p.ifnode else 0
        slice_end = -1 if new_loop_path[-2] == d2end3p.ifnode else None
        # pdb.set_trace()
        if slice_start == 1 or slice_end is not None:
            new_loop_path = new_loop_path[slice_start:slice_end]
    elif reaction_type is STACKING_INTERACTION:
        (h1end3p, h2end5p), (h2end3p, h1end5p) = elem1, elem2
        cmplx = h1end5p.domain.strand.complex
        # if h1end3p.domain.name in ("e", "C"):
        # pdb.set_trace()
        # Nomenclature: d: domain, h: Single helix, dh: double-helix
        # (in this case, duplex ends and double-helix ends are equivalent)
        #  DUPLEX/DOUBLE-HELIX 1
        #             dh1, d1         dh2, d4
        #             h1end3p         h1end5p
        # Helix 1   ----------3' : 5'----------
        # Helix 2   ----------5' : 3'----------
        #             h2end5p         h2end3p
        #             dh1, d2         dh2, d3
        # So, we could use alternative names:
        # dh1end3p = d1end3p = h1end3p
        # dh1end5p = d2end5p = h2end5p
        # dh2end3p = d3end3p = h2end3p
        # dh2end5p = d4end5p = h1end5p
        # and elem1, elem2 would then be:
        # assert (dh1end3p, dh1end5p), (dh2end3p, dh2end5p) == (elem1, elem2) # double-helix ends
        # assert (d1end3p, d2end5p), (d3end3p, d4end5p) == (elem1, elem2)     # domain ends

        dh1_delegate, dh2_delegate = h1end3p.ifnode.top_delegate(), h2end3p.ifnode.top_delegate()
        reactant_nodes = {dh1_delegate, dh2_delegate}
        ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP:
        assert dh1_delegate, dh2_delegate == (h2end5p.ifnode.top_delegate(), h1end5p.ifnode.top_delegate())

        new_loop_path = shortest_path(interface_graph, dh1_delegate, dh2_delegate)

        ## Adjusting new_loop_path shouldn't be required when using InterfaceGraph/Nodes:
        # slice_start = 1 if new_loop_path[1] == h2end5p else 0
        # slice_end = -1 if new_loop_path[-2] == h1end5p else None
        # if slice_start == 1 or slice_end is not None:
        #     new_loop_path = new_loop_path[slice_start:slice_end]
        # However, you MUST put the start and end nodes together as they would be when stacked,
        # (unless you choose to resolve this by increasing unstacking rate (together with throttle)..
        # - This is done AFTER getting the new_loop_path segments/elements...
    else:
        assert isinstance(elem1, DomainEnd) and isinstance(elem2, DomainEnd)
        d1end5p, d2end5p = elem1, elem2
        cmplx = d1end5p.domain.strand.complex
        new_loop_path = shortest_path(interface_graph,
                             d1end5p.ifnode.top_delegate(),
                             d2end5p.ifnode.top_delegate())
        slice_start, slice_end = 0, None


    # activity for shortest_path aka "a1" aka "new_loop_activity":
    # elem1, elem2 ifnodes have not yet been joined, so shortest path effectively already starts and ends
    # on the same ifnode, hence circularize_path=False. (This is NOT usually the case for existing loop paths!)
    activity = calculate_loop_activity(interface_graph, new_loop_path, simulate_reaction=reaction_type, circularize_path=False)
    if activity == 0:
        return 0, None

    new_loop_nodes = set(new_loop_path)
    processed_secondary_loopids = set()
    all_affected_loops = set()
    # changed_loops = {}   # old_id => [new_loop1_dict, new_loop2_dict, ...]
    changed_loop_paths = {} # # lookup of loopid: new_loop_path for all changed loops
    changed_loops_by_hash = {} # old_loop_hash: {dict with new loop info.}
    # It would be nice to also have the activities for the new loops created (not just the shortest-loop)
    # Maybe have a 'changed_loops' old_loopid => [list of dicts describing new loops]
    activity_change_ratios = []
    side_effects_factor = 1
    side_effects_activities = []
    side_effects_factors = []
    path_spec = cmplx.calculate_loop_path_spec(new_loop_path)
    ## TODO: Remove the checks below >>
    ifnode_by_hash_index_before = set(cmplx.ifnode_by_hash.items())
    cmplx.rebuild_ifnode_by_hash_index()
    cmplx.rebuild_loopid_by_hash_index()
    ifnode_by_hash_index_after = set(cmplx.ifnode_by_hash.items())
    assert path_spec == cmplx.calculate_loop_path_spec(new_loop_path)
    assert ifnode_by_hash_index_before == ifnode_by_hash_index_after
    print("\nnew loop (shortest_path):", new_loop_path)
    print(" - new_loop_path_spec:", path_spec)
    print(" - activity:", activity)
    ## TODO: Remove until here <<
    new_loop = {
        'loop_hash': cmplx.calculate_loop_hash(path_spec),
        'path_spec': tuple(path_spec),
        'activity': activity,
        'loop_entropy': ln(activity),
        'ddS': 0,  # Change in loop entropy caused by this reaction.
        'new_loop_debug_info': {
            'new_path_tuple': tuple(new_loop_path), # Use immutable tuples to avoid mutation bugs
            'description': "loop_formation_effects new loop"
        },
    }
    print(" - new_loop info: %s" % (new_loop,))
    changed_loop_paths[None] = new_loop_path

    changed_loops_hash_by_id = {}
    new_stacking_loops = {}
    # Keep track of loops that were considered by shared_loopids, but was found to be unchanged by the formation:
    touched_but_unchanged_loops = {}  # hash => path
    reactant_node_specs = {ifnode.state_fingerprint() for ifnode in reactant_nodes}
    # Note: "shared_nodes/loops" is between a particular new_changed loop and an existing, overlapping loop.
    # The "shared part" does NOT have to include the reactant ifnodes.
    # When doing ifnodes-are-merged loop path updates, it is better to go over all loops in the union of
    # set.union(cmplx.loops_by_ifnode[ifnode] for ifnode in reactant_node_specs)
    loop_effects = {
        # total loop formation activity; Product of shortest-path activity and any/all loop_change_factors:
        # 'activity': 0, # Renamed to 'total_loop_formation_activity' to avoid confusing individual loop's activity with reaction activity
        'total_loop_activity_change': activity,
        'is_forming': True,
        # 'shortest_path': tuple(new_loop_path), # use ['new_loop']['path'] instead.
        # 'shortest_path_spec': path_spec, # use ['new_loop']['path_spec'] instead.
        # 'shortest_path_activity': activity, # use ['new_loop']['activity'] instead.
        'new_loop': new_loop,
        # 'changed_loops': changed_loops, # old_loopid => [list of new loop dicts]
        'changed_loops_by_hash': changed_loops_by_hash, # old_loop_hash:  {dict with new loop info.}
        # 'shared_loop_hashes': None,  # shared/touched loops, from shared_loopids
        # 'touched_but_unchanged_loops': touched_but_unchanged_loops,
        'reactant_node_specs': reactant_node_specs,
        'directive_debug_info': {
            # Instance-specific debug info and other unused values
            'changed_loops_hash_by_id': changed_loops_hash_by_id,
            'reactant_nodes': reactant_nodes,
            'shared_loopids': None
        },
        # 'changed_loops_hash_by_id': changed_loops_hash_by_id,
        # changed_loops is a list of dicts, each representing the new loop formed, with keys 'path' and activity.
        # 'loops_considered': processed_secondary_loopids, # for debugging
        # 'loops_affected': all_affected_loops,  # Obsoleted by changed_loops
        'activity_change_ratios': activity_change_ratios,
        # loops not reached by "neighbor-loop-propagation" algorithm but which should still be considered.
        # There may be additional loops that are relevant due to stacking (network-connected ends)
        # These are collected in new_stacking_loops:
        'new_stacking_loops': new_stacking_loops,
        'stacking_side_effects_activities': side_effects_activities,
        'stacking_side_effects_factors': side_effects_factors,
        'complex_loop_ensemble_fingerprint': cmplx.loop_ensemble_fingerprint,
        'complex_state_fingerprint': cmplx.state_fingerprint(),
        'description': "loop_formation_effects directive",
    }

    if not cmplx.loops:
        # Complex has no existing loops, just return already:
        return activity, loop_effects


    ## Useful NetworkX functions:
    ## networkx.algorithms.cycles.cycle_basis

    ## TODO: I'm still contemplating if it's easier to check for side-effects if we
    ## temporarily modify the interface_graph to look as it would if the reactants were reacted
    ## and then undo the change afterwards. Right now we are simulating this in the methods:
    ## * join_edges_groups, join_two_edge_groups, and
    ## * closed_loop_segments
    ## These may be called repeatedly for each secondary loop that is checked for side-effects.
    ## ...but we are usually modifying the already-calculated shortest_path, so no.



    ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP:
    ## TODO: CHECK FOR SECONDARY LOOPS!
    ## We can easily do this simply by whether the two nodes are already in an existing loop.
    ## NOTE: If both reactant nodes are part of an existing loop, we don't strictly have to
    ## spend resources calculating the shortest path: The shortest path should be a sub-path
    ## of the existing loop. And, we also need the "other" part of the loop:
    ## We calculate reaction activity for "internal" loop node reactions as:
    ##      ab̌c = ab̌ + b̌c - ac
    ## where ab̌c is the internal loop closure, ac is closing the outer loop,
    ## and ab and bc is closing each of the inner loops.
    ## There might also be a term describing the change that e.g. hybridization
    ## would have on the existing loop (secondary side-effects).
    ## Note: The above expression is for shape energies. If we are talking activities, then
    ## we simply multiply (instead of add) and divide (instead of subtract).

    ## Q: When joining two ifnodes, are we sure that at most ONE existing loop is dividied into two?
    ## I would expect this to be the case, I can't currently imagine how joining two nodes
    ## should be able to split more than one existing loop...

    ## Q: Should cmplx.loops_by_interface contain *all* loops or only the current effective loops?
    ## A: Only the effective; determining other "outer" loops could be ambigious, I think.

    ## TODO: Short-circuit shortest_path calculation by checking if reactants are part of existing loop first
    ## (right now we are using it for assertions)

    ## TODO: There are cases when we are splitting an existing loop, but not by direct "pinching":
    ## In fact, IN MOST CASES, path[0] and path[-1] will NOT BOTH be part of an existing loop, even though
    ## they are effectively splitting a loop in two. See notes on Wang-Uhlenbeck entropy.
    # if False and (path[0] in cmplx.ifnode_loopids_index and len(cmplx.ifnode_loopids_index[path[0]]) > 0
    #     and path[-1] in cmplx.ifnode_loopids_index and len(cmplx.ifnode_loopids_index[path[-1]]) > 0):
    #     ## Uh, OK so both elems are part of one or more existing loops. But we don't know if they are
    #     ##  part of the SAME loop. Not the best way to do it.. In any case, it is probably better
    #     ##  to handle this case like any other, as it is done below.
    #     ## For the record, this is what you should have done:
    #     loop_with_both_nodes = cmplx.ifnode_loopids_index[path[0]] & cmplx.ifnode_loopids_index[path[0]]
    #     ## TODO: Need to test this - e.g. case 1(a) or 1(b)
    #     ## TODO: I feel like we should still cosider *all* shared loops and not just break out by this one case ^^
    #     ## TODO: THIS CASE NEEDS TO BE COMPLETELY RE-WORKED!
    #     ## Perfect "pinching" of existing loop (rare)
    #     # We have a secondary loop.
    #     # Find loops containing both the first and last path interface nodes:
    #     shared_loopids = cmplx.ifnode_loopids_index[path[0]] & cmplx.ifnode_loopids_index[path[-1]]
    #     assert len(shared_loopids) == 1           ## Not sure about this...
    #     loop0id = next(iter(shared_loopids))
    #     loop0 = cmplx.loops[loop0id]
    #     loop0_path = loop0['path'] # list of ifnodes from when the loop was made.
    #     # IF the reaction/node-merge splits more than ONE loop, we want to know!
    #     # Can we always use the smallest loop, or is there a particular order that must be obeyed?
    #     # Also, is the number of nodes "len(loop['nodes'])" the way to go?
    #     # It is not the same as the "shortest" path method..
    #     ## TODO: Check/fix this code:  - EDIT: Replaced by len(shared_loops) == 1 above.
    #     # loop0 = min(shared_loops, key=lambda loop: len(loop['nodes']))  # NOT GUARANTEED
    #     shortest_path_nodes = set(path)
    #     assert all(node in loop0['ifnodes'] for node in path)             # NOT GUARANTEED
    #     assert shortest_path_nodes <= loop0['ifnodes']
    #     # Using set operations is faster (path is subset of existing loop0)
    #     # Set operators: &: intersection, -: difference), ^: symmetric_difference, <=: is subset, >=: superset.
    #     path2_nodes = loop0['ifnodes'] - shortest_path_nodes
    #     loop0_path_idxs = (loop0_path.index(path[0]), loop0_path.index(path[-1]))
    #     if loop0_path_idxs[0] > loop0_path_idxs[1]:
    #         loop0_path_idxs = loop0_path_idxs[::-1]
    #     path1 = e1 = loop0_path[loop0_path_idxs[0]:loop0_path_idxs[1]+1]
    #     path2 = e2 = loop0_path[loop0_path_idxs[1]:] + loop0_path[:loop0_path_idxs[0]+1]
    #     # path is the shortest path and should equal either e1 or e2:
    #     assert all(node in path1 for node in path) or all(node in path2 for node in path)
    #     # a1 = calculate_loop_activity(path1, simulate_reaction=reaction_type)
    #     a2 = calculate_loop_activity(path2, simulate_reaction=reaction_type)
    #     a0 = loop0["activity"]
    #     ## THIS SHOULD JUST BE APPENDED TO activity_change_ratios:
    #     activity_change_ratios.append(a2/a0)
    #     processed_secondary_loopids.add(loop0id)
    #     pdb.set_trace()
    #     # If hybridization, special case where we account for the length of the newly-formed duplex.
    #     # TODO: This should not preclude processing of other shared loops, but rather be just an optimization.

    # elif cmplx.ifnode_loopids_index:

    ## TODO, WIP: Add bredth-first algorithm to search for all changed loops.
    # Algorithm:
    # Strategy:
    # * Initialize simply by the new loop changed_loops and new_loop_id=None to changed_loops_deque,
    # * Then loop over changed_loops_deque and process just like you are already doing.

    # Initialize bredth-first search for changed loops:
    changed_loops_deque = deque()
    changed_loops_deque.append(None) # None is used to temporarily identify the newly formed loop

    ## Note:
    ## After forming the new loop, Λ₁, if e3 < e2, then the old Λ₀ is updated to be the Λ₂ loop (aka Λ₀').
    ## If e3 is longer than e2, then we don't do anything. This is arguably an in-complete description
    ## of the shape entropy. But then, shape entropy should be described by the edges, not loops,
    ## since loops are not actually additive.

    ## Check for non-pinching loop-splitting cases: (partial loop sharing)
    ##           Λ₀              .------.
    ##            .----------,--´------. \
    ##           /   Λ₁     /          \  \ e₄ (not part of the original discussion)
    ##        e₁ \      e₃ /:/   Λ₂    /   \
    ##            \         /         / e₂  \
    ##             `----.--´---------´      |
    ##                   `-----------------´
    ## Where symbols Λ₀ Λ₁ Λ₂ e₁ e₂ e₃:
    ## Λ₁ = e₁+e₃ is the shortest path loop, for the reaction in question,
    ## Λ₀ = e₁+e₂ is an existing loop, overlapping the reaction shortest loop path,
    ## Λ₂ = e₂+e₃ is the path for the other loop, formed by splitting Λ₀.
    ## e₁ is the part shared by both Λ₁ and Λ₀ but not in Λ₂
    ## e₃ is the path in both Λ₁ and Λ₂ but not in Λ₀,
    ## e₂ is the part shared by both Λ₂ and Λ₀ but not in Λ1
    ## e2 > e1 because by shortest-path definition e₁+e₃ < e₂+e₃.
    ## The Wang-Uhlenbeck matrix can be written as:
    ##      ┌ e1+e3   e3  ┐   ┌ e1+e3   e1  ┐   ┌ e1+e2   e2  ┐
    ##  W = │             │ = │             │ = │             │
    ##      └  e3   e2+e3 ┘   └  e1   e1-e2 ┘   └  e2   e3+e2 ┘
    ## In general we have Wang-Uhlenbeck determinant (where dS = R α ln(det(W)))
    ##      Before forming e3: det(W₀) = Λ₀ = e₁+e₂
    ##      After  forming e3: det(W₁) = e₁e₂ + e₂e₃ + e₃e₁
    ##                         = (e1+e3) (e2+e3) - e3^2  = Λ₁ Λ₂ - e3^2
    ##                         = (e3+e1) (e2+e1) - e1^2  = Λ₁ Λ₀ - e1^2
    ## For e₃ == 0, we have the special case:
    ##  det(W₀) = e₁e₂       and det(W₀) = Λ₀ = e₁+e₂
    ##      a = a₁*a₂/a₀     where aᵢ is the activity calculated for loop Λᵢ alone.
    ## What about the case where e3 is a single, rigid element?
    ##  * This would be the stacking-reaction side-effects considered above..
    ##  * Means that Λ₁ and Λ₂ are independent. Same as for e3 == 0.
    ## Maybe all cases where Λ₁ and Λ₂ are independent can be calculated by multiplying a₁ with a₂/a₀.
    ## Can we consider the a₂/a₀ a "side-effect"?
    ## Thus, just go over all shared loops and ask:
    ## * Is e3 splitting the loop into two independent loops? - Calculate factor a₂/a₀.
    ## The criteria for e3 yielding independent loop splitting is:
    ## * If e3 is a single, rigid segment,
    ## * If e1 >> 0 and e3 < e1 (e1 < e2 by definition).
    ## For e1 ≈ 0:
    ##  det(W₁) = Λ₁ Λ₀ - e1^2 ≈ Λ₁ Λ₀
    ##  a ≈ a₁*a₀/a₀ = a₁
    ## This latter condition should also be applicable if e3 >> e2 > e1
    ## However, these case approximations assume e1, e2, e3 are made up of a big number of flexible kuhn links.
    ##
    ## Is there a way to look at these purely from a
    ##
    ## Can we view this as "the shortest-path activity modulated by the
    ## This is really non-obvious, but let's give it a try...
    ## First see if there is any overlap between the shortest path and existing loops:
    activity_change_ratio = 1
    while changed_loops_deque:
        # For the first round, the loop considered is Λ₁ with loop_id None.
        # It is not appropriate to use "loop0" to describe this; you can use "loop1" or "new_loop"
        # pdb.set_trace()
        loop1_id = changed_loops_deque.popleft()
        # loop1_info = changed_loops[loop1_id]
        # loop1_info = cmplx.loops[loop1_id]
        # loop1_shortest_path = tuple(loop1_info['path'])
        loop1_shortest_path = tuple(changed_loop_paths[loop1_id])


        # loops_with_both_nodes = cmplx.ifnode_loopids_index[path[0]] & cmplx.ifnode_loopids_index[path[-1]]
        loop1_path_nodes = set(loop1_shortest_path)
        # Quick way to get all nodes on the shortest path that are already on *ANY* existing loop.
        shared_nodes = loop1_path_nodes.intersection(cmplx.ifnode_loopids_index.keys())
        if shared_nodes:
            shared_loopids = set.union(*[
                cmplx.ifnode_loopids_index[ifnode] for ifnode in shared_nodes
            ]).difference(processed_secondary_loopids)
        else:
            shared_loopids = {}
        # assert loops_with_both_nodes <= shared_loopids
        # assert len(loops_with_both_nodes) <= 1
        ## We have a few special cases:
        ## If the new edge e₃ is completely rigid, then e₂ has no influence on Λ₁, nor does e₁ influence on Λ₂.
        ## - The "e₃ == 0" aka "loop pinching" above can be seen as a special case of this.
        ## Perhaps just simplify as:
        ## * If e₃ has fewer segments than e₁, e₂ then calculate activity as:
        ##      a = a₁*a₂/a₀, where a₁ = f(e₁+e₃)e₂, a₂ = f(e₂+e₃), and a₀ = f(e₁+e₂) is the outer loop activity.
        ##      and f is calculate_loop_activity function.
        ##   How is this different from what is done in the side-effects?
        ##   * Currently we only search for side-effects from stacking.
        ## * If e₃ has more segments than e₁, e₂ then calculate activity simply using the shortest path:
        ##      a = a₁ if e₁ < e₂ else a₂
        ## Symbols: Λ₀ Λ₁ Λ₂ e₁ e₂ e₃

        # if shared_loopids:
        # shared_loops = [ for loop0id in shared_loopids]
        ## TODO: Check that this does not overlap with the "side-effects" factors calculation..
        ## It actually seems like this could supplant the side-effects calculation.
        # Find nodes that are not part of any loops:
        # fully_unshared_nodes = shortest_path_nodes.difference(cmplx.ifnode_loopids_index.keys())
        for loop0_id in shared_loopids:
            loop0 = cmplx.loop_by_loopid[loop0_id]
            print("Processing possibly-affected loop0 %s: %s" % (loop0_id, loop0))
            # path is shortest-path loop Λ₁, shared loop is Λ₀ aka "loop0".
            # Makes lookup faster for long loop0 but with constant overhead.
            ## TODO (optimization): If path[0] and path[-1] are both in the loop, then this can be optimized,
            ##  because you don't have to group, but can simply calculate the e1, e2 paths directly.
            loop0_path = tuple(loop0['path'])
            # When stacking backbone-connected duplex-ends, we form (self-)loops with a single node and no edge.
            # Do not consider these for optimization.
            if len(loop0_path) < 2:
                continue
            # TODO: Remove loop_hash assertion:
            loop0_path_spec = cmplx.calculate_loop_path_spec(loop0_path)
            loop0_hash = cmplx.calculate_loop_hash(loop0_path_spec)
            assert loop0['loop_hash'] == loop0_hash  # TODO: Remove recalculation check of loop0_hash
            loop0_nodes = set(loop0_path)
            assert len(loop0_nodes) > 1  # If path has more than 1 elems, the set should as well.
            # Casting group_iter to list: (TODO: Check if that is really needed...)

            # Note: loop1_shortest_path (e1+e3) is *NOT* a candidate path for loop0.
            # loop1_shortest_path is used to find the overlap (e1) between loop0 and loop1 and determine if the
            # "non-overlapping subpath" (e3) is a shorter route for loop0.
            # Should probably check if e1 fully covers loop0.
            if len(loop0_nodes - loop1_path_nodes) == 0:
                print("loop1_path %s fully covers loop0 path %s, ignoring loop0..." % (loop1_shortest_path, loop0_path))
                touched_but_unchanged_loops[loop0_hash] = loop0_path
                continue
            e1, e3a, e3b, e3_groups, use_e3 = find_maybe_shorter_e1_e3_subpaths(
                interface_graph, loop0_nodes, loop1_shortest_path, reaction_type)

            # Calculate new activity for new loop0 path
            if use_e3:
                # New loops Λ₁ and Λ₂ are independent, multiply activity by factor a₁/a₀:
                # TODO: Account for cases where e3a == [] or e3b == [] or both (e3==0) !
                a0 = loop0['activity']

                # Find e2:
                e2_nodes = loop0_nodes - set(e1)
                e2_not_e1 = [(is_shared, list(group_iter)) for is_shared, group_iter
                            in groupby(loop0_path, lambda node: node in e2_nodes)]
                e2_subpaths = [sub_path_group[1] for sub_path_group in e2_not_e1 if sub_path_group[0]]
                # If Λ₀ loop path starts inside e2, we will have two sub-paths, otherwise only 1:
                if len(e2_subpaths) > 1:
                    e2 = e2_subpaths[1] + e2_subpaths[0]
                else:
                    e2 = e2_subpaths[0]
                # alternatively, you could have looked at whether loop0_path[0] and loop0_path[-1]
                # are both in e2_nodes
                # Figure out how e2 and e3 should be combined: e2+e3 or e3+e2.
                # Note: e1, e2 and e3a/e3b all share the node at the intersection, right? Edit: No, only e1, e3
                # Note: We are only concerned with the intersection nodes that are e3a[-1] and e3b[0]
                # It does not make sense to use e3a[0] or e3b[-1] for this
                # (which are the reaction pair and shortest_path start/end)
                if e2[0] in interface_graph[e3a[-1]]:
                    # Connection between first part of e3 and start of e2, and last part of e3 and end of e2:
                    assert e3a[-1] in interface_graph[e2[0]]
                    assert e2[-1] in interface_graph[e3b[0]]
                    assert e3b[0] in interface_graph[e2[-1]]
                    new_loop_path = e3a + e2 + e3b
                else:
                    # Connection between last node of first part of e3 and the last  node of e2,
                    # and between        first node of last part of e3 and the first node of e2:
                    assert e3a[-1] in interface_graph[e2[-1]] # last node of e3a + last of e2
                    assert e2[-1] in interface_graph[e3a[-1]]
                    assert e3b[0] in interface_graph[e2[0]]   # first node of e3a + first of e2
                    assert e2[0] in interface_graph[e3b[0]]
                    # Need to reverse direction of e3 or e2 before joining them:
                    # e3 is e3b extended by e3a (list performs best by extending at the end)
                    new_loop_path = e3a + e2[::-1] + e3b
                print(" - The new loop produced a shorter path for existing loop0 %s, new_loop_path = %s" % (
                    loop0_id, new_loop_path))
                # Should do we have to make sure new_loop_path starts and ends on the same ifnode by passing circularize_path=True?
                # I guess not: path starts and ends on elem1, elem2 ifnodes, which will be the same after the reaction has been done.
                # Again, this method works PREDICTIVELY before the reaction is performed.
                a2 = calculate_loop_activity(interface_graph, new_loop_path, simulate_reaction=reaction_type, circularize_path=False)
                if a2 == 0:
                    return 0, None
                if a2 == a0:
                    print("\n\nUNEXPECTED: a2 == old_activity: %s vs %s" % (a2, a0))
                    pdb.set_trace()
                activity_change_ratio = a2/a0
                ## Question: Do we need to add loop1 = e1+e3 as replacement loop?
                ## Current answer: No, loop1 will always be created as the shortest_path loop.
                ## The 'e1' and 'e2' edges may be different for different loop0s,
                ## but loop1 should always be the same (PROVIDED THAT THE SHORTEST PATH IS ON loop0)
                ## Note that this is NOT the always the case below where I'm considering loops on the full ds helix.
                # [ifnode.state_fingerprint() for ifnode in new_loop_path]
                new_loop_path_spec = cmplx.calculate_loop_path_spec(new_loop_path)
                new_loop_hash = cmplx.calculate_loop_hash(new_loop_path_spec)
                if new_loop_hash in changed_loops_by_hash:
                    print("New loop hash %s already in changed_loops_by_hash %s" %
                          (new_loop_hash, changed_loops_by_hash))
                    pdb.set_trace()
                if new_loop_hash in cmplx.loopid_by_hash:
                    print("\n\nNOTE: New loop hash %s already exists in cmplx.loopid_by_hash %s" %
                          (new_loop_hash, cmplx.loopid_by_hash))
                    # pdb.set_trace()
                # DEBUG: Test that the path_spec doesn't change if updating ifnode_by_hash index:
                cmplx.rebuild_ifnode_by_hash_index()
                cmplx.rebuild_loopid_by_hash_index()
                assert new_loop_path_spec == cmplx.calculate_loop_path_spec(new_loop_path)
                # TODO: Remove re-calculations and assertions above (or do at the end instead of in loop!)
                assert cmplx.calculate_loop_hash(new_loop_path_spec) == cmplx.calculate_loop_hash(tuple(new_loop_path_spec))
                # TODO: Remove index re-calculation and assertions above
                loop_update_info = {
                    'loop_hash': new_loop_hash, # Used to identify the loop instance
                    'path_spec': tuple(new_loop_path_spec), # Is used to re-create the loop path ifnodes
                    'activity': a2,
                    'activity_change_ratio': activity_change_ratio, # = a2/a0
                    'loop_entropy': ln(a2), # in units of R
                    'ddS': ln(activity_change_ratio), # Change in loop entropy
                    # Values for the original loop (mostly for debugging)
                    'loop_update_debug_info': {
                        'description': "loop_formation_effects updated loop found by breath-first search starting at new loop.",
                        'new_path_tuple': tuple(new_loop_path), # Only for debugging; the result is cached and uses between complexes
                        'old_activity': a0,
                        'old_loop_hash': loop0_hash, # so we can relate this to the original/parent loop.
                        'old_path_spec': loop0_path_spec,
                        'old_path_tuple': tuple(loop0_path),
                        'source_loop_id': loop0_id,
                    },
                }
                print(" - New info (path, hash, activity) for loop0 %s: %s\n" % (loop0_id, loop_update_info))
                assert loop0_id not in changed_loop_paths
                # changed_loops[loop0_id] = loop_update_info
                changed_loop_paths[loop0_id] = new_loop_path
                assert loop0_hash not in changed_loops_by_hash
                changed_loops_by_hash[loop0_hash] = loop_update_info
                changed_loops_hash_by_id[loop0_id] = loop0_hash
                activity_change_ratios.append(activity_change_ratio)
                # Add loop0_id to deque so we can check if loop0 has any affected neighboring loops:
                changed_loops_deque.append(loop0_id)
            # end if use_e3 (loop0 should be updated)
        # end for loop0_id in shared_loopids:
        ## TODO: THIS IS NOT DONE YET!
        # Union-assignment operator should modify in-place and be equilvalent to myset.update(other_set)
        processed_secondary_loopids |= shared_loopids
        processed_secondary_loopids.add(loop1_id)
    # end while changed_loops_deque


    #### Checking for further secondary side-effects: ####
    ## TODO: Make sure side-effects does not overlap with primary loop calculation
    ## Even if the reaction form a new loop, the reaction can still have structural effects
    ## that alters the energy of existing loops. The presense of these loops would equally
    ## influence the activity for forming this new loop.
    ## E.g. if stacking of two helices would require excessive stretching.
    ## Evaluate secondary effects of the reaction for existing loops:
    if reaction_type is STACKING_INTERACTION:
        # It should be enough to consider existing loops; if there are connections not part of
        # an existing loop, that connection should be part of the new loop created by this reaction.
        # Any loop with nodes on dh1/dh2:
        # What are double-helices? DomainEnds? Or IfNodes? Or an object instance with both?
        # dh1 = cmplx.helix_by_domainend[h1end3p]  # double-helix, duplexed and fully-stacked
        # dh2 = cmplx.helix_by_domainend[h2end3p]
        # dh1_ifnodes, dh2_ifnodes = ({end.ifnode for end in dh1}, {end.ifnode for end in dh2})
        # Actually, it's probably easier to calculate it as-needed using:
        # dh1_upstream_ifnodes = list(h1end3p.upstream_stacked_top_ifnodes_generator())
        # dh2_upstream_ifnodes = list(h2end3p.upstream_stacked_top_ifnodes_generator())
        # Treat the two stacking double-helices as arms, pivoting at the stacking junctiong:
        # Both are in direction inside-out away from the junction and includes the first reactant.
        dh1_upstream_ifnodes = h1end3p.upstream_stacked_top_ifnodes_list([dh1_delegate])
        dh2_upstream_ifnodes = h2end3p.upstream_stacked_top_ifnodes_list([dh2_delegate])
        # dh_ifnodes_list = (dh1_upstream_ifnodes[::-1]
        #                    + [top_ifnode for top_ifnode in
        #                       (end.ifnode.top_delegate() for end in (h1end3p, h2end3p))
        #                       if (top_ifnode != dh1_upstream_ifnodes[0] and
        #                           top_ifnode != dh2_upstream_ifnodes[0])]
        #                    + dh2_upstream_ifnodes)
        # dh_arm_idx_by_ifnode = {sign*i: ifnode
        #                         for sign, nodes in ((1, dh1_upstream_ifnodes), (-1, dh2_upstream_ifnodes))
        #                         for i, ifnode in enumerate(nodes)}
        dh1_ifnodes = set(dh1_upstream_ifnodes)
        dh2_ifnodes = set(dh2_upstream_ifnodes)
        dh_ifnodes = dh1_ifnodes | dh2_ifnodes
        dh1_nodes_in_loops = dh1_ifnodes.intersection(cmplx.ifnode_loopids_index.keys())
        dh2_nodes_in_loops = dh2_ifnodes.intersection(cmplx.ifnode_loopids_index.keys())
        # Edit: Loops are currently stored as mutable dicts, you cannot make a set of loops.
        # You must either use loop IDs or create an object for each loop.
        dh1_loopids = {loopid for ifnode in dh1_nodes_in_loops if ifnode in cmplx.ifnode_loopids_index
                       for loopid in cmplx.ifnode_loopids_index[ifnode]}
        dh2_loopids = {loopid for ifnode in dh2_nodes_in_loops if ifnode in cmplx.ifnode_loopids_index
                       for loopid in cmplx.ifnode_loopids_index[ifnode]}
        # Alternatively, go through all loops and check if any of the loop's nodes are on the doublehelix:
        # dh1_loops = {loopid for loopid in cmplx.loops.keys() if any(...) }
        ## NOTE: We should not consider loops that are only on one of the dh arms and not the other,
        ## hence the two set comprehensions which are then followed by the intersection.
        ## Affected loops are loops which goes through both double-helix 1 and 2:
        affected_loops = dh1_loopids & dh2_loopids
        all_affected_loops |= affected_loops  # equivalent to in-place update(other_set)
        ## By definition, the shortest path for loops with nodes on both dh1 and dh2 should now
        ## go through the created stacked junction.
        ## In general, a stacking interaction will cause affected loops to split into two loops.
        ## However, there is a special but frequent case where one of the new loops will simply be the
        ## phosphate backbone between two adjacent duplexes. In that case we don't have to split the loop in two,
        ## but simply re-calculate it.
        ## Loop over all affected loops and re-calculate the activity for them.
        ## We then take the product: side_effects_factor = ∏(a_new/a_old for loop in affected loops)
        ## However, we must either modify the path to ensure it goes through the stacked junction,
        ## or re-calculate the shortest path for the loop, expecting it to go through the junction.
        # Optimization:
        affected_loops -= processed_secondary_loopids
        for loop0id in affected_loops:
            print("After updating the initial search for changed loops, this should never happen anymore!")
            pdb.set_trace()
            if loop0id in processed_secondary_loopids:
                # pass
                continue
            ## TODO: Add another (optional; expensive) check to look for "network-connected helix stacking branches"

            loop0 = cmplx.loop_by_loopid[loop0id]
            # What is stored in a loop? Let's just start by letting loop be a dict with whatever info we need:
            loop0_path = loop0["path"]   # List of ifnodes
            assert loop0_path == [ifnode.top_delegate() for ifnode in loop0_path] # make sure we have top delegates
            loop0_path_spec = cmplx.calculate_loop_path_spec(loop0_path)
            loop0_hash = cmplx.calculate_loop_hash(loop0_path_spec)
            assert loop0['loop_hash'] == loop0_hash  # TODO: Remove recalculation check of loop0_hash
            loop0_nodes = set(loop0_path)
            assert len(loop0_path) == len(loop0_nodes)
            # Make sure the shortest path is not a sub-path of loop0, this should've been processed above.
            assert not new_loop_nodes <= loop0_nodes
            # new_paths = []
            a0 = loop0['activity']
            a1a2 = []
            # changed_loops[loop0id] = {
            #     'new_paths': new_paths,
            #     'a1a2': a1a2,   # Primary loop closing activity for each new path
            #     'a0': a0,       # Store this, for good measure...
            #     #'loop_pinching_activity': loop_pinching_activity  # Combined effect of splitting old loop in two
            # }
            replacement_loop1 = {}
            # list of dicts. Two dicts, unless the two ends are currently adjacent.
            # Note: Any loop with nodes present on the shortest-path should be considered in the check above.
            # If the two ends are adjacent in the interface graph, then they are certainly on the shortest path!
            if dh1_delegate in loop0["path"] and dh2_delegate in loop0["path"]: #and \
                ## TODO: MAKE SURE WE ARE NOT DOUBLE-COUNTING; C.F. CALCULATION USING SECONDARY LOOP SEARCH BELOW.
                print("WEIRD: Reaction ifnodes %s and %s are both on an existing loop, but somehow this loop "
                      "was not included when processing shared_loops??")
                # abs(loop0["path"].index(dh1_delegate) - loop0["path"].index(dh2_delegate)) == 1:
                # We can just modify the path, joining it together at the new junction:
                idx1, idx2 = loop0["path"].index(dh1_delegate), loop0["path"].index(dh2_delegate)
                if idx1 > idx2:
                    idx1, idx2 = idx2, idx1
                new_path1 = loop0["path"][idx2:] + loop0["path"][:idx1+1]  # +1 to include the node.
                # new_paths.append(new_path1)  # new_paths = [new_path1, new_path2] aka [e1, e2]
                a1 = calculate_loop_activity(interface_graph, new_path1, simulate_reaction=reaction_type, circularize_path=False)
                if a1 == 0:
                    return a1, None
                a1a2.append(a1)
                replacement_loop1['path'] = new_path1
                replacement_loop1['path_spec'] = cmplx.calculate_loop_path_spec(new_path1)
                replacement_loop1['loop_hash'] = cmplx.calculate_loop_hash(replacement_loop1['path_spec'])
                replacement_loop1['a0'] = a0
                replacement_loop1['loop0_hash'] = loop0_hash
                replacement_loop1['old_loop0_id'] = loop0id
                replacement_loop1['activity'] = a1
                replacement_loop1['activity_change_ratio'] = a1/a0
                replacement_loop1['description'] = "loop_formation_effects updated loop (stacked helices).",

                if idx2-idx1 == 1:
                    # Case 1a: Special case: The two joined interface nodes are adjacent; just re-calculate current loop0.
                    # No extra paths
                    new_stacking_loops[loop0_hash] = [replacement_loop1]
                    # Case 0 is shortest_path is a complete subset of loop0, i.e. e3=0.
                    # Case 1 is loop0 is partially overlapping with shortest path,
                    # Case 3 is no overlap with shortest_path; loop0 detected through stacked double-helix.
                    replacement_loop1['description'] = ("Case 0: shared loop affected by duplex stacking "
                                                        "of adjacent ifnodes.")
                    print("WEIRD: ifnode %s and %s are adjacent for STACKING reaction, "
                          "yet they were not processed when considering loops on the shortest path!?")
                    pdb.set_trace()
                else:
                    # Case 1b: Split the existing path into two new loops:
                    # We have already split out one of the loops, just need the other:
                    new_path2 = loop0["path"][idx1:] + loop0["path"][:idx2+1]
                    # new_paths.append(new_path2)
                    a2 = calculate_loop_activity(interface_graph, new_path2, simulate_reaction=reaction_type, circularize_path=False)
                    if a2 == 0:
                        return 0, None
                    a1a2.append(a2)
                    replacement_loop2 = {}
                    replacement_loop2['path'] = new_path2
                    replacement_loop2['path_spec'] = cmplx.calculate_loop_path_spec(new_path2)
                    replacement_loop2['loop_hash'] = cmplx.calculate_loop_hash(replacement_loop2['path_spec'])
                    replacement_loop2['a0'] = a0
                    replacement_loop2['loop0_hash'] = loop0_hash
                    replacement_loop2['old_loop0_id'] = loop0id
                    replacement_loop2['activity'] = a2
                    replacement_loop2['activity_change_ratio'] = a2/a0
                    replacement_loop2['description'] = "Case 0: shared loop splitted by duplex stacking."
                    new_stacking_loops[loop0_hash] = [replacement_loop1, replacement_loop2]
            else:
                # Case 2: Stacking forms a genuine new loop paths: loop0 -> loop1, loop2
                # (I think this will also take care of case 1b and even 1a)
                # We have to generate two new loops, partially composed of the old path
                # and partially of the newly-stacked double-helix:
                # OBS OBS: One of the two loops you are creating will already have been
                # First, split the existing loop0 in two on the nodes that touch the newly-stacked double-helix:
                grouped_path = [(is_in_dh, list(group_iter)) for is_in_dh, group_iter in
                                itertools.groupby(loop0_path, lambda ifnode: ifnode in dh_ifnodes)]
                # if len(grouped_path) in (2, 3): print("Case 1a")
                # should be either 4 or 5 elements.
                # 4 if the loop0 path happens to start/end right between the double-helix and the other parts.
                assert 4 <= len(grouped_path) <= 5
                if len(grouped_path) == 5:
                    # Remove the first element and append it to the last:
                    grouped_path[-1].extend(grouped_path.pop(0))
                    assert grouped_path[0][0][0] == grouped_path[0][-1][0]
                assert len(grouped_path) == 4
                ## TODO: Having done the grouping above, we could probably treat all cases (1a/1b/2a/2b) the same
                if grouped_path[0][0]: # The first group is on the double-helix
                    e1, t1, e2, t2 = [path_group for is_in_dh, path_group in grouped_path]
                else:
                    t1, e2, t2, e1 = [path_group for is_in_dh, path_group in grouped_path]

                # Note: The variable nomenclature here is different from the "loops on shortest path" search above;
                # e1, e2 are parts of the newly-formed arm, i.e. e3a and e3b above.
                new_stacking_loops[loop0id] = []
                for dh_before, t_group, dh_after in ((e1, t1, e2), (e2, t2, e1)):
                    # double-helix before the t-loop, then t-loop, then the double-helix after.
                    # ARGH: One of the loops just formed will already have been formed when processing/"splitting"
                    # the previous loop.
                    # These are signed, + for dh1 and - for dh2:
                    # e1_start, e1_end, e2_start, e2_end = (
                    #     dh_arm_idx_by_ifnode[dh_before[0]], dh_arm_idx_by_ifnode[dh_before[-1]],
                    #     dh_arm_idx_by_ifnode[dh_after[0]], dh_arm_idx_by_ifnode[dh_after[-1]])
                    # dh_start_idx = dh_ifnodes_list.index(dh_before[-1]) # dh_arm_idx_by_ifnode[dh_before[-1]]
                    # dh_end_idx = dh_ifnodes_list.index(dh_after[0]) # dh_arm_idx_by_ifnode[dh_after[0]]
                    if dh_before[-1] in dh1_ifnodes:
                        assert dh_after[0] in dh2_ifnodes
                        arm1, arm2 = dh1_upstream_ifnodes, dh2_upstream_ifnodes
                    else:
                        assert dh_before[-1] in dh2_ifnodes
                        assert dh_after[0] in dh1_ifnodes
                        arm1, arm2 = dh2_upstream_ifnodes, dh1_upstream_ifnodes
                    arm1_idx = arm1.index(dh_before[-1])
                    arm2_idx = arm2.index(dh_after[0])
                    # if arm1_idx == 0 and arm2_idx == 0: print("Case 1(b)")
                    # Revert arm2 so order is outside-in:
                    new_path = arm1[:arm1_idx+1] + t_group + arm2[arm2_idx::-1]  # start-at:stop-before:step-by
                    assert not set(new_path) == new_loop_nodes
                    # This "new_path" is already considered as the shortest path.
                    # Note that this MUSTN'T happen: The present loop0 should have been included in
                    # processed_secondary_loopids when checking shared_loops above.
                    # new_paths.append(new_path)
                    print("Adding new stacking-induced path:\n", new_path)
                    #pdb.set_trace()
                    # ... so, is the idea to just take the product of all of them ?
                    # There is a chance that one of the loops is going to be the primary Λ₁ loop with activity a₁.
                    # So if you do it this way, then make sure you are not including this twice. ₀₁₂
                    # Also, usually I would say a = a₁ * prod(a₂/a₀ for a₀, a₂ in independent loops)
                    # However, here I am not including a₀ in my side-effect factor,
                    # but only including it at the end. Not sure that is right.
                    a2 = calculate_loop_activity(interface_graph, new_path, simulate_reaction=reaction_type, circularize_path=False)
                    if a2 == 0:
                        return 0, None
                    # side_effects_new_loop_activities.append(a2)
                    side_effects_factors.append(a2/a0)
                    a1a2.append(a2)
                    replacement_loop = {}
                    replacement_loop['path'] = new_path
                    replacement_loop['a0'] = a0
                    replacement_loop['old_loop_id'] = loop0id
                    replacement_loop['activity'] = a2
                    replacement_loop['activity_change_ratio'] = a2/a0
                    replacement_loop['description'] = "Case 1+stack: shared loop on shortest_path from stacking."
                    new_stacking_loops[loop0id].append(replacement_loop)
                side_effects_activities.append(a1a2[0]*a1a2[1]/a0)

            # pdb.set_trace()
            # if len(new_paths) == 1: # Optimization for this special case (excessive optimization)
            #     path_activities = [calculate_loop_activity(new_paths[0], simulate_reaction=reaction_type)]
            #     loop_pinching_activity = path_activities[0]/loop0["activity"]
            # else:
            # path_activities = [calculate_loop_activity(new_path, simulate_reaction=reaction_type)
            #                    for new_path in new_paths]
            #loop_pinching_activity = prod(path_activities)/loop0["activity"]
            #side_effects_factors.append(loop_pinching_activity)
            # For the moment, just: old_loopid => dict with new loop paths
            # changed_loops[loop0id] = {
            #     'new_paths': new_paths,
            #     'a1a2': a1a2,   # Primary loop closing activity for each new path
            #     'a0': a0,       # Store this, for good measure...
            #     #'loop_pinching_activity': loop_pinching_activity  # Combined effect of splitting old loop in two
            # }
            processed_secondary_loopids.add(loop0id)
        # end for loop in affected_loops

        # side_effects_factors = [(calculate_loop_activity(loop["path"], simulate_reaction=reaction_type)/
        #                          loop["activity"]) for loop in affected_loops]

        # if any(factor == 0 for factor in side_effects_factors):  # Edit: Is checked as early as possible
        #     return 0
        side_effects_factor = prod(side_effects_factors)


    ## TODO: Compare stacking side_effect_factors and activity_change_ratios
    # if any(factor == 0 for factor in activity_change_ratios):
    #     return 0
    if activity_change_ratios:
        activity_change_ratio = prod(activity_change_ratios)
        print("New loop shortests-path activity: %0.05g" % activity)
        print("Changed loops' activity change ratios and total:",
              " * ".join("%0.03g" % v for v in activity_change_ratios),
              " = %0.03g" % activity_change_ratio)
        print("Activity incl. changed loops' activity change ratios: %0.03g * %0.03g = %0.03g" %
              (activity, activity_change_ratio, activity*activity_change_ratio))
    if side_effects_factors:
        # Stacking side-effect factors include a0 for all loops; a = a1/a0 * a2/a0
        print("Stacking side-effects activities: ",
              " * ".join("%0.02g" % v for v in side_effects_activities),
              " = %0.03g" % prod(side_effects_activities))
        print("Stacking side-effects factors:    ",
              " * ".join("%0.02g" % v for v in side_effects_factors),
              " = %0.03g" % side_effects_factor)
        print("Note: Stacking side-effect factors include a0 for all loops; a = a1/a0 * a2/a0. "
              "I.e. they divide by a0 once too many.")
        #pdb.set_trace()


    """
    ### ALTERNATIVE IDEA: ###
    Assume you are always gonna form exactly one new loop.
    Everything else is side-effects.
    * You form the primary loop (shortest path).
    * Then you look at secondary loops = loops touching the primary loop / shortest path.
    * For secondary loops which were modified:
    * Then you look at tertiary loops = loops touching the modified secondary loops
    * For the tertiary loops that were modified:
    * Repeat, etc.
    At every level you define:
        e1 = the edge which this loop have in common with the previous loop.
        e3 = the edge which split the previous loop in two (e.g. a stacking duplex).
        e2 = the edge in this loop that is not e1

    This should probably be done with a breath-first algorithm using a deque or similar.

    For loops between two stacking double helices, since all existing loops will be split, this should work.
    What if you have more complex networks, where the loops are not just


    ### TODO - a better way to determine the effect of stacking double-helices: ###
    For stacking interactions of domains, there is no guarantee that looking at (primary/shortest loops)
    will detect if the helices being stacked are connected by other means upstream of the stacking ends.
    A better approach would be to detect branch-points:
        1. Find nodes on each arm that are branching off, and order them by distance to the stacking ends,
            with the nodes furthest away first.
        2. For each pair of furthest away from the stacking ends, find the shortest path between the two nodes.
            If the shortest-path goes through a node downstream (closer to the stacking ends),
            make that node the current node and discart the more distant nodes.
        3. Calculate loop activity for path = (arm1 + current_shortest_path + arm2)
        4. Proceed to the next pair in the list until the list is empty.

    """



    # loop_effects = {
    #     'activity': 0,
    #     'shortest_path': path,
    #     'shortest_path_spec': path-of-ifnode-fingerprints
    #     'shortest_path_activity': activity,
    #     'changed_loops': changed_loops,  # old_id => [new_loop1_dict, (new_loop2_dict]
    #     'changed_loops_specs': changed_loops_specs, # loop_state_hash => [path of ifnodes_state_fingerprints]
    #     'loops_considered': processed_secondary_loopids,
    #     'loops_affected': all_affected_loops,
    #     'activity_change_ratios': activity_change_ratios,
    #     'stacking_side_effects_activities': side_effects_activities,
    #     'stacking_side_effects_factors': side_effects_factors,
    # }

    # "activity" variable is loop activity for the newly-created loop.
    # How to get the actual path and side-effects when performing the reaction?
    # Return activity *and* loop info with path and side_effects (dict?) - for caching?
    if activity_change_ratios:
        activity *= activity_change_ratio
    # loop_effects['activity'] = activity # Obsoleted, caused confusion between individual loop's activity and
    # "reaction effects" activity; use total_loop_activity_change instead (or "total_loop_formation_activity")
    loop_effects['total_loop_activity_change'] = (
        loop_effects['new_loop']['activity'] *
        prod(loop['activity_change_ratio'] for loop in loop_effects['changed_loops_by_hash'].values()))
    loop_effects['total_loop_entropy_change'] = ln(loop_effects['total_loop_activity_change'])
    # assert loop_effects['total_loop_activity_change'] == loop_effects['activity']

    return activity, loop_effects


def find_maybe_shorter_e1_e3_subpaths(interface_graph, loop0_nodes, loop1_path, reaction_type):
    r"""
    Arguments:
        loop0_nodes: Nodes in a current loop for which we are considering there might be shorter path.
        loop1_path: A recently created or updated loop forming a new shortest path.
    Note: loop1_path is *NOT* the "possibly shorter candidate path for loop0.
    It is a path that overlaps with loop0 and which is the shortest path for another
    new or recently updated loop, but it is not a path that we are considering for loop0.
    Instead, what we are considering is a new loop path for loop0 formed partially by
    loop0 and partially by loop1_path, specifically the part of loop1 which does not overlap with loop0.
    If you want to compare a current loop with a candidate loop, there is another method for doing that.

    Consider the system below. It currently has a single loop, Λ₀. We have are forming a new
    edge connection, e₃, producing a new loop (Λ₁) with shortest path e₁+e₃.
    We want to know if Λ₀ (e₂+e₁) can be made shorter by using e₃ instead of e₁, forming the updated
    shortest-path e₂+e₃ (Λ₂).
    To do that, we must first determine what e₁ and e₃ sub-paths actually are (and by subtraction e₂).
        Λ₀  .-------,--------.
           /   Λ₁  /  e₃     \  e₂
        e₁ \      /:/    Λ₂  /
            `------´--------´
    This method will find e1 and e3 sub-paths for loop1, where
    e1 is the sub-path that is on Λ₀ and Λ₁ but not Λ₂,
    e3 is the sub-path that is on Λ₁ and Λ₂ but not Λ₀,
    e2 is the sub-path that is on Λ₀ and Λ₂ but not Λ₁.

    Note: The nomenclature of e1, e3 is relative to the loop being considered defined as Λ₀ and the
    recently updated shortest-path for another loop Λ₁.

    That is, after creating the initial "new loop" (= Λ₁) and determining that we should update the old loop,
    Λ₀ --> Λ₂, then we may want to look at a another existing loop, Λ₃ = (e₂+e₄), and determine if this loop
    should use the possibly-shorter path e₃+e₄ instead.
    We would then call find_maybe_shorter_loop2_subpaths(loop0_path=Λ₃=e₂+e₄, loop1_path=Λ₂=e₂+e₃, ...)
    Within this method, Λ₃=e₂+e₄ --> Λ₀ (e₂+e₁)  and  Λ₂=e₂+e₃ --> Λ₁=e₁+e₃.

    Question: Should the two "intersection nodes" where e1, e2 and e3 meet be part of e1, e2, e3?
    It makes sense that they are part of e2 because obviously they are shared.
    However, when comparing the original e1 subpath length with the new e3 subpath, it is nice if both sub-paths
    contain the intersection nodes because we need to include the edge-length to/from the intersection nodes.
    So, although it is not 100% intuitive, we define the intersection nodes
    to be part of e1 and e3 and NOT PART of e2.

    Also note that loop1_path must start and end on e3, i.e. at the new connection being formed.

    """
    ### 1. Find e1 and e3 subpaths: ###
    ### TODO: DOES THIS WORK if loop1_path does not start and end on e3??
    ### This may be the case when working through all the "possibly changed" loops..
    ### Although I think for changed_loops, we always set the path to start at e3.

    e1_not_e3_sub_paths = [(is_shared, list(group_iter)) for is_shared, group_iter in
                          groupby(loop1_path, lambda ifnode: ifnode in loop0_nodes)]
    # grouped_path consists of the first part of e3, then e1, then the last part of e3.
    # Do we ever encounter e2? Well, not when we are grouping shortest path (e3+e1 by definition)..
    # Except if e3a or e3b have zero length (full or half pinching).
    # If e3 == 0 (e3_nodes = Ø), then we'll only have a single group!
    # This can happen both for stacking and hybridization.
    # There are also cases where len(e1_not_e3_sub_paths) == 2!
    assert len(e1_not_e3_sub_paths) <= 3
    if len(e1_not_e3_sub_paths) == 1:
        # If we only have one sub-path, it is e1; e3 is just the two nodes at the start/end of e1.
        assert e1_not_e3_sub_paths[0][0] is True  # verify that e₁ is on Λ₀
        e1 = e1_not_e3_sub_paths[0][1]
        e3a, e3b = [e1[0]], [e1[-1]]
    elif len(e1_not_e3_sub_paths) == 2:
        if e1_not_e3_sub_paths[0][0] is True:
            # We have e3b > 0; e3a starts at e1 (has zero length)
            e1 = e1_not_e3_sub_paths[0][1]
            e3a, e3b = [e1[0]], [e1[-1]]
            e3b.extend(e1_not_e3_sub_paths[1][1])
        else:
            # We have e3a > 0; e3b starts and ends on e1 (has zero length)
            e3a, e1 = e1_not_e3_sub_paths[0][1], e1_not_e3_sub_paths[1][1]
            e3a, e3b = [e1[0]], [e1[-1]]
            e3a.append(e1[0])
    else:
        assert tuple(tup[0] for tup in e1_not_e3_sub_paths) == (False, True, False)
        # Note: grouped_path is grouped by whether node is on Λ₀ or not, not stiffness:
        # The "intersection" nodes are usually considered part of both e1, e2 and e3:
        # intersect_nodes = (e1[0], e1[-1])
        e3a, e1 = e1_not_e3_sub_paths[0][1], e1_not_e3_sub_paths[1][1]
        e3a.append(e1[0])
        # e3b starts with the last node in e1 and then extends the rest of the sub-path:
        e3b = [e1[-1]]
        e3b.extend(e1_not_e3_sub_paths[2][1])

    ### 2. Compare the lengths of e1 and e3 subpaths: ###
    # TODO: Use a Blist to store path-nodes (you will be doing lots of arbitrary inserts/deletions)
    # grouped by stiffness, i.e. similar to segments
    e1_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
                 in group_interfaces_path_by_stiffness(interface_graph, e1)]
    e3a_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
                  in group_interfaces_path_by_stiffness(interface_graph, e3a)] # First part of e3
    e3b_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
                  in group_interfaces_path_by_stiffness(interface_graph, e3b)] # Last part of e3
    ## If e3a and e3b are just single nodes, we won't have any edges!
    ## Maybe it is better to join the groups first before joining the edges?
    ## Or maybe we need to include the intersection node in both e1 and e3a/e3b?

    if e3a_groups and e3b_groups:
        e3_groups = join_two_edge_groups(interface_graph, e3a_groups, e3b_groups, simulate_reaction=reaction_type)
    # elif e3a_groups:
    #     # start of e3a stacks directly with end of e1
    #     e3_groups = join_two_edge_groups(interface_graph, e3a_groups, e1, simulate_reaction=reaction_type)
    else:
        # This is new (but not unexpected; will happen for empty e3a/e3b) - seems to work.
        # TODO: CHECK THIS MORE.
        e3_groups = [group for segment in (e3a_groups, e3b_groups) if segment for group in segment]
    # Uh, it would, perhaps be nice to be able to join/react un-grouped paths,
    # and then group afterwards...
    # edge groups is a list of (stiffness, [(length, len_sq, stiffness, source, target), ...]) tuples.
    # Use edge groups if you need to process something.
    # (it is also slightly cheaper because it does not calculate segment length/length_sq sums).
    if len(e3_groups) <= 1:
        use_e3 = True
    else:
        # Just summing length and see if e3 is shorter than e1.
        e1_length = sum(edge_tup[0] for stiffness, edge_tuples in e1_groups for edge_tup in edge_tuples)
        e3_length = sum(edge_tup[0] for stiffness, edge_tuples in e3_groups for edge_tup in edge_tuples)
        # summing segment-length squared:
        e1_len_sq = sum(sum(etup[0] for etup in edge_tuples)**2 if stiffness > 0 else
                        sum(etup[1] for etup in edge_tuples)
                        for stiffness, edge_tuples in e1_groups)
        e3_len_sq = sum(sum(etup[0] for etup in edge_tuples)**2 if stiffness > 0 else
                        sum(etup[1] for etup in edge_tuples)
                        for stiffness, edge_tuples in e3_groups)
        if e3_length < e1_length and e3_len_sq < e1_len_sq:
            use_e3 = True
        else:
            use_e3 = False
    return e1, e3a, e3b, e3_groups, use_e3






class LoopTracker(object):
    """
    Inheritance vs composition?
    - It is tempting to make this inherit from InterfaceMultiGraph. In that case, it would basically operate on itself.
    On the other hand, if I have other data (e.g. configuration parameters) then the init call signature would change.
    Also: If I want I can move merge/split/delegate/undelegate() methods to this class.
    But it is nice to be able to test these separately and have the functionality available in a separate class.
    Thought materials:
    * https://www.thoughtworks.com/insights/blog/composition-vs-inheritance-how-choose

    System-level vs complex-level?
    (a) Both LoopTracker and InterfaceGraph are system-level singletons.
    (b) System-level InterfaceGraph singleton and per-complex LoopTracker with a reference to the global InterfaceGraph (and a reference to the parent complex).
    (c) Per-complex InterfaceGraph AND LoopTracker instance.
    (d) No class-encapsulation at all, functions all the way (Julia style!)

    Thought food:
    * http://stackoverflow.com/questions/8108688/in-python-when-should-i-use-a-function-instead-of-a-method
    * https://www.thoughtworks.com/insights/blog/composition-vs-inheritance-how-choose
    * http://www.javaworld.com/article/2076814/core-java/inheritance-versus-composition--which-one-should-you-choose-.html

    Questions and criteria:
    * "Do we need to access private attributes/methods?" No -- +1 for not being a Complex member method.
    * "Do we otherwise want to abstract and encapsulate in order to isolate calling code form implementation?"
        - Not really. +1 for not being a member method.
    * "Is the operation being done ON the complex or BY the complex?" Uh...
    * "Do we want to take advantage of polymorphism or inheritance?" No -- +1 for not being a Complex member method.


    """
    def __init__(self, data, **attr):
        self.loop_graph = InterfaceMultiGraph(data, **attr)


    def merge_nodes(self, node1, node2):
        """
        """
        # TODO: Add a deterministic way to determine delegator, delegatee given two arbitrary nodes.
        delegator, delegatee = node1, node2
        return self.loop_graph.merge(delegator, delegatee)

    def split_nodes(self, node1, node2):
        """
        """
        return self.loop_graph.split(node1, node2)




