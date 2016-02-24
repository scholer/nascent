

"""



## How fingerprinting is currently done:   (2016/02/22)

1. First a global "state fingerprint" is calculated for the complex as a whole.
    The Complex' state_fingerprint is hash of
        (strands_fingerprint, hybridization_fingerprint, stacking_fingerprint)
    When doing this, we are also counting all domain's "in_complex_identifier" and invoking
    adjust_icid_radius_or_use_instance if any count is larger than 1.
    Q1: How does the domain's in_complex_identifier currently affect the complex state fingerprint?

    strands_fingerprint is currently calculated by the count of each strand.name,
    returning a frozenset of tuples (strand.name: Number of strands of that name)

    hybridization_fingerprint is currently calculated from the frozenset of
"""
{frozenset(((d1.domain_strand_specie, d1.in_complex_identifier),
            (d2.domain_strand_specie, d2.in_complex_identifier))
           for d1, d2 in hybridization_edges())}
"""
    So yes, we are using domain.in_complex_identifier.

    For stacking_fingerprint, we do essentially the same as for hybridization_fingerprint,
    although including the directionality of a stacking interaction (defined to go *from*
    the 5' upstream domain *to* the 3' downstream domain).


2. In order to make domains unique, I'm currently doing this by invoking adjust_icid_radius_or_use_instance():
    This will double Complex.icid_radius up to 5 times,
    each time counting d.in_complex_identifier() for all domain in the complex.

   This *may* give a unique, state-invariant icid labelling of all domains in the complex,
   but that is *not* guaranteed by the current implementation.

   This will *always* fail for symmetric complexes.

   It may also be very expensive to perform these calculations.


"""

"""
        I would venture to say that:
            (a) Complex state fingerprints can be for symmetric complexes. That's fine.
            (b) Domain/DomainEnd nodes within a complex must be made sufficiently unique that
                reaction_spec_pairs do not erroneously match another un-identical pair.
                Consider the following symmetric complex:
                    A -------- B -------- C
                    |                     |
                    C -------- B -------- A
                We must be able to distinguish {A, C} for the case where C is the closest C vs the more distant C.
            TODO: Consider if we have to perform a separate check that DomainEnds are also
                  asymmetric provided an in-complex-identifier for the parent domain.
        How about:
            (a) Complex state is initially based on the connections from each Domain/DomainEnd.
                Q1: If it is just a matter of comparing edges, why is graph isomorphism such a hard problem?
                A1: Because the graph nodes are not labelled, so you have to perform permutations to see:
                    "If we exchange node 1 and 2 in the adjacency list, are the graph adjacency lists then identical?"
                Q2: Do we need nodes to be uniquely labelled before we can make a proper fingerprint based
                    only on the adjacency list?
                Q2a: Can we find a way to have two different complexes (in different states) give the same
                      fingerprint


"""

"""

## Propoal for a more efficient node-labelling scheme for symmetric complexes: ##  (2016/02/22)

Nodes can be either Domains or DomainEnds.

Make a set of all nodes:
non_unique_nodes = all_nodes
For each node, clear node.icid by setting it to node.name.
Loop while non_unique_nodes and ring_number < max_icid_radius:
    Remove all nodes with unique icid.

    # For each node make a "ring0" hash based on the node's connection with other nodes (by their type/name not instance).
    # Then, for each node make a ring1 hash as a frozenset of (interaction, neighbor_ring0_hash) for all neighbors.

See Complex.label_complex_nodes() method.

"""
