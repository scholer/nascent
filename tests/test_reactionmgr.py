test_throttles0 = {
    'h+*': {'a': 1.0, 'b': 1.0},
    'h+ ': {'a': 1.0, 'b': 1.0},
    's+ ': {'a': 1.0, 'b': 1.0},
}
# Actual, precise throttles after a run:
test_throttles1 = {
    'h+*': {'a': 0.6895, 'b': 0.076233},
    'h+ ': {'a': 0.00349, 'b': 0.02856},
    # 's+ ': {'a': 0.6, 'b': 0.06, 'A': 0.6, 'B': 0.06},  # Original
    's+ ': {'a': 0.01, 'b': 0.01, 'A': 0.01, 'B': 0.01},
}
# Approximated valus:
test_throttles2 = {
    'h+*': {'a': 0.7, 'b': 0.8},
    'h+ ': {'a': 0.003, 'b': 0.03},
    's+ ': {'a': 0.6, 'b': 0.06},
}
# Throttling intra-complex hyb and stack;
# Making inter-strand reactions un-throttled: -- This shifts equilibrium.
test_throttles3 = {
    'h+*': {'a': 1.0, 'b': 1.0},
    'h+ ': {'a': 0.003, 'b': 0.03},
    's+ ': {'a': 0.6, 'b': 0.06},
}
# Making one side of the graph high, but keeping the other side low:
test_throttles4 = {
    'h+*': {'a': 1.0, 'b': 0.1},
    'h+ ': {'a': 0.003, 'b': 0.1},
    's+ ': {'a': 0.6, 'b': 0.06},
}
# Only throttle stacking: -- Does NOT shift equilibria
test_throttles5 = {
    'h+*': {'a': 1.0, 'b': 1.0},
    'h+ ': {'a': 1.0, 'b': 1.0},
    's+ ': {'a': 0.06, 'b': 0.06},
}
# Only throttle inter-complex hybridization and stacking: - This does NOT shift equilibrium.
# Well, maybe, it does shift *in the other direction* (more dehybridized)
test_throttles6 = {
    'h+*': {'a': 0.1, 'b': 0.1},  # IntER throttled
    'h+ ': {'a': 1.0, 'b': 1.0},  # intra un-throttled
    's+ ': {'a': 0.06, 'b': 0.06},# stack throttled
}
# Only throttle inter-complex hybridization
test_throttles7 = {
    'h+*': {'a': 0.1, 'b': 0.1},  # IntER throttled
    'h+ ': {'a': 1.0, 'b': 1.0},  # intra un-throttled
    's+ ': {'a': 1.0, 'b': 1.0},  # stack un-throttled
}
# Throttling intra-complex hyb and stack;
# Making inter-strand reactions un-throttled: -- This shifts equilibrium.
test_throttles8 = {
    'h+*': {'a': 1.0, 'b': 1.0},
    'h+ ': {'a': 0.03, 'b': 0.03},
    's+ ': {'a': 1.0, 'b': 1.0},
}

