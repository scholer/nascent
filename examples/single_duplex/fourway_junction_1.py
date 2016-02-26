import os, yaml
desc = """
####          Fourway Junction 1        ####
#         A       .        B
# 5' ------------' `--------------- 3'
# 3' -----a------. .-------b------- 5'
#                | |                  .
# 5' -----d------' `-------c------- 3'
# 3' ------------...--------------- 5'
#         D                C
# See also: fourway1
# Domains:
# - A, a	CGACTGGAAAGCGGGC, GCCCGCTTTCCAGTCG
# - T, TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
# - B, b	TGACCATGATTACGAA, TTCGTAATCATGGTCA
# - C, c	TGGCGCCCAATACGCA, TGCGTATTGGGCGCCA
# - D	TGATTTATAAGGGATTTTGCCGATTTCGGAAC
# - E, h2b	GGCACGACAGGTTTCC, GGAAACCTGTCGTGCC
#### Fourway Junction 1 conformation AB/CD     ####
#              A         .        B
#  A 5' CGACTGGAAAGCGGGC' `TGACCATGATTACGAA 3'  B
#  a 3' GCTGACCTTTCGCCCG. .ACTGGTACTAATGCTT 5'  b
#                       | |
#  d 5' GGAAACCTGTCGTGCC' `TGCGTATTGGGCGCCA 3' c
#  D 3' CCTTTGGACAGCACGG...ACGCATAACCCGCGGT 5' C
#              D                  C
#### Fourway Junction 1 conformation AD/CB     ####
#              b         .        c
#  b 5' TTCGTAATCATGGTCA' `TGCGTATTGGGCGCCA 3'  c
#  B 3' AAGCATTAGTACCAGT. .ACGCATAACCCGCGGT 5'  C
#                       | |
#  A 5' CGACTGGAAAGCGGGC' `GGCACGACAGGTTTCC 3' D
#  a 3' GCTGACCTTTCGCCCG...CCGTGCTGTCCAAAGG 5' d
#              a                  d
#
# Most prevalent is AD/CB, but CB not stacked.

"""
seqs = """
    A: CGACTGGAAAGCGGGC
    a: GCCCGCTTTCCAGTCG
    B: TGACCATGATTACGAA
    b: TTCGTAATCATGGTCA
    C: TGGCGCCCAATACGCA
    c: TGCGTATTGGGCGCCA
    D: GGCACGACAGGTTTCC
    d: GGAAACCTGTCGTGCC
"""
strands = """
    s1: [A, B]
    s2: [b, c]
    s3: [C, D]
    s4: [d, a]
"""
seqs = yaml.load(seqs)
strands = yaml.load(strands)
headers = "Strand	Domains	Sequence".split()
sep = "\t"
WC = dict(zip("ATGC", "TACG"))
rcompl = lambda seq: "".join(WC[b] for b in seq[::-1])
for name, seq in seqs.items():
    cname = name.lower() if name == name.upper() else name.upper()
    if cname in seqs:
        assert seq == rcompl(seqs[cname])
fn = os.path.splitext(__file__)[0]+".txt"
print("Writing output to file:", fn)
with open(fn, 'w') as fp:
    output = "\n".join((
        sep.join(headers),
        "\n".join(sep.join((k, ", ".join(dnames), ", ".join(seqs[domain] for domain in dnames)))
                  for k, dnames in strands.items()),
        desc))
    print(output)
    fp.write(output)
