import os, yaml
desc = """
####          Circfb design (Circular FoldBack):        ####
#                            E          A          T
# 5' .---------------------------3'5'---------------------.
#    |                 .----------------------.            |
#    |                (                        )           |
#    |                 `---------5' 3'--------'            |
#    |_____________________________________________________|
#            D               C          B
# Domains:
# - A, a	CGACTGGAAAGCGGGC, GCCCGCTTTCCAGTCG
# - T, TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
# - B, b	TGACCATGATTACGAA, TTCGTAATCATGGTCA
# - C, c	TGGCGCCCAATACGCA, TGCGTATTGGGCGCCA
# - D	TGATTTATAAGGGATTTTGCCGATTTCGGAAC
# - E, h2b	GGCACGACAGGTTTCC, GGAAACCTGTCGTGCC
"""
seqs = """
    A: CGACTGGAAAGCGGGC
    a: GCCCGCTTTCCAGTCG
    T: TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    B: TGACCATGATTACGAA
    b: TTCGTAATCATGGTCA
    C: TGGCGCCCAATACGCA
    c: TGCGTATTGGGCGCCA
    D: TGATTTATAAGGGATTTTGCCGATTTCGGAAC
    E: GGCACGACAGGTTTCC
    e: GGAAACCTGTCGTGCC
"""
strands = """
    templ: [A, T, B, C, D, E]
    splint: [c, e, a, b]
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
