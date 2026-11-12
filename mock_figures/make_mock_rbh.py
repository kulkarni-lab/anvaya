#!/usr/bin/env python3
"""Builds mock .rbh datasets with a designed, known structure mimicking the
real analysis: parasitiform ingroup, acariform outgroup, two arachnopulmonates
carrying a shared WGD. Used only to render figure templates."""

import os

REF = "cALGconsensus"

# 12 linkage groups with a realistic, decaying size distribution
LG_SIZES = {"cALG1": 26, "cALG2": 22, "cALG3": 19, "cALG4": 17, "cALG5": 15, "cALG6": 13,
            "cALG7": 11, "cALG8": 10, "cALG9": 8, "cALG10": 7, "cALG11": 6, "cALG12": 5}
LGS = list(LG_SIZES)

# True tree: (Sarc,(Ixod,(Varr,Neoseiu)))
#   cALG1+cALG2  fused clean  -> Ixod, Varr, Neoseiu   (parasitiform synapomorphy)
#   cALG3+cALG4  fused clean  -> Varr, Neoseiu         (mesostigmatid synapomorphy)
#   cALG5+cALG6  fused MIXED  -> Neoseiu               (shows the mixed state)
#   cALG7+cALG8  fused clean  -> Ixod and Neoseiu      (deliberate homoplasy)
FUSIONS = {
    "Sarc":    [],
    "Ixod":    [(("cALG1", "cALG2"), "clean"), (("cALG7", "cALG8"), "clean")],
    "Varr":    [(("cALG1", "cALG2"), "clean"), (("cALG3", "cALG4"), "clean")],
    "Neoseiu": [(("cALG1", "cALG2"), "clean"), (("cALG3", "cALG4"), "clean"),
                (("cALG5", "cALG6"), "mixed"), (("cALG7", "cALG8"), "clean")],
    "Tvan":    [],
    "Atypus":  [],
}
WGD_TAXA = {"Tvan", "Atypus"}          # arachnopulmonates: every LG retained twice

COLS = ["rbh", "gene_group", "color", f"{REF}_gene", f"{REF}_scaf", f"{REF}_pos",
        "{sp}_gene", "{sp}_scaf", "{sp}_pos", "{sp}_breakchrom", "{sp}_ix", "{sp}_break_ix",
        f"{REF}_ix", f"{REF}_break_ix", f"{REF}_breakchrom", "whole_FET", "break_FET",
        f"{REF}_D", "{sp}_D", "{sp}_plotindex", "{sp}_plotpos",
        f"{REF}_plotindex", f"{REF}_plotpos", "orthology_type"]


def layout(sp, wgd=False, rescue=True):
    """-> list of (alg, ref_pos, scaffold, pos, orthology_type)"""
    fus = {}
    for (a, b), kind in FUSIONS.get(sp, []):
        fus[a] = (f"{sp}_{a}_{b}", kind, a, b)
        fus[b] = (f"{sp}_{a}_{b}", kind, a, b)

    rows, pos = [], {}
    copies = ["cpA", "cpB"] if wgd else [None]
    for copy in copies:
        for lg in LGS:
            n = LG_SIZES[lg]
            if lg in fus:
                scaf, kind, a, b = fus[lg]
            else:
                scaf, kind = f"{sp}_{lg}", None
            if copy:
                scaf = f"{scaf}_{copy}"
            for i in range(n):
                # mixed fusion: interleave the two groups by alternating offsets
                if kind == "mixed":
                    order_key = i * 2 + (0 if lg == fus[lg][2] else 1)
                else:
                    order_key = (0 if (kind and lg == fus[lg][2]) else n) + i
                pos.setdefault(scaf, {})
                pos[scaf][order_key] = pos[scaf].get(order_key, 0)
                otype = "direct"
                if wgd and copy == "cpB":
                    # under plain RBH the second copy is lost; rescue brings it back
                    otype = "wgd_inferred" if rescue else None
                if otype is None:
                    continue
                rows.append((lg, i, scaf, order_key, otype))
    # assign monotonically increasing coordinates within each scaffold
    out, counters = [], {}
    for lg, ref_pos, scaf, order_key, otype in sorted(rows, key=lambda r: (r[2], r[3])):
        counters[scaf] = counters.get(scaf, 0) + 10000
        out.append((lg, ref_pos, scaf, counters[scaf], otype))
    return out


def write(sp, rows, outdir):
    header = [c.format(sp=sp) for c in COLS]
    with open(os.path.join(outdir, f"{sp}.rbh"), "w") as f:
        f.write("\t".join(header) + "\n")
        for i, (lg, ref_pos, scaf, p, otype) in enumerate(rows):
            gp = LGS.index(lg) * 1000 + ref_pos
            rec = {"rbh": f"rbh_{i}", "gene_group": lg, "color": "",
                   f"{REF}_gene": f"{lg}_{ref_pos}", f"{REF}_scaf": lg, f"{REF}_pos": gp,
                   f"{sp}_gene": f"{sp}_{lg}_{ref_pos}_{scaf}", f"{sp}_scaf": scaf, f"{sp}_pos": p,
                   f"{sp}_breakchrom": f"{scaf}:0-end", f"{sp}_ix": i, f"{sp}_break_ix": i,
                   f"{REF}_ix": i, f"{REF}_break_ix": i, f"{REF}_breakchrom": f"{lg}:0-end",
                   "whole_FET": 1e-14, "break_FET": 1e-14, f"{REF}_D": "", f"{sp}_D": "",
                   f"{sp}_plotindex": i, f"{sp}_plotpos": p,
                   f"{REF}_plotindex": i, f"{REF}_plotpos": gp,
                   "orthology_type": otype}
            f.write("\t".join(str(rec[c.format(sp=sp)]) for c in COLS) + "\n")


def build(outdir, rescue=True, include=None):
    os.makedirs(outdir, exist_ok=True)
    taxa = include or list(FUSIONS)
    for sp in taxa:
        write(sp, layout(sp, wgd=sp in WGD_TAXA, rescue=rescue), outdir)
    return taxa


if __name__ == "__main__":
    t = build("mockfigs/rbh_main", rescue=True)
    print("main dataset:", t)
    build("mockfigs/rbh_norescue", rescue=False, include=["Sarc", "Ixod", "Tvan", "Atypus"])
    build("mockfigs/rbh_rescue", rescue=True, include=["Sarc", "Ixod", "Tvan", "Atypus"])
    print("WGD before/after datasets written")
