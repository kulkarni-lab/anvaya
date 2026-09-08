#!/usr/bin/env python3
"""
genomes_to_rbh.py
==================

Bridges raw genome assemblies - the exact CDS/GFF3/[pep] input layout used
by run_wgd_pipeline.py - to the .rbh input format alg_parsimony.py needs,
then runs alg_parsimony.py on the result.

--------------------------------------------------------------------------
WHY THIS EXISTS / WHAT IT ASSUMES (READ THIS FIRST)
--------------------------------------------------------------------------
alg_parsimony.py needs, for every taxon, which orthology group ("linkage
group") each gene belongs to, plus where that gene sits in the taxon's own
genome. In the original odp-style workflow that group comes from an
external, independently-reconstructed reference (e.g. Simakov et al. 2022's
BCnS ancestral linkage groups, or an odp-computed "odp_LGs_from_<taxa>" set).

*** This script assumes you do NOT have any such reference. *** Instead it
manufactures ALG-proxy groups from your OWN input genomes, in one of two
modes:

  --reference-species NAME
      Use one genome's own chromosomes as the "linkage groups" - fast,
      simple, but only as good as that one genome's assembly quality and
      how ancestral its karyotype actually is.

  --consensus-species SP1 SP2 [SP3 ...] | --consensus-species all
      Reconstruct a de novo consensus linkage-group set from SEVERAL
      genomes at once, the way odp itself does it (in simplified form):
      compute reciprocal best hits among every pair of the given "core"
      species, cluster genes into orthogroups (connected components of the
      RBH graph), then group orthogroups into linkage groups wherever they
      stay co-located on the same chromosome across (by default) every
      core genome. This is the more defensible option when you have 3+
      reasonably-conserved genomes to spare for building the reference -
      it isn't tied to any single genome's idiosyncrasies. EVERY input
      genome (including the ones used to build the consensus set) gets its
      own .rbh file and becomes a taxon in the downstream tree - unlike
      --reference-species mode, where the reference genome itself is
      necessarily excluded (there's nothing to compare it to).

Both modes are real, defensible ways to ask "how has genome architecture
diverged from a chosen reference" - but neither is the same claim as
testing against deep, independently-validated ancestral linkage groups:
  * --reference-species: fusion/mixing calls are only as good as the
    reference genome you pick. A poorly-assembled or already
    highly-rearranged reference will manufacture spurious "fusions" (and
    miss real ones) everywhere. Pick the most contiguous, least-rearranged
    genome you have - often the outgroup taxon.
  * --consensus-species: only as good as the "core" genomes you feed it,
    and requiring co-location in every core genome (the default) throws
    out any orthogroup missing from even one of them - more core species
    means better-supported linkage groups but a stricter, smaller set (and
    more diamond runs: pairwise among the core, then core-vs-everyone-else).
  * Either way, say explicitly in any write-up that these are "fusions
    relative to <reference/consensus set>," not ancestral linkage-group
    fusions with independent validation behind them.
  * In --reference-species mode, the reference species itself is never
    written as one of the output .rbh files (nothing to compare it to) and
    so never appears as a taxon in the downstream tree - exactly as
    BCnSSimakov2022 itself never appeared as a taxon in the original
    odp-style workflow. In --consensus-species mode every input genome,
    core or not, DOES get an .rbh file and a place in the tree.

--------------------------------------------------------------------------
OPTIONAL: WHOLE-GENOME/SEGMENTAL DUPLICATION (--detect-wgd)
--------------------------------------------------------------------------
Reciprocal best hit is, by construction, a 1:1 mapping. If the species being
compared to the reference/core set carries a lineage-specific whole-genome
or large segmental duplication that the reference does not share (the
textbook case: comparing a spider, descended from the ancestral
arachnopulmonate WGD, to a non-duplicated outgroup such as a tick), every
ancestral gene generally has TWO surviving copies (ohnologs) in that
species - but plain RBH can only ever recover one of them. The other is
silently absent from the orthology table, which can make a real ALG
disappear from the significance calling entirely (diluted below the
--min-genes/--alpha threshold on BOTH ohnologous chromosomes) rather than
show up as an obvious doubled hit - which is exactly what makes this
failure mode easy to miss: an odp-style ribbon plot between a non-
duplicated and a duplicated genome will generally NOT show a clean "1:2"
pattern, for this same reason.

--detect-wgd fixes this at the source: for each requested species (default:
every species compared to the reference/core set; restrict with
--wgd-species), an intragenomic self-comparison finds reciprocal-best-hit
PARALOG pairs, aggregates them per scaffold-pair, and confirms genuine
ohnologous scaffold pairs (as opposed to incidental one-off paralogy) via a
minimum-pair-count threshold (--wgd-min-paralog-pairs) and, where there is
other paralogy elsewhere in the genome to compare against, a Fisher's exact
enrichment test (--wgd-alpha). Once a scaffold pair is confirmed
ohnologous, any accepted RBH pair whose species-side gene has a paralog on
the confirmed partner scaffold has that paralog added back in as a second,
explicitly flagged ortholog of the same reference gene - recovering exactly
the copy plain RBH structurally drops. Detected scaffold pairs are written
to duplication_blocks.tsv, and every .rbh file gains an orthology_type
column ('direct' vs 'wgd_inferred') so the rescued rows stay identifiable.

Scope: this only rescues ohnologs on the "being compared" (non-reference,
non-core) side of each search - it does not detect or correct for
duplication private to the reference genome or the consensus core species
themselves, since orthogroup construction there assumes the core set is not
independently duplicated.

--------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------
Same directory layout as run_wgd_pipeline.py:
    input_dir/
        A.cds     (required)  nucleotide CDS FASTA
        A.gff3    (required, or A.gtf)
        A.pep     (optional)  protein FASTA; translated from .cds if absent
        B.cds
        B.gff3
        ...
Species are auto-discovered from every "*.cds" file, exactly like
run_wgd_pipeline.py's discover_species().

If you already ran run_wgd_pipeline.py on this same --input-dir, pass its
--output-dir here as --wgd-output-dir and this script will reuse the
protein FASTAs it already translated (skipping that step). Nothing else
from that run is reusable - its diamond step is a self-blast within each
genome, not cross-species orthology.

--------------------------------------------------------------------------
WHAT IT DOES
--------------------------------------------------------------------------
--reference-species R mode, for every other discovered species S:
  1. Translate CDS -> protein for R and S if no .pep was supplied (or found
     under --wgd-output-dir).
  2. diamond makedb + blastp R-vs-S and S-vs-R (two one-way searches - a
     genuine two-genome reciprocal-best-hit, not a self-blast).
  3. Take the reciprocal best hits and map each gene to its own genome's
     scaffold/chromosome + coordinate via its GFF3/GTF (reusing the same
     permissive multi-key ID-matching logic as run_wgd_pipeline.py, so
     NCBI/AUGUSTUS/Ensembl-style annotations are all handled the same way).
  4. For each (S-scaffold, R-chromosome) pair, run a Fisher's exact test
     for enrichment (mirrors odp's whole_FET column) so alg_parsimony.py's
     significance filtering still means something; falls back to leaving
     whole_FET blank if scipy isn't installed (--min-genes-only filtering
     in alg_parsimony.py still works either way).
  5. Write S.rbh in the column layout alg_parsimony.py expects, with
     gene_group = the reference chromosome/scaffold each ortholog sits on.

--consensus-species mode, given a "core" set of species:
  1. Translate CDS -> protein for every discovered species.
  2. diamond blastp (both directions) between every PAIR of core species,
     and take reciprocal best hits from each pair.
  3. Cluster genes into orthogroups = connected components of the RBH
     graph across all core-vs-core comparisons. Any component containing
     more than one gene from the same core species is dropped as
     ambiguous (probable leaked paralogy).
  4. Keep orthogroups present in at least --consensus-min-species core
     species (default: all of them).
  5. For each pair of surviving orthogroups, check whether they sit on the
     same chromosome in the core genomes that have both: if that holds in
     at least --consensus-min-frac of those genomes (default 1.0 = every
     one of them), link the two orthogroups. Connected components of that
     linkage graph are the reconstructed linkage groups ("cALG1", "cALG2",
     ... - largest first).
  6. Every core species gets an .rbh file directly from its own orthogroup
     memberships and coordinates. Every non-core species gets RBH'd
     against each core species in turn; a gene's linkage-group call is
     whichever orthogroup/ALG the majority of those hits agree on.
  7. whole_FET and file-writing are the same as --reference-species mode,
     with gene_group = the reconstructed cALG rather than a real
     chromosome name, and the "reference" position axis anchored to the
     first-listed core species (for dot-plot x-axis comparability).

Finally (unless --rbh-only), invokes alg_parsimony.py on the resulting .rbh
directory. Anything you put after --forward is passed straight through to
alg_parsimony.py, e.g. --forward --min-genes 3 --constraint "(...)".

--------------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------------
On PATH: diamond. Python: biopython (only needed if some .pep is missing
and must be translated), scipy (only for the whole_FET enrichment column -
skipped gracefully if absent), pandas, and alg_parsimony.py itself
(auto-detected next to this script, in the current directory, or on PATH;
override with --parsimony-script).

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------
    python3 genomes_to_rbh.py --input-dir ./species_data \\
        --reference-species Ixod \\
        --rbh-outdir ./rbh_out --parsimony-outdir ./parsimony_out \\
        --forward --min-genes 3 --alpha 0.05

    # reuse translations from a prior run_wgd_pipeline.py run
    python3 genomes_to_rbh.py --input-dir ./species_data --wgd-output-dir ./wgd_run \\
        --reference-species Ixod --rbh-outdir ./rbh_out --parsimony-outdir ./parsimony_out

    # only build the .rbh files, don't run the parsimony search
    python3 genomes_to_rbh.py --input-dir ./species_data --reference-species Ixod \\
        --rbh-outdir ./rbh_out --rbh-only

    # multi-genome consensus linkage groups instead of a single reference
    python3 genomes_to_rbh.py --input-dir ./species_data \\
        --consensus-species Avul Ixod Nlon Varr \\
        --rbh-outdir ./rbh_out --parsimony-outdir ./parsimony_out \\
        --forward --min-genes 3 --plot-dotplots --plot-trees

    # a species known/suspected to carry a lineage-specific whole-genome
    # duplication (e.g. a spider being compared to a non-duplicated tick):
    # rescue the ohnolog reciprocal-best-hit structurally drops
    python3 genomes_to_rbh.py --input-dir ./species_data \\
        --reference-species Ixod --detect-wgd --wgd-species Pbig \\
        --rbh-outdir ./rbh_out --parsimony-outdir ./parsimony_out
"""

import argparse
import re
import shutil
import subprocess
import sys
import zlib
from collections import defaultdict, Counter
from itertools import combinations
from pathlib import Path

import pandas as pd

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    SeqIO = None
    Seq = None

try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None


DEFAULT_PALETTE = [
    "#AB7E26", "#47957F", "#7B7B7B", "#3B69B2", "#764A9B", "#C76E00",
    "#8B0000", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
]

GFF_EXTENSIONS = (".gff3", ".gff")
GTF_EXTENSIONS = (".gtf",)


def log(msg):
    print(f"[genomes-to-rbh] {msg}", flush=True)


def run_cmd(cmd, **kwargs):
    log("RUN: " + " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kwargs)


def check_tool(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required external tool '{name}' not found on PATH.")


# ==========================================================================
# Species discovery + GFF/GTF gene->chromosome mapping
# (same conventions/robustness as run_wgd_pipeline.py, reused verbatim)
# ==========================================================================

def discover_species(input_dir: Path):
    species = {}
    for cds_file in sorted(input_dir.glob("*.cds")):
        name = cds_file.stem
        gff_file, gff_is_gtf = None, False
        for ext in GFF_EXTENSIONS:
            candidate = input_dir / f"{name}{ext}"
            if candidate.exists():
                gff_file, gff_is_gtf = candidate, False
                break
        if gff_file is None:
            for ext in GTF_EXTENSIONS:
                candidate = input_dir / f"{name}{ext}"
                if candidate.exists():
                    gff_file, gff_is_gtf = candidate, True
                    break
        if gff_file is None:
            log(f"WARNING: skipping '{name}' - no matching .gff3/.gff/.gtf found next to {cds_file.name}")
            continue
        pep_file = None
        for ext in (".pep", ".faa", ".fasta.pep"):
            candidate = input_dir / f"{name}{ext}"
            if candidate.exists():
                pep_file = candidate
                break
        species[name] = {"cds": cds_file, "gff": gff_file, "gff_is_gtf": gff_is_gtf, "pep": pep_file}
    return species


_GFF3_ID_KEYS = ("ID", "Name", "Parent", "protein_id", "transcript_id", "gene_id")
_GTF_ID_KEYS = ("ID", "Name", "gene_id", "transcript_id", "protein_id", "gene_name", "transcript_name")
_GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_gff3_line_ids(attr_field):
    ids = []
    for kv in attr_field.split(";"):
        kv = kv.strip()
        if "=" not in kv:
            continue
        key, val = kv.split("=", 1)
        if key in _GFF3_ID_KEYS:
            for v in val.split(","):
                v = v.strip()
                if v:
                    ids.append(v)
    return ids


def parse_gtf_line_ids(attr_field):
    ids = []
    for key, val in _GTF_ATTR_RE.findall(attr_field):
        if key in _GTF_ID_KEYS and val:
            ids.append(val)
    return ids


def build_gene_to_chr(gff_file, is_gtf):
    id_to_chr = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= 8:
                continue
            chrom = parts[0]
            ids = parse_gtf_line_ids(parts[8]) if is_gtf else parse_gff3_line_ids(parts[8])
            for gid in ids:
                id_to_chr.setdefault(gid, chrom)
    log(f"  -> extracted {len(id_to_chr)} ID-to-chromosome mappings from {gff_file.name}")
    return id_to_chr


_ID_SUFFIXES_TO_STRIP = (".cds", ".CDS", ".protein", ".pep", ".p1", "-RA", "-mRNA-1")


def lookup_chr(gene_id, gene_to_chr):
    if gene_id in gene_to_chr:
        return gene_to_chr[gene_id]
    for suf in _ID_SUFFIXES_TO_STRIP:
        if gene_id.endswith(suf) and gene_id[: -len(suf)] in gene_to_chr:
            return gene_to_chr[gene_id[: -len(suf)]]
    m = re.match(r"^(.*)\.\d+$", gene_id)
    if m and m.group(1) in gene_to_chr:
        return gene_to_chr[m.group(1)]
    return None


def lookup_pos(gene_id, gene_to_pos):
    if gene_id in gene_to_pos:
        return gene_to_pos[gene_id]
    for suf in _ID_SUFFIXES_TO_STRIP:
        if gene_id.endswith(suf) and gene_id[: -len(suf)] in gene_to_pos:
            return gene_to_pos[gene_id[: -len(suf)]]
    m = re.match(r"^(.*)\.\d+$", gene_id)
    if m and m.group(1) in gene_to_pos:
        return gene_to_pos[m.group(1)]
    return None


def build_gene_to_pos(gff_file, is_gtf):
    """Like build_gene_to_chr, but returns the feature's start coordinate."""
    id_to_pos = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= 8:
                continue
            try:
                start = int(parts[3])
            except ValueError:
                continue
            ids = parse_gtf_line_ids(parts[8]) if is_gtf else parse_gff3_line_ids(parts[8])
            for gid in ids:
                id_to_pos.setdefault(gid, start)
    return id_to_pos


# ==========================================================================
# Translation (reused from run_wgd_pipeline.py)
# ==========================================================================

def translate_cds(cds_file, protein_file, genetic_code=1):
    if SeqIO is None:
        raise RuntimeError("biopython is required for translation (pip install biopython)")
    records = []
    for record in SeqIO.parse(str(cds_file), "fasta"):
        seq = str(record.seq).upper().replace(" ", "").replace("\n", "")
        protein = Seq(seq).translate(table=genetic_code, to_stop=False)
        if protein.endswith("*"):
            protein = protein[:-1]
        record.seq = protein
        record.description = f"translated_from={record.id}"
        records.append(record)
    SeqIO.write(records, str(protein_file), "fasta")
    log(f"  translated {len(records)} CDS sequences -> {protein_file}")


def ensure_protein_fasta(name, info, work_dir, wgd_output_dir, genetic_code):
    if info["pep"] is not None:
        return info["pep"]
    if wgd_output_dir is not None:
        candidate = wgd_output_dir / name / f"{name}.faa"
        if candidate.exists():
            log(f"  reusing translated protein FASTA from wgd run: {candidate}")
            return candidate
    faa_file = work_dir / f"{name}.faa"
    if not faa_file.exists():
        log(f"  translating CDS -> protein for {name}")
        translate_cds(info["cds"], faa_file, genetic_code=genetic_code)
    return faa_file


# ==========================================================================
# Minimum scaffold size filter (--min-scaffold-length)
# ==========================================================================

def compute_scaffold_lengths(gff_file, is_gtf):
    """Scaffold/chromosome lengths in bp. Prefers GFF3 '##sequence-region
    seqid start end' pragma lines (the authoritative source, when present);
    falls back to the maximum feature end coordinate observed on that
    scaffold for any scaffold without one (GTF never has these pragmas, and
    not every GFF3 does either) - an underestimate, but the only length
    information available without a separate genome assembly FASTA."""
    lengths = {}
    max_end_seen = defaultdict(int)
    with open(gff_file) as f:
        for line in f:
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        lengths[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= 4:
                continue
            try:
                end = int(parts[4])
            except ValueError:
                continue
            if end > max_end_seen[parts[0]]:
                max_end_seen[parts[0]] = end
    approximated = 0
    for seqid, end in max_end_seen.items():
        if seqid not in lengths:
            lengths[seqid] = end
            approximated += 1
    if approximated:
        log(f"  -> {approximated} scaffold(s) had no '##sequence-region' pragma in the GFF3/GTF; "
            f"using the longest annotated feature's end coordinate as an approximate (minimum) length")
    return lengths


def filter_fasta_by_scaffold(faa_in, faa_out, gene_to_chr, scaffold_lengths, min_length):
    """Write a copy of faa_in to faa_out, dropping sequences whose gene sits
    on a scaffold shorter than min_length (or of unknown length). Filtering
    happens here, before diamond ever sees these sequences, rather than
    after the fact, so a small-scaffold gene can never win - and so cost a
    real ortholog - a reciprocal-best-hit slot it would otherwise be
    discarded from anyway. Returns (n_kept, n_excluded, excluded_scaffolds)."""
    n_kept, n_excluded = 0, 0
    excluded_scaffolds = set()
    write_this = False
    with open(faa_in) as fin, open(faa_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                gid = line[1:].strip().split()[0]
                chrom = lookup_chr(gid, gene_to_chr)
                length = scaffold_lengths.get(chrom) if chrom else None
                write_this = length is not None and length >= min_length
                if write_this:
                    n_kept += 1
                else:
                    n_excluded += 1
                    if chrom:
                        excluded_scaffolds.add(chrom)
            if write_this:
                fout.write(line)
    return n_kept, n_excluded, excluded_scaffolds


def get_search_ready_faa(name, info, work_dir, wgd_output_dir, genetic_code, gene_to_chr, args):
    """ensure_protein_fasta(), then apply --min-scaffold-length filtering if requested."""
    faa_path = ensure_protein_fasta(name, info, work_dir, wgd_output_dir, genetic_code)
    if not args.min_scaffold_length or args.min_scaffold_length <= 0:
        return faa_path
    scaffold_lengths = compute_scaffold_lengths(info["gff"], info["gff_is_gtf"])
    filtered_path = work_dir / f"{name}.minscaf.faa"
    n_kept, n_excluded, excluded_scaffolds = filter_fasta_by_scaffold(
        faa_path, filtered_path, gene_to_chr, scaffold_lengths, args.min_scaffold_length)
    if n_excluded:
        log(f"  -> --min-scaffold-length {args.min_scaffold_length}: excluded {n_excluded} gene(s) on "
            f"{len(excluded_scaffolds)} scaffold(s) < threshold for {name} ({n_kept} kept)")
    if n_kept == 0:
        longest = max(scaffold_lengths.values()) if scaffold_lengths else 0
        sys.exit(f"--min-scaffold-length {args.min_scaffold_length} excludes every scaffold in {name} "
                  f"(longest computed scaffold length: {longest} bp) - nothing is left to search. Lower "
                  f"--min-scaffold-length, or omit it for this species, or (if scaffold lengths are being "
                  f"approximated from annotated features rather than read from a '##sequence-region' "
                  f"pragma - check the log above) supply real scaffold lengths by adding "
                  f"'##sequence-region' pragma lines to {name}'s GFF3.")
    return filtered_path


# ==========================================================================
# Cross-species diamond RBH
# ==========================================================================

_BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def run_diamond_oneway(query_faa, db_faa, work_dir, tag, evalue, max_target_seqs, threads,
                        reuse_blast=True):
    db_prefix = work_dir / f"db_{db_faa.stem}"
    if not Path(f"{db_prefix}.dmnd").exists():
        run_cmd(["diamond", "makedb", "--in", str(db_faa), "-d", str(db_prefix), "--quiet"])
    # The search parameters are baked into the cache filename, so runs that differ in
    # --evalue or --max-target-seqs never collide on the same cached file. This makes
    # repeated invocations over the same genomes (e.g. a consensus-threshold parameter
    # sweep, where only post-alignment options change) skip the expensive alignment
    # step entirely rather than recomputing an identical result each time.
    sig = f"e{evalue}_m{max_target_seqs}"
    out_tsv = work_dir / f"{tag}.{sig}.blast.tsv"
    if reuse_blast and out_tsv.exists() and out_tsv.stat().st_size > 0:
        log(f"  reusing cached alignment: {out_tsv.name}")
        return out_tsv
    run_cmd([
        "diamond", "blastp",
        "--db", f"{db_prefix}.dmnd", "--query", str(query_faa),
        "--outfmt", "6", *_BLAST_COLS,
        "--evalue", str(evalue), "--max-target-seqs", str(max_target_seqs),
        "--threads", str(threads), "--out", str(out_tsv), "--quiet",
    ])
    return out_tsv


def best_hits(blast_tsv):
    """query -> best (highest-bitscore) subject."""
    df = pd.read_csv(blast_tsv, sep="\t", names=_BLAST_COLS)
    if df.empty:
        return {}
    best = df.sort_values(["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid")
    return dict(zip(best["qseqid"], best["sseqid"]))


def compute_rbh(fwd_blast_tsv, rev_blast_tsv):
    """fwd: reference(query) -> species(subject). rev: species(query) -> reference(subject).
    Returns list of (ref_gene, sp_gene) reciprocal best hit pairs."""
    fwd = best_hits(fwd_blast_tsv)   # ref_gene -> sp_gene
    rev = best_hits(rev_blast_tsv)   # sp_gene -> ref_gene
    pairs = []
    for ref_gene, sp_gene in fwd.items():
        if rev.get(sp_gene) == ref_gene:
            pairs.append((ref_gene, sp_gene))
    return pairs


def compute_rbh_pairs_between(sp_a, faa_a, sp_b, faa_b, work_root, args):
    """Genuine two-genome reciprocal best hit between any two species (not
    tied to a 'reference' role for either one). Returns list of
    (gene_in_a, gene_in_b) pairs."""
    pair_dir = work_root / f"{sp_a}_vs_{sp_b}"
    pair_dir.mkdir(parents=True, exist_ok=True)
    fwd_tsv = run_diamond_oneway(faa_a, faa_b, pair_dir, f"{sp_a}_vs_{sp_b}",
                                  args.evalue, args.max_target_seqs, args.threads,
                                  reuse_blast=not args.no_reuse_blast)
    rev_tsv = run_diamond_oneway(faa_b, faa_a, pair_dir, f"{sp_b}_vs_{sp_a}",
                                  args.evalue, args.max_target_seqs, args.threads,
                                  reuse_blast=not args.no_reuse_blast)
    return compute_rbh(fwd_tsv, rev_tsv)


# ==========================================================================
# Whole-genome/segmental duplication (ohnolog) detection and rescue
# (--detect-wgd)
#
# Reciprocal BEST hit is, by construction, a 1:1 mapping: a query gene
# points to exactly one subject (best_hits() keeps only the top-bitscore row
# per query). If the species being compared to the reference/core set has an
# extra whole-genome or segmental duplication that the reference does not
# share (e.g. the ancestral arachnopulmonate WGD when comparing a spider to
# a non-duplicated outgroup such as a tick), every ancestral gene generally
# has TWO surviving copies (ohnologs) in that species, but cross-species RBH
# can only ever recover one of them - the other is silently absent from the
# orthology table, not merged or miscounted. Downstream this can either
# dilute an ALG's ortholog count below the significance threshold on BOTH
# ohnologous chromosomes (a false negative - the duplication goes
# undetected) or, less often, contribute to a spurious fusion call if the
# "wrong" ohnolog wins for different genes in a way that happens to line up
# with a different linkage group on the same chromosome.
#
# --detect-wgd fixes this at the source rather than patching the symptom:
# for each requested species, an intragenomic self-comparison finds
# reciprocal-best-hit PARALOG pairs; pairs landing on two different
# scaffolds are aggregated per scaffold-pair and tested for enrichment
# (same Fisher's-exact logic as the ALG-scaffold association test) to call
# genuine ohnologous chromosome/scaffold pairs, as distinct from incidental
# one-off paralogy. Once a scaffold pair is confirmed ohnologous, any
# accepted cross-species RBH pair whose species-side gene has a paralog on
# the confirmed partner scaffold has that paralog added back in as a second,
# explicitly flagged ortholog of the same reference gene - recovering
# exactly the copy plain RBH structurally drops.
#
# Scope: rescue is only applied on the "being compared" side of each
# reference/core-vs-species search (i.e. it assumes the reference genome or
# core set itself is not independently duplicated). A duplication private to
# the reference/core genomes is not handled by this mechanism; see the
# script docstring.
# ==========================================================================

def best_hits_excluding_self(blast_tsv):
    """Like best_hits(), but for a self-vs-self search: drops the trivial
    self-hit (qseqid == sseqid) before picking each query's top subject."""
    df = pd.read_csv(blast_tsv, sep="\t", names=_BLAST_COLS)
    if df.empty:
        return {}
    df = df[df["qseqid"] != df["sseqid"]]
    if df.empty:
        return {}
    best = df.sort_values(["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid")
    return dict(zip(best["qseqid"], best["sseqid"]))


def compute_self_paralog_pairs(self_blast_tsv):
    """Intragenomic reciprocal-best-hit paralog pairs from a single self-vs-
    self diamond search (one search suffices: reciprocity is checked within
    the same best-hit table). Returns list of (gene1, gene2)."""
    best = best_hits_excluding_self(self_blast_tsv)
    pairs = []
    seen = set()
    for g1, g2 in best.items():
        if best.get(g2) == g1:
            pair = tuple(sorted((g1, g2)))
            if pair not in seen:
                seen.add(pair)
                pairs.append(pair)
    return pairs


def detect_ohnologous_chromosomes(paralog_pairs, gene_to_chr, min_pairs, alpha):
    """Aggregate intragenomic paralog pairs into candidate ohnologous
    scaffold/chromosome pairs. Same-scaffold paralogs are excluded (more
    likely tandem duplicates than whole-genome/segmental duplication
    survivors). A scaffold pair is confirmed if it is supported by at least
    `min_pairs` paralog pairs and, where scipy is available, a Fisher's
    exact enrichment test (analogous to the ALG-scaffold whole_FET test)
    does not reject it at `alpha`.
    Returns (ohnolog_of: dict[scaffold] -> set[scaffold], report: list of
    (scaffold_a, scaffold_b, n_paralog_pairs, p_value_or_None))."""
    pair_counts = Counter()
    for g1, g2 in paralog_pairs:
        c1 = lookup_chr(g1, gene_to_chr)
        c2 = lookup_chr(g2, gene_to_chr)
        if c1 is None or c2 is None or c1 == c2:
            continue
        pair_counts[tuple(sorted((c1, c2)))] += 1

    total = sum(pair_counts.values())
    chrom_totals = Counter()
    for (c1, c2), n in pair_counts.items():
        chrom_totals[c1] += n
        chrom_totals[c2] += n

    ohnolog_of = defaultdict(set)
    report = []
    for (c1, c2), n in sorted(pair_counts.items(), key=lambda kv: -kv[1]):
        if n < min_pairs:
            continue
        pval = None
        # The enrichment test needs some OTHER paralog-pair evidence elsewhere in the
        # genome to serve as a statistical background; when this candidate pair is the
        # only paralogy signal at all (total == n, i.e. a degenerate table with a zero
        # margin), scipy correctly reports p=1.0 for that table, which would wrongly
        # reject the single strongest (and only) candidate in the genome. In that case
        # fall back to the min_pairs count alone, which the pair already cleared above.
        if fisher_exact is not None and total > n:
            a = n
            b = chrom_totals[c1] - n
            c = chrom_totals[c2] - n
            d = total - a - b - c
            if min(b, c, d) >= 0:
                try:
                    _, pval = fisher_exact([[a, b], [c, d]])
                except Exception:
                    pval = None
        if pval is not None and pval > alpha:
            continue
        ohnolog_of[c1].add(c2)
        ohnolog_of[c2].add(c1)
        report.append((c1, c2, n, pval))
    return ohnolog_of, report


def run_wgd_detection(sp_name, faa, gene_to_chr, work_root, args):
    """Self-blast a species against itself and detect confirmed ohnologous
    scaffold pairs. Returns (paralog_partners: dict[gene]->list[gene],
    ohnolog_of: dict[scaffold]->set[scaffold], report list)."""
    log(f"=== WGD/ohnolog detection: {sp_name} (self-vs-self) ===")
    self_dir = work_root / f"{sp_name}_self"
    self_dir.mkdir(parents=True, exist_ok=True)
    self_tsv = run_diamond_oneway(faa, faa, self_dir, f"{sp_name}_self",
                                   args.evalue, args.max_target_seqs + 1, args.threads,
                                   reuse_blast=not args.no_reuse_blast)
    paralog_pairs = compute_self_paralog_pairs(self_tsv)
    log(f"  -> {len(paralog_pairs)} intragenomic reciprocal-best-hit paralog pairs")

    paralog_partners = defaultdict(list)
    for g1, g2 in paralog_pairs:
        paralog_partners[g1].append(g2)
        paralog_partners[g2].append(g1)

    ohnolog_of, report = detect_ohnologous_chromosomes(
        paralog_pairs, gene_to_chr, args.wgd_min_paralog_pairs, args.wgd_alpha)
    if report:
        log(f"  -> {len(report)} confirmed ohnologous scaffold pair(s):")
        for c1, c2, n, pval in report:
            pstr = f"{pval:.2e}" if pval is not None else "n/a"
            log(f"       {c1} <-> {c2}  ({n} paralog pairs, p={pstr})")
    else:
        log(f"  -> no scaffold pair reached --wgd-min-paralog-pairs "
            f"({args.wgd_min_paralog_pairs}) / --wgd-alpha ({args.wgd_alpha}); "
            f"no evidence of retained whole-genome/segmental duplication found")
    return dict(paralog_partners), ohnolog_of, report


def rescue_ohnolog_pairs(pairs, sp_paralog_partners, sp_gene_to_chr, sp_ohnolog_of):
    """pairs: list of (ref_gene, sp_gene) accepted RBH pairs for one
    reference/core-vs-species comparison. For each accepted pair, if the
    species-side gene has a paralog on a CONFIRMED ohnologous scaffold that
    is not already claimed as some other pair's direct hit, add it back in
    as a second, explicitly flagged ortholog of the same reference gene.
    Returns (extended: list of (ref_gene, sp_gene, is_inferred), n_rescued)."""
    claimed = {sp_gene for _, sp_gene in pairs}
    extended = [(ref_gene, sp_gene, False) for ref_gene, sp_gene in pairs]
    n_rescued = 0
    for ref_gene, sp_gene in pairs:
        c1 = lookup_chr(sp_gene, sp_gene_to_chr)
        for partner in sp_paralog_partners.get(sp_gene, []):
            if partner in claimed:
                continue
            c2 = lookup_chr(partner, sp_gene_to_chr)
            if c1 is None or c2 is None or c2 not in sp_ohnolog_of.get(c1, set()):
                continue
            extended.append((ref_gene, partner, True))
            claimed.add(partner)
            n_rescued += 1
    return extended, n_rescued


# ==========================================================================
# Multi-genome consensus linkage-group reconstruction (--consensus-species)
# ==========================================================================

class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def build_core_orthogroups(core_species, faa_of, gene_to_chr_of, gene_to_pos_of, work_root, args):
    """All-vs-all pairwise RBH among the core species -> orthogroups
    (connected components of the RBH graph). Returns dict
    og_id -> {species: (gene, chrom, pos)}, restricted to components with
    exactly one gene per species that has one (components with 2+ genes
    from the same core species are dropped as ambiguous/paralogous)."""
    uf = UnionFind()
    n_edges = 0
    for sp_a, sp_b in combinations(core_species, 2):
        log(f"=== core RBH: {sp_a} vs {sp_b} ===")
        pairs = compute_rbh_pairs_between(sp_a, faa_of[sp_a], sp_b, faa_of[sp_b], work_root, args)
        log(f"  -> {len(pairs)} reciprocal best hits")
        for gene_a, gene_b in pairs:
            uf.union((sp_a, gene_a), (sp_b, gene_b))
            n_edges += 1
    if n_edges == 0:
        sys.exit("No reciprocal best hits found among any pair of --consensus-species; "
                 "cannot build orthogroups.")

    components = defaultdict(list)
    for node in uf.parent:
        components[uf.find(node)].append(node)

    orthogroups = {}
    ambiguous = 0
    unresolved = 0
    for members in components.values():
        by_species = defaultdict(list)
        for sp, gene in members:
            by_species[sp].append(gene)
        if any(len(genes) > 1 for genes in by_species.values()):
            ambiguous += 1
            continue
        record = {}
        for sp, genes in by_species.items():
            gene = genes[0]
            chrom = lookup_chr(gene, gene_to_chr_of[sp])
            pos = lookup_pos(gene, gene_to_pos_of[sp])
            if chrom is None:
                continue
            record[sp] = (gene, chrom, pos)
        if not record:
            unresolved += 1
            continue
        og_id = f"og{len(orthogroups)}"
        orthogroups[og_id] = record

    log(f"built {len(orthogroups)} orthogroups from {len(components)} RBH-graph components "
        f"({ambiguous} dropped as ambiguous/multi-copy, {unresolved} dropped for unresolvable "
        f"GFF coordinates)")
    return orthogroups


def reconstruct_consensus_algs(core_species, orthogroups, min_species, min_frac, min_alg_orthogroups=1):
    """Group orthogroups into linkage groups wherever they stay co-located
    on the same chromosome across (>= min_frac fraction of) the core
    genomes that have both. Returns dict og_id -> alg_id string.

    Gene-level --min-scaffold-length filtering (applied upstream, before
    orthology search) already keeps genes on short/fragmentary scaffolds
    out of the orthogroups entirely. That alone does not bound the number
    of *reconstructed linkage groups*, though: even with every gene coming
    from a scaffold that individually clears the length filter, imperfect
    cross-genome agreement (default: identical chromosome assignment in
    every core genome) can still split what should be one ancestral unit
    into many small groups, especially with more than a couple of core
    genomes. min_alg_orthogroups is a direct, final filter on the
    reconstruction's own output: any linkage group with fewer than that
    many member orthogroups is dropped outright (its member genes are
    simply absent from every downstream .rbh file, exactly like a gene on
    a too-short scaffold), which is the most reliable way to collapse a
    reconstruction that came out as hundreds of small, likely-spurious
    linkage groups down to the handful of well-supported ones."""
    kept = {og: recs for og, recs in orthogroups.items() if len(recs) >= min_species}
    log(f"{len(kept)}/{len(orthogroups)} orthogroups present in >= {min_species} core species")
    if not kept:
        sys.exit("No orthogroups survive --consensus-min-species; lower it or add more/better-"
                 "conserved --consensus-species.")

    if min_species == len(core_species) and min_frac >= 1.0:
        # fast path: every kept orthogroup has all core species, so "co-located in
        # every genome that has both" reduces to "identical per-genome chromosome tuple"
        # - a direct grouping by key, with no transitivity concern.
        groups = defaultdict(list)
        for og, recs in kept.items():
            key = tuple(recs[sp][1] for sp in core_species)
            groups[key].append(og)
    else:
        # general path: agglomerative clustering with cluster-level re-validation at
        # every merge, NOT naive pairwise union-find. This distinction matters: two
        # orthogroups whose *individual* pairwise agreement clears min_frac can still
        # end up transitively chained, via a shared "bridge" orthogroup, into the same
        # cluster as a third orthogroup they have essentially nothing to do with -
        # plain union-find on pairwise votes has no protection against this, and with
        # thousands of orthogroups and a permissive-enough threshold it reliably
        # collapses everything into one giant, meaningless cluster (observed directly:
        # min_frac=0.6 on a real 6-genome, ~4100-orthogroup dataset produced exactly
        # one "linkage group" containing every orthogroup). Requiring every merge to be
        # re-justified against each cluster's own growing majority vote, not just the
        # single original pair that triggered it, is what actually prevents this: a
        # cluster's consensus gets harder to accidentally match as it grows, rather
        # than easier, which is the opposite of what plain transitive closure gives you.
        per_species_by_chrom = defaultdict(lambda: defaultdict(list))
        for og, recs in kept.items():
            for sp, (gene, chrom, pos) in recs.items():
                per_species_by_chrom[sp][chrom].append(og)

        pair_same = defaultdict(int)
        for sp, chrom_groups in per_species_by_chrom.items():
            for chrom, ogs in chrom_groups.items():
                if len(ogs) > 2000:
                    log(f"WARNING: {sp} chromosome {chrom} has {len(ogs)} orthogroups; "
                        f"pairwise consensus checking may be slow")
                for og1, og2 in combinations(sorted(set(ogs)), 2):
                    pair_same[(og1, og2)] += 1

        presence = {og: set(recs.keys()) for og, recs in kept.items()}
        uf = UnionFind()
        for og in kept:
            uf.find(og)

        # cluster_votes[cluster_root][species] = Counter(chromosome -> member count),
        # i.e. each cluster's own per-species majority chromosome, updated incrementally
        # as clusters merge rather than recomputed from scratch each time.
        cluster_votes = defaultdict(lambda: defaultdict(Counter))
        for og, recs in kept.items():
            for sp, (gene, chrom, pos) in recs.items():
                cluster_votes[og][sp][chrom] += 1

        def cluster_majority(root, sp):
            votes = cluster_votes[root].get(sp)
            return votes.most_common(1)[0][0] if votes else None

        def cluster_agreement(root1, root2):
            shared = set(cluster_votes[root1]) & set(cluster_votes[root2])
            if len(shared) < 2:
                return None
            same = sum(1 for sp in shared if cluster_majority(root1, sp) == cluster_majority(root2, sp))
            return same / len(shared)

        # Process candidate pairs most-confident-first: only pairs that already share a
        # chromosome in at least one genome ever appear here at all (pair_same is built
        # only from such co-occurrences), so this is already restricted to plausible
        # candidates, not a full O(n^2) scan.
        ordered_pairs = sorted(
            pair_same.keys(),
            key=lambda pr: -(pair_same[pr] / max(1, len(presence[pr[0]] & presence[pr[1]]))))

        for og1, og2 in ordered_pairs:
            r1, r2 = uf.find(og1), uf.find(og2)
            if r1 == r2:
                continue
            agreement = cluster_agreement(r1, r2)
            if agreement is None or agreement < min_frac:
                continue
            uf.union(r1, r2)
            new_root = uf.find(r1)
            old_root = r2 if new_root == r1 else r1
            for sp, counts in cluster_votes[old_root].items():
                cluster_votes[new_root][sp].update(counts)
            del cluster_votes[old_root]

        groups = defaultdict(list)
        for og in kept:
            groups[uf.find(og)].append(og)

    n_groups_before = len(groups)
    # full size distribution, before filtering, so the threshold can be picked from data
    # rather than guessed - same idea as duplication_blocks.tsv for --wgd-min-paralog-pairs
    size_report = sorted((len(ogs) for ogs in groups.values()), reverse=True)

    if min_alg_orthogroups > 1:
        small = {key: ogs for key, ogs in groups.items() if len(ogs) < min_alg_orthogroups}
        if small:
            n_dropped_ogs = sum(len(ogs) for ogs in small.values())
            log(f"  -> --consensus-min-alg-orthogroups {min_alg_orthogroups}: dropping "
                f"{len(small)} linkage group(s) with fewer than that many orthogroups "
                f"({n_dropped_ogs} orthogroups excluded)")
            groups = {key: ogs for key, ogs in groups.items() if len(ogs) >= min_alg_orthogroups}
        if not groups:
            sys.exit(f"--consensus-min-alg-orthogroups {min_alg_orthogroups} drops every "
                     f"reconstructed linkage group (largest had "
                     f"{max((len(v) for v in small.values()), default=0)} orthogroups); lower it.")

    alg_of = {}
    for rank, ogs in enumerate(sorted(groups.values(), key=len, reverse=True), start=1):
        for og in ogs:
            alg_of[og] = f"cALG{rank}"
    log(f"reconstructed {len(set(alg_of.values()))} consensus linkage groups "
        f"({n_groups_before} before size filtering) from {len(alg_of)} orthogroups")
    return alg_of, size_report

def run_consensus_mode(species, core_species, faa_of, gene_to_chr_of, gene_to_pos_of,
                        work_root, args, rbh_outdir, wgd_report=None):
    orthogroups = build_core_orthogroups(core_species, faa_of, gene_to_chr_of, gene_to_pos_of,
                                          work_root, args)
    alg_of, lg_size_report = reconstruct_consensus_algs(core_species, orthogroups, args.consensus_min_species,
                                                         args.consensus_min_frac, args.consensus_min_alg_orthogroups)

    sizes_path = rbh_outdir / "consensus_linkage_group_sizes.tsv"
    with open(sizes_path, "w") as f:
        f.write("candidate_rank\tn_member_orthogroups\tkept\n")
        for rank, n in enumerate(lg_size_report, start=1):
            kept = "yes" if n >= args.consensus_min_alg_orthogroups else "no"
            f.write(f"{rank}\t{n}\t{kept}\n")
    log(f"wrote {sizes_path} ({len(lg_size_report)} candidate linkage groups, before size filtering)")

    anchor = core_species[0]
    ref_pos_of_og = {}
    running = 0
    for og, recs in orthogroups.items():
        if og not in alg_of:
            continue
        if anchor in recs:
            ref_pos_of_og[og] = recs[anchor][2]
        else:
            running += 1
            ref_pos_of_og[og] = running
    ref_label = "consensusALG_from_" + "_".join(core_species)

    # reverse lookup: (core species, its own gene id) -> orthogroup id, for mapping
    # non-core species' RBH hits back into an orthogroup/ALG
    core_gene_to_og = {sp: {} for sp in core_species}
    for og, recs in orthogroups.items():
        if og not in alg_of:
            continue
        for sp, (gene, chrom, pos) in recs.items():
            if sp in core_gene_to_og:
                core_gene_to_og[sp][gene] = og

    written = []
    for sp_name in species:
        rows = []
        if sp_name in core_species:
            for og, recs in orthogroups.items():
                if og not in alg_of or sp_name not in recs:
                    continue
                gene, chrom, pos = recs[sp_name]
                rows.append({
                    "ref_gene": og, "ref_chrom": alg_of[og], "ref_pos": ref_pos_of_og[og],
                    "sp_gene": gene, "sp_scaf": chrom, "sp_pos": pos,
                })
        else:
            log(f"=== mapping {sp_name} onto consensus linkage groups ===")
            sp_gene_to_chr = gene_to_chr_of[sp_name]
            sp_gene_to_pos = gene_to_pos_of[sp_name]

            sp_paralog_partners, sp_ohnolog_of = {}, {}
            if args.detect_wgd and (not args.wgd_species or sp_name in args.wgd_species):
                sp_paralog_partners, sp_ohnolog_of, report = run_wgd_detection(
                    sp_name, faa_of[sp_name], sp_gene_to_chr, work_root, args)
                if wgd_report is not None:
                    wgd_report[sp_name] = report

            votes = defaultdict(list)  # noncore_gene -> [(alg, og), ...] one per core species hit
            is_inferred = {}
            n_rescued_total = 0
            for core_sp in core_species:
                pairs = compute_rbh_pairs_between(core_sp, faa_of[core_sp], sp_name, faa_of[sp_name],
                                                   work_root, args)
                if sp_paralog_partners:
                    extended, n_rescued = rescue_ohnolog_pairs(pairs, sp_paralog_partners,
                                                                 sp_gene_to_chr, sp_ohnolog_of)
                    n_rescued_total += n_rescued
                    for core_gene, nc_gene, inf in extended:
                        is_inferred[nc_gene] = is_inferred.get(nc_gene, False) or inf
                    pairs = [(cg, ng) for cg, ng, _ in extended]
                for core_gene, nc_gene in pairs:
                    og = core_gene_to_og[core_sp].get(core_gene)
                    if og is None or og not in alg_of:
                        continue
                    votes[nc_gene].append((alg_of[og], og))
            if n_rescued_total:
                log(f"  -> rescued {n_rescued_total} ohnolog hit(s) dropped by plain reciprocal-best-hit")
            conflicts = 0
            for nc_gene, vs in votes.items():
                algs = [a for a, _ in vs]
                chosen_alg, n_votes = Counter(algs).most_common(1)[0]
                if len(set(algs)) > 1:
                    conflicts += 1
                chosen_og = next(og for a, og in vs if a == chosen_alg)
                chrom = lookup_chr(nc_gene, sp_gene_to_chr)
                pos = lookup_pos(nc_gene, sp_gene_to_pos)
                if chrom is None:
                    continue
                rows.append({
                    "ref_gene": chosen_og, "ref_chrom": chosen_alg, "ref_pos": ref_pos_of_og[chosen_og],
                    "sp_gene": nc_gene, "sp_scaf": chrom, "sp_pos": pos,
                    "orthology_type": "wgd_inferred" if is_inferred.get(nc_gene) else "direct",
                })
            if conflicts:
                log(f"  -> {conflicts}/{len(votes)} genes had conflicting core-species hits; "
                    f"kept the majority-vote linkage group for each")
            log(f"  -> {len(rows)} orthologs mapped for {sp_name}")

        if not rows:
            log(f"  -> WARNING: no rows for {sp_name}; skipping (no .rbh file written).")
            continue

        df = pd.DataFrame(rows)
        fet_lookup = compute_whole_fet(df[["sp_scaf", "ref_chrom"]])
        for r in rows:
            r["fet"] = fet_lookup.get((r["sp_scaf"], r["ref_chrom"]))

        out_path = rbh_outdir / f"{sp_name}.rbh"
        write_rbh(out_path, ref_label, sp_name, rows)
        log(f"  -> wrote {out_path} ({len(rows)} orthologs)")
        written.append(out_path)

    return written


# ==========================================================================
# Enrichment (whole_FET) + .rbh writing
# ==========================================================================

def compute_whole_fet(df):
    """df has columns ref_chrom, sp_scaf. Returns dict (sp_scaf, ref_chrom) -> p-value."""
    if fisher_exact is None:
        log("  scipy not installed - skipping whole_FET enrichment test "
            "(alg_parsimony.py's --min-genes filtering still applies)")
        return {}
    n_total = len(df)
    scaf_counts = df["sp_scaf"].value_counts().to_dict()
    chrom_counts = df["ref_chrom"].value_counts().to_dict()
    pair_counts = df.groupby(["sp_scaf", "ref_chrom"]).size().to_dict()

    result = {}
    for (scaf, chrom), a in pair_counts.items():
        b = scaf_counts[scaf] - a
        c = chrom_counts[chrom] - a
        d = n_total - a - b - c
        try:
            _, p = fisher_exact([[a, b], [c, d]])
        except Exception:
            p = float("nan")
        result[(scaf, chrom)] = p
    return result


def color_for(name, palette=DEFAULT_PALETTE):
    """Deterministic across processes, platforms and Python versions (unlike
    the built-in hash() for strings, which is randomized per-process by
    default) - the SAME linkage-group name always maps to the SAME palette
    index, in genomes_to_rbh.py and in alg_parsimony.py alike, so colors
    stay comparable across separate analyses and across the two scripts."""
    idx = zlib.crc32(name.encode("utf-8")) % len(palette)
    return palette[idx]


def write_rbh(path, ref_name, sp_name, rows):
    """rows: list of dicts with keys ref_gene, ref_chrom, ref_pos, sp_gene, sp_scaf, sp_pos, fet,
    and optionally orthology_type ('direct' or 'wgd_inferred'; defaults to 'direct' if absent)."""
    cols = [
        "rbh", "gene_group", "color",
        f"{ref_name}_gene", f"{ref_name}_scaf", f"{ref_name}_pos",
        f"{sp_name}_gene", f"{sp_name}_scaf", f"{sp_name}_pos",
        f"{sp_name}_breakchrom", f"{sp_name}_ix", f"{sp_name}_break_ix",
        f"{ref_name}_ix", f"{ref_name}_break_ix", f"{ref_name}_breakchrom",
        "whole_FET", "break_FET",
        f"{ref_name}_D", f"{sp_name}_D", f"{sp_name}_plotindex", f"{sp_name}_plotpos",
        f"{ref_name}_plotindex", f"{ref_name}_plotpos", "orthology_type",
    ]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for i, r in enumerate(rows):
            record = {
                "rbh": f"rbh_{ref_name}_{sp_name}_{i}",
                "gene_group": r["ref_chrom"],
                "color": color_for(r["ref_chrom"]),
                f"{ref_name}_gene": r["ref_gene"],
                f"{ref_name}_scaf": r["ref_chrom"],
                f"{ref_name}_pos": r["ref_pos"] if r["ref_pos"] is not None else "",
                f"{sp_name}_gene": r["sp_gene"],
                f"{sp_name}_scaf": r["sp_scaf"],
                f"{sp_name}_pos": r["sp_pos"] if r["sp_pos"] is not None else "",
                f"{sp_name}_breakchrom": f"{r['sp_scaf']}:0-end",
                f"{sp_name}_ix": i, f"{sp_name}_break_ix": i,
                f"{ref_name}_ix": i, f"{ref_name}_break_ix": i,
                f"{ref_name}_breakchrom": f"{r['ref_chrom']}:0-end",
                "whole_FET": r["fet"] if r["fet"] is not None else "",
                "break_FET": r["fet"] if r["fet"] is not None else "",
                f"{ref_name}_D": "", f"{sp_name}_D": "",
                f"{sp_name}_plotindex": i, f"{sp_name}_plotpos": r["sp_pos"] or "",
                f"{ref_name}_plotindex": i, f"{ref_name}_plotpos": r["ref_pos"] or "",
                "orthology_type": r.get("orthology_type", "direct"),
            }
            f.write("\t".join(str(record[c]) for c in cols) + "\n")


def write_duplication_report(path, all_reports):
    """all_reports: dict species -> list of (scaf_a, scaf_b, n_pairs, pval)."""
    with open(path, "w") as f:
        f.write("species\tscaffold_a\tscaffold_b\tn_paralog_pairs\tp_value\n")
        for sp, report in all_reports.items():
            for c1, c2, n, pval in report:
                pstr = f"{pval:.6g}" if pval is not None else "NA"
                f.write(f"{sp}\t{c1}\t{c2}\t{n}\t{pstr}\n")


# ==========================================================================
# Per-species driver
# ==========================================================================

def process_species(sp_name, sp_info, ref_name, ref_faa, ref_gene_to_chr, ref_gene_to_pos,
                     work_root, wgd_output_dir, args, rbh_outdir, wgd_report=None):
    log(f"=== {ref_name} vs {sp_name} ===")
    sp_dir = work_root / sp_name
    sp_dir.mkdir(parents=True, exist_ok=True)

    sp_gene_to_chr = build_gene_to_chr(sp_info["gff"], sp_info["gff_is_gtf"])
    sp_gene_to_pos = build_gene_to_pos(sp_info["gff"], sp_info["gff_is_gtf"])
    sp_faa = get_search_ready_faa(sp_name, sp_info, sp_dir, wgd_output_dir, args.genetic_code,
                                   sp_gene_to_chr, args)

    fwd_tsv = run_diamond_oneway(ref_faa, sp_faa, sp_dir, f"{ref_name}_vs_{sp_name}",
                                  args.evalue, args.max_target_seqs, args.threads,
                                  reuse_blast=not args.no_reuse_blast)
    rev_tsv = run_diamond_oneway(sp_faa, ref_faa, sp_dir, f"{sp_name}_vs_{ref_name}",
                                  args.evalue, args.max_target_seqs, args.threads,
                                  reuse_blast=not args.no_reuse_blast)
    pairs = compute_rbh(fwd_tsv, rev_tsv)
    log(f"  -> {len(pairs)} reciprocal best hits")
    if not pairs:
        log(f"  -> WARNING: no RBH pairs for {sp_name}; skipping (no .rbh file written).")
        return None

    is_inferred = {}
    if args.detect_wgd and (not args.wgd_species or sp_name in args.wgd_species):
        paralog_partners, ohnolog_of, report = run_wgd_detection(sp_name, sp_faa, sp_gene_to_chr,
                                                                   work_root, args)
        if wgd_report is not None:
            wgd_report[sp_name] = report
        extended, n_rescued = rescue_ohnolog_pairs(pairs, paralog_partners, sp_gene_to_chr, ohnolog_of)
        if n_rescued:
            log(f"  -> rescued {n_rescued} ohnolog(s) dropped by plain reciprocal-best-hit")
        pairs = [(rg, sg) for rg, sg, _ in extended]
        is_inferred = {sg: inf for _, sg, inf in extended}

    rows = []
    dropped = 0
    for ref_gene, sp_gene in pairs:
        ref_chrom = lookup_chr(ref_gene, ref_gene_to_chr)
        sp_scaf = lookup_chr(sp_gene, sp_gene_to_chr)
        if ref_chrom is None or sp_scaf is None:
            dropped += 1
            continue
        rows.append({
            "ref_gene": ref_gene, "ref_chrom": ref_chrom, "ref_pos": lookup_pos(ref_gene, ref_gene_to_pos),
            "sp_gene": sp_gene, "sp_scaf": sp_scaf, "sp_pos": lookup_pos(sp_gene, sp_gene_to_pos),
            "orthology_type": "wgd_inferred" if is_inferred.get(sp_gene) else "direct",
        })
    if dropped:
        log(f"  -> WARNING: {dropped} RBH pairs dropped (gene ID missing from GFF/GTF)")
    if not rows:
        log(f"  -> WARNING: no pairs with resolvable coordinates for {sp_name}; skipping.")
        return None

    df = pd.DataFrame(rows)
    fet_lookup = compute_whole_fet(df[["sp_scaf", "ref_chrom"]])
    for r in rows:
        r["fet"] = fet_lookup.get((r["sp_scaf"], r["ref_chrom"]))

    out_path = rbh_outdir / f"{sp_name}.rbh"
    write_rbh(out_path, ref_name, sp_name, rows)
    log(f"  -> wrote {out_path} ({len(rows)} orthologs)")
    return out_path


# ==========================================================================
# Main
# ==========================================================================

def find_parsimony_script(explicit):
    if explicit:
        return Path(explicit)
    here = Path(__file__).resolve().parent / "alg_parsimony.py"
    if here.exists():
        return here
    cwd = Path.cwd() / "alg_parsimony.py"
    if cwd.exists():
        return cwd
    which = shutil.which("alg_parsimony.py")
    if which:
        return Path(which)
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Build .rbh files from raw genomes - using either a single reference "
                    "genome or a multi-genome consensus as the ALG proxy - then run "
                    "alg_parsimony.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input-dir", required=True, type=Path,
                        help="Directory with {name}.cds, {name}.gff3|gtf, [{name}.pep] per species")
    parser.add_argument("--wgd-output-dir", type=Path, default=None,
                        help="Optional: --output-dir of a prior run_wgd_pipeline.py run on the same "
                             "--input-dir, reused to skip re-translating CDS -> protein")

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--reference-species", default=None,
                        help="Single-reference mode: species to use as the ALG-proxy reference "
                             "(its chromosomes become the 'linkage groups'). If neither this nor "
                             "--consensus-species is given, the species with the fewest distinct "
                             "scaffolds is auto-picked as the reference and a warning is printed.")
    mode_group.add_argument("--consensus-species", nargs="+", default=None,
                        help="Multi-genome consensus mode: reconstruct linkage groups from "
                             "reciprocal best hits among these species (give >= 2 names, or the "
                             "single word 'all' to use every discovered species as the core set). "
                             "Every input genome - core or not - ends up with its own .rbh file.")
    parser.add_argument("--consensus-min-species", type=int, default=None,
                        help="Consensus mode: minimum number of core species an orthogroup must "
                             "appear in to be kept (default: all core species)")
    parser.add_argument("--consensus-min-frac", type=float, default=1.0,
                        help="Consensus mode: minimum fraction of core genomes (that have both) "
                             "in which two orthogroups must be co-located to link them into the "
                             "same linkage group (default: 1.0 = must hold in every one of them)")
    parser.add_argument("--consensus-min-alg-orthogroups", type=int, default=1,
                        help="Consensus mode: minimum number of member orthogroups a reconstructed "
                             "linkage group must have to be kept (default: 1 = no filtering). Raise "
                             "this to collapse a reconstruction that came out as hundreds of small "
                             "linkage groups - e.g. from fragmented core-genome assemblies where many "
                             "small-but-length-filter-passing scaffolds each land in their own group - "
                             "down to the well-supported ones; dropped groups' genes are simply absent "
                             "from every downstream .rbh file, same as a gene on a too-short scaffold.")

    parser.add_argument("--work-dir", type=Path, default=None,
                        help="Scratch directory for translations/diamond dbs/blast output "
                             "(default: <rbh-outdir>/_work)")
    parser.add_argument("--rbh-outdir", required=True, type=Path, help="Where to write .rbh files")
    parser.add_argument("--rbh-only", action="store_true",
                        help="Only build .rbh files; don't invoke alg_parsimony.py")
    parser.add_argument("--parsimony-script", default=None,
                        help="Path to alg_parsimony.py (auto-detected next to this script / cwd / "
                             "PATH if omitted)")
    parser.add_argument("--parsimony-outdir", type=Path, default=None,
                        help="--outdir to pass to alg_parsimony.py (default: <rbh-outdir>/parsimony_out)")
    parser.add_argument("--evalue", default="1e-3")
    parser.add_argument("--max-target-seqs", type=int, default=5)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--genetic-code", type=int, default=1)
    parser.add_argument("--no-reuse-blast", action="store_true",
                        help="Recompute every diamond alignment even if a cached result from a "
                             "previous run with the same --evalue/--max-target-seqs is present in "
                             "the work directory. By default such results ARE reused, which makes "
                             "repeated runs over the same genomes (e.g. sweeping the consensus "
                             "thresholds, which only affect post-alignment steps) much faster. Use "
                             "this if the input FASTA/annotation files themselves have changed but "
                             "the work directory was not cleared.")
    parser.add_argument("--min-scaffold-length", type=int, default=0,
                        help="Exclude genes on scaffolds shorter than this many bp from being considered "
                             "for reciprocal-best-hit search at all (default: 0 = no filtering). Scaffold "
                             "length is read from the GFF3 '##sequence-region' pragma when present, "
                             "otherwise approximated from the longest annotated feature's end coordinate "
                             "on that scaffold (a lower bound, since there's no assembly FASTA to measure "
                             "true length from). Applies to every species (reference/core and target alike).")
    parser.add_argument("--species", nargs="*", default=None,
                        help="Restrict to these species (plus the reference/core species); "
                             "default: all discovered")
    parser.add_argument("--detect-wgd", action="store_true",
                        help="Detect retained whole-genome/segmental duplication (ohnologs) in the "
                             "species being compared to the reference/core set via an intragenomic "
                             "self-blast, and rescue orthologs that plain reciprocal-best-hit would "
                             "otherwise silently drop (RBH is 1:1, so only one of two ohnologs ever "
                             "survives an ordinary cross-species search). Writes duplication_blocks.tsv "
                             "listing confirmed ohnologous scaffold pairs per species, and tags rescued "
                             "rows in each .rbh file's orthology_type column ('direct' vs 'wgd_inferred'). "
                             "Does not detect duplication in the reference genome or core species "
                             "themselves - see the script docstring.")
    parser.add_argument("--wgd-species", nargs="+", default=None,
                        help="Restrict --detect-wgd to these species only (default: every species "
                             "being compared to the reference/core set). Use this once you know which "
                             "taxa are the suspected polyploids/WGD lineages, to skip the extra "
                             "self-blast cost for everything else.")
    parser.add_argument("--wgd-min-paralog-pairs", type=int, default=5,
                        help="Minimum intragenomic paralog-pair count to call two scaffolds ohnologous "
                             "(default: 5)")
    parser.add_argument("--wgd-alpha", type=float, default=0.05,
                        help="Significance threshold for the ohnologous-scaffold-pair enrichment test "
                             "(default: 0.05)")
    parser.add_argument("--forward", nargs=argparse.REMAINDER, default=[],
                        help="Everything after this flag is passed through verbatim to "
                             "alg_parsimony.py, e.g. --forward --min-genes 3 --constraint '(...)'")
    args = parser.parse_args()

    args.rbh_outdir.mkdir(parents=True, exist_ok=True)
    work_root = args.work_dir or (args.rbh_outdir / "_work")
    work_root.mkdir(parents=True, exist_ok=True)

    species = discover_species(args.input_dir)
    if args.species:
        keep = set(args.species)
        if args.reference_species:
            keep.add(args.reference_species)
        if args.consensus_species and args.consensus_species != ["all"]:
            keep |= set(args.consensus_species)
        species = {k: v for k, v in species.items() if k in keep}
    if not species:
        sys.exit(f"No species discovered in {args.input_dir}")
    log(f"discovered {len(species)} species: {', '.join(species)}")

    try:
        check_tool("diamond")
    except RuntimeError as e:
        sys.exit(str(e))

    if args.consensus_species:
        core_species = list(species) if args.consensus_species == ["all"] else args.consensus_species
        unknown = set(core_species) - set(species)
        if unknown:
            sys.exit(f"--consensus-species mentions species not discovered in {args.input_dir}: "
                     f"{sorted(unknown)}. Discovered: {list(species)}")
        if len(core_species) < 2:
            sys.exit("--consensus-species needs at least 2 species")
        min_species = args.consensus_min_species or len(core_species)
        args.consensus_min_species = min_species
        log(f"consensus mode: core species = {core_species} "
            f"(min_species={min_species}, min_frac={args.consensus_min_frac})")

        faa_of = {}
        gene_to_chr_of = {}
        gene_to_pos_of = {}
        for sp_name, sp_info in species.items():
            sp_dir = work_root / sp_name
            sp_dir.mkdir(parents=True, exist_ok=True)
            gene_to_chr_of[sp_name] = build_gene_to_chr(sp_info["gff"], sp_info["gff_is_gtf"])
            gene_to_pos_of[sp_name] = build_gene_to_pos(sp_info["gff"], sp_info["gff_is_gtf"])
            faa_of[sp_name] = get_search_ready_faa(sp_name, sp_info, sp_dir, args.wgd_output_dir,
                                                     args.genetic_code, gene_to_chr_of[sp_name], args)

        wgd_report = {}
        written = run_consensus_mode(species, core_species, faa_of, gene_to_chr_of, gene_to_pos_of,
                                      work_root, args, args.rbh_outdir, wgd_report=wgd_report)
        if wgd_report:
            dup_path = args.rbh_outdir / "duplication_blocks.tsv"
            write_duplication_report(dup_path, wgd_report)
            log(f"wrote {dup_path}")
    else:
        ref_name = args.reference_species
        if ref_name is None:
            ref_name = min(species, key=lambda n: len(set(
                v for v in build_gene_to_chr(species[n]["gff"], species[n]["gff_is_gtf"]).values()
            )))
            log(f"WARNING: --reference-species not given; auto-picked '{ref_name}' (fewest distinct "
                f"scaffolds -> presumed most contiguous assembly). Override with --reference-species "
                f"(or use --consensus-species for a multi-genome reconstruction instead) if this "
                f"isn't the genome you want as your ALG-proxy reference.")
        if ref_name not in species:
            sys.exit(f"--reference-species '{ref_name}' not found among discovered species {list(species)}")

        ref_info = species[ref_name]
        ref_work_dir = work_root / ref_name
        ref_work_dir.mkdir(parents=True, exist_ok=True)
        ref_gene_to_chr = build_gene_to_chr(ref_info["gff"], ref_info["gff_is_gtf"])
        ref_gene_to_pos = build_gene_to_pos(ref_info["gff"], ref_info["gff_is_gtf"])
        ref_faa = get_search_ready_faa(ref_name, ref_info, ref_work_dir, args.wgd_output_dir,
                                        args.genetic_code, ref_gene_to_chr, args)
        log(f"reference species: {ref_name} ({len(set(ref_gene_to_chr.values()))} distinct scaffolds/chromosomes)")

        written = []
        wgd_report = {}
        for sp_name, sp_info in species.items():
            if sp_name == ref_name:
                continue
            try:
                out_path = process_species(sp_name, sp_info, ref_name, ref_faa, ref_gene_to_chr,
                                            ref_gene_to_pos, work_root, args.wgd_output_dir, args,
                                            args.rbh_outdir, wgd_report=wgd_report)
                if out_path:
                    written.append(out_path)
            except subprocess.CalledProcessError as e:
                log(f"  -> ERROR running diamond for {sp_name}: {e}. Skipping species.")
            except Exception as e:
                log(f"  -> ERROR processing {sp_name}: {e}. Skipping species.")
        if wgd_report:
            dup_path = args.rbh_outdir / "duplication_blocks.tsv"
            write_duplication_report(dup_path, wgd_report)
            log(f"wrote {dup_path}")

    if not written:
        sys.exit("No .rbh files were written for any species - nothing to run alg_parsimony.py on.")
    log(f"wrote {len(written)} .rbh file(s) to {args.rbh_outdir}")

    if args.rbh_only:
        return

    script_path = find_parsimony_script(args.parsimony_script)
    if script_path is None or not script_path.exists():
        sys.exit("Could not find alg_parsimony.py (looked next to this script, in the current "
                 "directory, and on PATH). Pass --parsimony-script /path/to/alg_parsimony.py, "
                 "or rerun with --rbh-only and invoke it yourself.")

    parsimony_outdir = args.parsimony_outdir or (args.rbh_outdir / "parsimony_out")
    cmd = [sys.executable, str(script_path), "--indir", str(args.rbh_outdir),
           "--outdir", str(parsimony_outdir), *args.forward]
    log("handing off to alg_parsimony.py:")
    run_cmd(cmd)


if __name__ == "__main__":
    main()
