#!/usr/bin/env python

import argparse

import numpy as np
import snapatac2 as snap
import utils
from __init__ import ASSAY

MAX_CELL = 10**6


def get_umi_count(bcs, umis, cell_barcodes, sample):
    """
    Args:
        bcs: raw barcodes
        umis: umi count
        cell_barcodes: cell barcodes
    """
    a = [(umi, bc) for umi, bc in zip(umis, bcs) if umi > 0]
    a.sort(reverse=True)
    cell_barcodes = set(cell_barcodes)
    plot_data = {}
    n = len(a)
    first_noncell = n - 1
    for i, (umi, bc) in enumerate(a):
        if bc not in cell_barcodes:
            first_noncell = i
            break
    print(f"first non-cell barcode rank: {first_noncell}")
    last_cell = 0
    for i in range(min(n - 1, MAX_CELL), -1, -1):
        bc = a[i][1]
        if bc in cell_barcodes:
            last_cell = i
            break
    pure = sample + ".cells.pure" + f"({first_noncell}/{first_noncell}, 100%)"
    bg_cells = n - first_noncell
    bg = sample + ".cells.background" + f"(0/{bg_cells}, 0%)"
    plot_data[pure] = {}
    plot_data[bg] = {}
    for i in range(first_noncell):
        plot_data[pure][i + 1] = int(a[i][0])

    n_mix = last_cell - first_noncell + 1
    if n_mix != 0:
        n_total = len(cell_barcodes)
        n_mix_cell = n_total - first_noncell
        mix_rate = round(n_mix_cell / n_mix * 100, 2)
        mix = sample + ".cells.mix" + f"({n_mix_cell}/{n_mix}, {mix_rate}%)"
        plot_data[mix] = {}
        for i in range(first_noncell, last_cell + 1):
            plot_data[mix][i + 1] = int(a[i][0])

    for i in range(last_cell + 1, min(MAX_CELL, n), 10):
        plot_data[bg][i + 1] = int(a[i][0])
    # do not record every umi count
    for i in range(MAX_CELL, n, 1000):
        plot_data[bg][i + 1] = int(a[i][0])
    return plot_data


def run(args):
    sample = args.sample
    stats = {}
    # genome
    genome = snap.genome.Genome(fasta=args.fasta, annotation=args.gtf)

    # read fragment
    filter_h5ad = f"{sample}.filtered.h5ad"
    data = snap.pp.import_data(
        args.fragment,
        chrom_sizes=genome,
        file=filter_h5ad,  # Optional
        sorted_by_barcode=False,
        min_num_fragments=0,
        chrM=["chrM", "M", "MT", "Mt"],
    )

    # tss enrich
    snap.metrics.tsse(data, genome)

    # umi count
    bcs = list(data.obs_names)
    umis = list(data.obs["n_fragment"])

    # filter cell
    frag = list(data.obs["n_fragment"])
    is_test = len(frag) < 10**6
    frag.sort(reverse=True)
    min_counts = frag[100] // 10
    min_tsse = 10
    if is_test:
        min_tsse = 2
    snap.pp.filter_cells(data, min_counts=min_counts, min_tsse=min_tsse)
    cell_barcodes = list(data.obs_names)

    # metrics
    n_cell = data.shape[0]
    stats["Estimated Number of Cells"] = n_cell
    median_fragment = np.median(data.obs["n_fragment"])
    stats["Median Fragments per Cell"] = median_fragment
    median_tsse = np.median(data.obs["tsse"])
    stats["Median TSS Enrichment Score per Cell"] = median_tsse
    frac_dup = utils.get_frac(np.median(data.obs["frac_dup"]))
    stats["Median Fraction of Duplicated Reads per Cell"] = frac_dup
    print(stats)

    # umi count json
    plot_data = get_umi_count(bcs, umis, cell_barcodes, sample)
    utils.write_multiqc(plot_data, sample, ASSAY, "umi_count")

    snap.pp.add_tile_matrix(data)

    snap.pp.select_features(data)
    weighted_by_sd = True
    if is_test:
        weighted_by_sd = False
    snap.tl.spectral(data, weighted_by_sd=weighted_by_sd)
    snap.tl.umap(data)
    snap.pp.knn(data)
    snap.tl.leiden(data)
    umap_fn = f"{sample}.umap.png"
    snap.pl.umap(data, color="leiden", interactive=False, show=False, height=500, out_file=umap_fn)

    # gene matrix
    snap.pp.make_gene_matrix(data, genome, file=f"{sample}.gene_mat.h5ad")

    # cluster peaks
    snap.tl.macs3(data, groupby="leiden")
    peaks = snap.tl.merge_peaks(data.uns["macs3"], genome)
    snap.pp.make_peak_matrix(data, use_rep=peaks["Peaks"], file=f"{sample}.peak_mat.h5ad")

    # peak metrics
    snap.metrics.frip(data, {"peaks_frac": peaks["Peaks"]})
    median_frac = utils.get_frac(np.median(data.obs["peaks_frac"]))
    stats["Median Fraction of Reads in Peaks per Cell"] = median_frac
    stats["Number of Peaks"] = peaks.shape[0]

    # write multiqc
    utils.write_multiqc(stats, sample, ASSAY, "snapatac2" + ".stats")


def main():
    """
    fragment,fasta,gtf
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument("--fragment", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    args = parser.parse_args()

    run(args)


if __name__ == "__main__":
    main()
