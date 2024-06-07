#!/usr/bin/env python

import argparse

import numpy as np
import snapatac2 as snap
import utils
from __init__ import ASSAY


def run(args):
    sample = args.sample
    stats = {}
    # genome
    genome = snap.genome.Genome(fasta=args.fasta, annotation=args.gtf)

    # read fragment
    h5ad = f"{sample}.h5ad"
    data = snap.pp.import_data(
        args.fragment,
        chrom_sizes=genome,
        file=h5ad,  # Optional
        sorted_by_barcode=False,
        min_num_fragments=0,
        chrM=["chrM", "M", "MT", "Mt"],
    )

    # tss enrich
    snap.metrics.tsse(data, genome)

    # filter cell
    frag = list(data.obs["n_fragment"])
    is_test = len(frag) < 10**6
    frag.sort(reverse=True)
    min_counts = frag[100] // 10
    min_tsse = 10
    if is_test:
        min_tsse = 2
    snap.pp.filter_cells(data, min_counts=min_counts, min_tsse=min_tsse)

    # metrics
    n_cell = data.shape[0]
    stats["Estimated Number of Cells"] = n_cell
    median_fragment = np.median(data.obs["n_fragment"])
    stats["Median Fragments per Cell"] = median_fragment
    median_tsse = np.median(data.obs["tsse"])
    stats["Median TSS Enrichment Score per Cell"] = median_tsse
    print(stats)

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

    # cluster peaks
    snap.tl.macs3(data, groupby="leiden")
    peaks = snap.tl.merge_peaks(data.uns["macs3"], genome)
    snap.pp.make_peak_matrix(data, use_rep=peaks["Peaks"], file=f"{sample}.peak.h5ad")

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
