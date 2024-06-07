def multiqc_sgr_config():
    from multiqc import config

    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "scrna/stats": {"fn": "*scrna.*stats.json"},
        "scrna/umi_count": {"fn": "*scrna.umi_count.json"},
        "scrna/saturation": {
            "fn": "*scrna.saturation.json",
        },
        "scrna/median_gene": {
            "fn": "*scrna.median_gene.json",
        },
        "scsnp/stats": {"fn": "*scsnp.*stats.json"},
        "scsnp/gene": {"fn": "*scsnp.gene.json"},
        "scsnp/count": {"fn": "*scsnp.count.json"},
        "scsnp/meta": {"fn": "*scsnp.meta.json"},
        "scatac/stats": {"fn": "*scatac.*stats.json"},
        "scatac/umi_count": {"fn": "*scatac.umi_count.json"},
    }
    config.update_dict(config.sp, sgr_search_patterns)

    config.update_dict(
        config.table_columns_visible,
        {
            "QualiMap": {
                "total_reads": False,
                "mapped_reads": False,
                "general_error_rate": False,
                "mean_coverage": False,
                "percentage_aligned": True,
                "median_coverage": False,
                "median_insert_size": False,
                "avg_gc": False,
                "1_x_pc": False,
                "5_x_pc": False,
                "10_x_pc": False,
                "30_x_pc": False,
                "50_x_pc": False,
            },
        },
    )
