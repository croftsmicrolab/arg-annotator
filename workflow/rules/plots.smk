rule plot_insert_rank_curve:
    input:
        data=output_path("insert_cluster_summary"),
        script=script("plot_rank_curve.py")
    output:
        data=output_path("insert_cluster_distribution_data"),
        png=output_path("insert_cluster_distribution_png"),
        svg=output_path("insert_cluster_distribution_svg"),
        pdf=output_path("insert_cluster_distribution_pdf")
    log:
        output_path("insert_distribution_log")
    params:
        python=tool("python"),
        cache=output_path("matplotlib_cache_dir"),
        dpi=config["plots"]["dpi"],
        color=config["plots"]["line_color"],
        title=config["project"]["sample"] + " raw sequence cluster size curve"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})" {params.cache:q}
        MPLCONFIGDIR={params.cache:q} {params.python:q} {input.script:q} \
            --input {input.data:q} --id-column cluster_id --value-column size \
            --title {params.title:q} --x-label "Cluster rank (sorted by size)" \
            --y-label "Raw sequences per cluster (log scale)" \
            --data-output {output.data:q} --png-output {output.png:q} \
            --svg-output {output.svg:q} --pdf-output {output.pdf:q} \
            --dpi {params.dpi} --color {params.color:q} > {log:q} 2>&1
        """

rule plot_weighted_reporf_rank_curve:
    input:
        data=output_path("reporf_support"),
        script=script("plot_rank_curve.py")
    output:
        data=output_path("reporf_support_distribution_data"),
        png=output_path("reporf_support_distribution_png"),
        svg=output_path("reporf_support_distribution_svg"),
        pdf=output_path("reporf_support_distribution_pdf")
    log:
        output_path("support_distribution_log")
    params:
        python=tool("python"),
        cache=output_path("matplotlib_cache_dir"),
        dpi=config["plots"]["dpi"],
        color=config["plots"]["line_color"],
        title=config["project"]["sample"] + " abundance-weighted RepORF rank curve"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})" {params.cache:q}
        MPLCONFIGDIR={params.cache:q} {params.python:q} {input.script:q} \
            --input {input.data:q} --id-column reporf_id --value-column reporf_support \
            --title {params.title:q} --x-label "RepORF rank (sorted by weighted support)" \
            --y-label "Abundance-weighted RepORF support (log scale)" \
            --data-output {output.data:q} --png-output {output.png:q} \
            --svg-output {output.svg:q} --pdf-output {output.pdf:q} \
            --dpi {params.dpi} --color {params.color:q} > {log:q} 2>&1
        """
