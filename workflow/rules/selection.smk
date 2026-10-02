rule select_reporf_candidates:
    input:
        support=output_path("reporf_support"),
        fna=output_path("predicted_orfs_fna"),
        faa=output_path("predicted_proteins_faa"),
        script=script("select_tie_aware_candidates.py")
    output:
        summary=output_path("selection_summary"),
        primary_table=output_path("primary_table"),
        primary_ids=output_path("primary_ids"),
        primary_fna=output_path("primary_fna"),
        primary_faa=output_path("primary_faa"),
        exploratory_table=output_path("exploratory_support1_table"),
        exploratory_ids=output_path("exploratory_support1_ids"),
        exploratory_fna=output_path("exploratory_support1_fna"),
        exploratory_faa=output_path("exploratory_support1_faa")
    log:
        output_path("selection_log")
    params:
        python=tool("python"),
        fraction=config["selection"]["nominal_fraction"],
        mode=config["selection"]["mode"],
        minimum=config["selection"]["minimum_support"]
    shell:
        r"""
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {log:q})"
        {params.python:q} {input.script:q} \
            --table {input.support:q} \
            --mode {params.mode:q} \
            --fraction {params.fraction} \
            --minimum-primary-support {params.minimum} \
            --orf-fna {input.fna:q} \
            --orf-faa {input.faa:q} \
            --primary-table {output.primary_table:q} \
            --primary-ids {output.primary_ids:q} \
            --primary-fna {output.primary_fna:q} \
            --primary-faa {output.primary_faa:q} \
            --exploratory-table {output.exploratory_table:q} \
            --exploratory-ids {output.exploratory_ids:q} \
            --exploratory-fna {output.exploratory_fna:q} \
            --exploratory-faa {output.exploratory_faa:q} \
            --summary {output.summary:q} > {log:q} 2>&1
        """


rule sanitize_primary_proteins_for_rgi:
    input:
        data=output_path("primary_faa"),
        script=script("sanitize_protein_fasta.py")
    output:
        output_path("primary_rgi_faa")
    params:
        python=tool("python"),
    shell:
        "{params.python:q} {input.script:q} --input {input.data:q} --output {output:q} --allow-empty"
