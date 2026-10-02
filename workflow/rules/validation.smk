rule validate_complete_pipeline:
    input:
        input_fasta=input_path("fasta"),
        normalized=output_path("normalized_fasta"),
        id_mapping=output_path("id_mapping"),
        sequence_qc=output_path("sequence_qc"),
        insert_members=output_path("insert_cluster_members"),
        insert_summary=output_path("insert_cluster_summary"),
        insert_reps=output_path("insert_representatives"),
        orf_fna=output_path("predicted_orfs_fna"),
        orf_faa=output_path("predicted_proteins_faa"),
        orf_members=output_path("orf_cluster_members"),
        orf_summary=output_path("orf_cluster_summary"),
        representative_orfs=output_path("representative_orfs_fna"),
        mapping=output_path("orf_to_insert_mapping"),
        missing=output_path("missing_mappings"),
        support=output_path("reporf_support"),
        selection_summary=output_path("selection_summary"),
        primary_table=output_path("primary_table"),
        primary_fna=output_path("primary_fna"),
        primary_faa=output_path("primary_faa"),
        exploratory_table=output_path("exploratory_support1_table"),
        exploratory_fna=output_path("exploratory_support1_fna"),
        exploratory_faa=output_path("exploratory_support1_faa"),
        annotations=output_path("candidate_annotations"),
        script=script("validate_pipeline.py")
    output:
        output_path("validation_checks")
    log:
        output_path("validation_log")
    params:
        python=tool("python"),
        expected_input=config["validation"]["expected_input_sequences"],
        expected_insert_clusters=config["validation"]["expected_insert_clusters"],
        expected_orfs=config["validation"]["expected_predicted_orfs"],
        expected_reporfs=config["validation"]["expected_reporfs"],
        expected_primary=config["validation"]["expected_primary_candidates"]
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        {params.python:q} {input.script:q} \
            --input-fasta {input.input_fasta:q} \
            --normalized-fasta {input.normalized:q} \
            --id-mapping {input.id_mapping:q} \
            --sequence-qc {input.sequence_qc:q} \
            --insert-members {input.insert_members:q} \
            --insert-summary {input.insert_summary:q} \
            --insert-representatives {input.insert_reps:q} \
            --orf-fna {input.orf_fna:q} \
            --orf-faa {input.orf_faa:q} \
            --orf-members {input.orf_members:q} \
            --orf-summary {input.orf_summary:q} \
            --representative-orfs {input.representative_orfs:q} \
            --orf-mapping {input.mapping:q} \
            --missing-mappings {input.missing:q} \
            --support-table {input.support:q} \
            --selection-summary {input.selection_summary:q} \
            --primary-table {input.primary_table:q} \
            --primary-fna {input.primary_fna:q} \
            --primary-faa {input.primary_faa:q} \
            --exploratory-table {input.exploratory_table:q} \
            --exploratory-fna {input.exploratory_fna:q} \
            --exploratory-faa {input.exploratory_faa:q} \
            --annotations {input.annotations:q} \
            --expected-input {params.expected_input} \
            --expected-insert-clusters {params.expected_insert_clusters} \
            --expected-orfs {params.expected_orfs} \
            --expected-reporfs {params.expected_reporfs} \
            --expected-primary {params.expected_primary} \
            --output {output:q} > {log:q} 2>&1
        """
