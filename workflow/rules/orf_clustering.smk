rule create_orf_database:
    input:
        output_path("predicted_orfs_fna")
    output:
        dbtype("orf_db_prefix")
    log:
        output_path("orf_createdb_log")
    params:
        mmseqs=tool("mmseqs"),
        prefix=output_path("orf_db_prefix")
    shell:
        r"""
        mkdir -p "$(dirname {params.prefix:q})" "$(dirname {log:q})"
        {params.mmseqs:q} createdb {input:q} {params.prefix:q} > {log:q} 2>&1
        """


rule cluster_orfs:
    input:
        dbtype("orf_db_prefix")
    output:
        dbtype("orf_cluster_db_prefix")
    threads:
        config["resources"]["mmseqs_threads"]
    log:
        output_path("orf_cluster_log")
    params:
        mmseqs=tool("mmseqs"),
        database=output_path("orf_db_prefix"),
        cluster=output_path("orf_cluster_db_prefix"),
        tmp=output_path("orf_tmp_dir"),
        identity=config["orf_clustering"]["min_sequence_identity"],
        coverage=config["orf_clustering"]["coverage"],
        cov_mode=config["orf_clustering"]["coverage_mode"],
        cluster_mode=config["orf_clustering"]["cluster_mode"]
    shell:
        r"""
        mkdir -p {params.tmp:q} "$(dirname {log:q})"
        {params.mmseqs:q} cluster {params.database:q} {params.cluster:q} {params.tmp:q} \
            --min-seq-id {params.identity} -c {params.coverage} \
            --cov-mode {params.cov_mode} --cluster-mode {params.cluster_mode} \
            --threads {threads} > {log:q} 2>&1
        """


rule write_orf_cluster_membership:
    input:
        database=dbtype("orf_db_prefix"),
        cluster=dbtype("orf_cluster_db_prefix")
    output:
        output_path("orf_cluster_members")
    threads:
        config["resources"]["mmseqs_threads"]
    log:
        output_path("orf_members_log")
    params:
        mmseqs=tool("mmseqs"),
        database=output_path("orf_db_prefix"),
        cluster=output_path("orf_cluster_db_prefix")
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        {params.mmseqs:q} createtsv {params.database:q} {params.database:q} {params.cluster:q} \
            {output:q} --threads {threads} > {log:q} 2>&1
        """


rule summarize_orf_clusters:
    input:
        data=output_path("orf_cluster_members"),
        script=script("summarize_clusters.py")
    output:
        output_path("orf_cluster_summary")
    params:
        python=tool("python"),
    shell:
        "{params.python:q} {input.script:q} --input {input.data:q} --output {output:q}"


rule create_orf_representative_database:
    input:
        database=dbtype("orf_db_prefix"),
        cluster=dbtype("orf_cluster_db_prefix")
    output:
        dbtype("orf_rep_db_prefix")
    params:
        mmseqs=tool("mmseqs"),
        database=output_path("orf_db_prefix"),
        cluster=output_path("orf_cluster_db_prefix"),
        representative=output_path("orf_rep_db_prefix")
    shell:
        "{params.mmseqs:q} createsubdb {params.cluster:q} {params.database:q} {params.representative:q}"


rule write_orf_representatives:
    input:
        dbtype("orf_rep_db_prefix")
    output:
        output_path("representative_orfs_fna")
    log:
        output_path("orf_representatives_log")
    params:
        mmseqs=tool("mmseqs"),
        representative=output_path("orf_rep_db_prefix")
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        {params.mmseqs:q} convert2fasta {params.representative:q} {output:q} > {log:q} 2>&1
        """


rule compute_reporf_support:
    input:
        orfs=output_path("predicted_orfs_fna"),
        orf_clusters=output_path("orf_cluster_members"),
        insert_clusters=output_path("insert_cluster_summary"),
        script=script("compute_orf_support.py")
    output:
        support=output_path("reporf_support"),
        mapping=output_path("orf_to_insert_mapping"),
        missing=output_path("missing_mappings")
    log:
        output_path("support_log")
    params:
        python=tool("python"),
    shell:
        r"""
        mkdir -p "$(dirname {output.support:q})" "$(dirname {log:q})"
        {params.python:q} {input.script:q} \
            --orfs {input.orfs:q} \
            --orf-clusters {input.orf_clusters:q} \
            --insert-clusters {input.insert_clusters:q} \
            --output {output.support:q} \
            --mapping-output {output.mapping:q} \
            --missing-output {output.missing:q} > {log:q} 2>&1
        """
