rule create_insert_database:
    input:
        output_path("normalized_fasta")
    output:
        dbtype("insert_db_prefix")
    log:
        output_path("insert_createdb_log")
    params:
        mmseqs=tool("mmseqs"),
        prefix=output_path("insert_db_prefix")
    shell:
        r"""
        mkdir -p "$(dirname {params.prefix:q})" "$(dirname {log:q})"
        {params.mmseqs:q} createdb {input:q} {params.prefix:q} > {log:q} 2>&1
        """


rule cluster_inserts:
    input:
        dbtype("insert_db_prefix")
    output:
        dbtype("insert_cluster_db_prefix")
    threads:
        config["resources"]["mmseqs_threads"]
    log:
        output_path("insert_cluster_log")
    params:
        mmseqs=tool("mmseqs"),
        database=output_path("insert_db_prefix"),
        cluster=output_path("insert_cluster_db_prefix"),
        tmp=output_path("insert_tmp_dir"),
        identity=config["insert_clustering"]["min_sequence_identity"],
        coverage=config["insert_clustering"]["coverage"],
        cov_mode=config["insert_clustering"]["coverage_mode"],
        cluster_mode=config["insert_clustering"]["cluster_mode"]
    shell:
        r"""
        mkdir -p {params.tmp:q} "$(dirname {log:q})"
        {params.mmseqs:q} cluster {params.database:q} {params.cluster:q} {params.tmp:q} \
            --min-seq-id {params.identity} -c {params.coverage} \
            --cov-mode {params.cov_mode} --cluster-mode {params.cluster_mode} \
            --threads {threads} > {log:q} 2>&1
        """


rule write_insert_cluster_membership:
    input:
        database=dbtype("insert_db_prefix"),
        cluster=dbtype("insert_cluster_db_prefix")
    output:
        output_path("insert_cluster_members")
    threads:
        config["resources"]["mmseqs_threads"]
    log:
        output_path("insert_members_log")
    params:
        mmseqs=tool("mmseqs"),
        database=output_path("insert_db_prefix"),
        cluster=output_path("insert_cluster_db_prefix")
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        {params.mmseqs:q} createtsv {params.database:q} {params.database:q} {params.cluster:q} \
            {output:q} --threads {threads} > {log:q} 2>&1
        """


rule summarize_insert_clusters:
    input:
        data=output_path("insert_cluster_members"),
        script=script("summarize_clusters.py")
    output:
        output_path("insert_cluster_summary")
    params:
        python=tool("python"),
    shell:
        "{params.python:q} {input.script:q} --input {input.data:q} --output {output:q}"


rule create_insert_representative_database:
    input:
        database=dbtype("insert_db_prefix"),
        cluster=dbtype("insert_cluster_db_prefix")
    output:
        dbtype("insert_rep_db_prefix")
    params:
        mmseqs=tool("mmseqs"),
        database=output_path("insert_db_prefix"),
        cluster=output_path("insert_cluster_db_prefix"),
        representative=output_path("insert_rep_db_prefix")
    shell:
        "{params.mmseqs:q} createsubdb {params.cluster:q} {params.database:q} {params.representative:q}"


rule write_insert_representatives:
    input:
        dbtype("insert_rep_db_prefix")
    output:
        output_path("insert_representatives")
    log:
        output_path("insert_representatives_log")
    params:
        mmseqs=tool("mmseqs"),
        representative=output_path("insert_rep_db_prefix")
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        {params.mmseqs:q} convert2fasta {params.representative:q} {output:q} > {log:q} 2>&1
        """
