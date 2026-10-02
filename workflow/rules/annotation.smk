rule diamond_swissprot:
    input:
        query=output_path("primary_faa"),
        database=[database_path("swissprot_diamond")] if SWISS_ENABLED else []
    output:
        output_path("diamond_hits")
    threads:
        config["resources"]["diamond_threads"]
    log:
        output_path("diamond_log")
    params:
        diamond=tool("diamond") if SWISS_ENABLED else "diamond",
        enabled=int(SWISS_ENABLED),
        database=database_path("swissprot_diamond"),
        mode="" if config["annotation"]["swissprot"]["mode"] == "fast" else "--" + config["annotation"]["swissprot"]["mode"],
        evalue=config["annotation"]["swissprot"]["evalue"],
        targets=config["annotation"]["swissprot"]["max_target_sequences"],
        fields=" ".join(config["annotation"]["swissprot"]["output_fields"])
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        if [ {params.enabled} -eq 0 ] || [ ! -s {input.query:q} ]; then
            : > {output:q}
            echo "Swiss-Prot skipped: disabled or no candidates" > {log:q}
        else
        {params.diamond:q} blastp \
            --db {params.database:q} \
            --query {input.query:q} \
            --out {output:q} \
            --outfmt 6 {params.fields} \
            --max-target-seqs {params.targets} \
            --evalue {params.evalue} \
            {params.mode} \
            --threads {threads} > {log:q} 2>&1
        fi
        """


rule rgi_card_protein:
    input:
        query=output_path("primary_rgi_faa"),
        card=[str(RGI_DATABASE / "card.json")] if RGI_ENABLED else [],
        script=script("run_rgi.py")
    output:
        txt=output_path("rgi_txt"),
        json=output_path("rgi_json")
    threads:
        config["resources"]["rgi_threads"]
    log:
        stdout=output_path("rgi_stdout_log"),
        stderr=output_path("rgi_stderr_log")
    params:
        python=tool("python"),
        rgi=tool("rgi") if RGI_ENABLED else "rgi",
        enabled=int(RGI_ENABLED),
        database=str(RGI_DATABASE),
        loose=int(config["annotation"]["rgi"]["include_loose"]),
        low_quality=int(config["annotation"]["rgi"]["low_quality"]),
        clean=int(config["annotation"]["rgi"]["clean"])
    shell:
        r"""
        mkdir -p "$(dirname {log.stdout:q})"
        {params.python:q} {input.script:q} --query {input.query:q} \
            --database {params.database:q} --txt {output.txt:q} --json {output.json:q} \
            --rgi {params.rgi:q} --enabled {params.enabled} --threads {threads} \
            --include-loose {params.loose} --low-quality {params.low_quality} --clean {params.clean} \
            > {log.stdout:q} 2> {log.stderr:q}
        """


rule merge_candidate_annotations:
    input:
        table=output_path("primary_table"),
        faa=output_path("primary_faa"),
        mapping=output_path("orf_to_insert_mapping"),
        diamond=output_path("diamond_hits"),
        rgi=output_path("rgi_txt"),
        script=script("merge_candidate_annotations.py")
    output:
        annotations=output_path("candidate_annotations"),
        summary=output_path("annotation_summary")
    log:
        output_path("merge_log")
    params:
        python=tool("python"),
        rgi_status="searched" if RGI_ENABLED else "disabled",
        swiss_status="searched" if SWISS_ENABLED else "disabled",
        medium_identity=config["annotation"]["evidence_thresholds"]["medium_identity"],
        medium_coverage=config["annotation"]["evidence_thresholds"]["medium_query_coverage"],
        low_identity=config["annotation"]["evidence_thresholds"]["low_identity"],
        low_coverage=config["annotation"]["evidence_thresholds"]["low_query_coverage"]
    shell:
        r"""
        mkdir -p "$(dirname {output.annotations:q})" "$(dirname {log:q})"
        {params.python:q} {input.script:q} \
            --candidate-table {input.table:q} \
            --candidate-faa {input.faa:q} \
            --orf-mapping {input.mapping:q} \
            --diamond {input.diamond:q} \
            --rgi {input.rgi:q} \
            --rgi-status {params.rgi_status:q} \
            --swissprot-status {params.swiss_status:q} \
            --medium-identity {params.medium_identity} \
            --medium-query-coverage {params.medium_coverage} \
            --low-identity {params.low_identity} \
            --low-query-coverage {params.low_coverage} \
            --output {output.annotations:q} \
            --summary {output.summary:q} > {log:q} 2>&1
        """
