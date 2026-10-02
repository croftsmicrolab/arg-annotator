rule predict_orfs:
    input:
        output_path("insert_representatives")
    output:
        fna=output_path("predicted_orfs_fna"),
        faa=output_path("predicted_proteins_faa"),
        gff=output_path("predicted_orfs_gff")
    log:
        stdout=output_path("prodigal_stdout_log"),
        stderr=output_path("prodigal_stderr_log")
    params:
        prodigal=tool("prodigal"),
        mode=config["orf_prediction"]["mode"],
        output_format=config["orf_prediction"]["output_format"]
    shell:
        r"""
        mkdir -p "$(dirname {output.fna:q})" "$(dirname {log.stdout:q})"
        {params.prodigal:q} \
            -i {input:q} \
            -d {output.fna:q} \
            -a {output.faa:q} \
            -o {output.gff:q} \
            -f {params.output_format:q} \
            -p {params.mode:q} > {log.stdout:q} 2> {log.stderr:q}
        """
