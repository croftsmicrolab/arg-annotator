rule prepare_input_fasta:
    input:
        data=input_path("fasta"),
        script=script("prepare_fasta.py")
    output:
        fasta=output_path("normalized_fasta"),
        mapping=output_path("id_mapping"),
        qc=output_path("sequence_qc")
    log:
        output_path("prepare_log")
    params:
        python=tool("python"),
        symbols=config["qc"]["allowed_nucleotide_symbols"],
        width=config["qc"]["fasta_line_width"],
        invalid_flag="--fail-on-invalid-bases" if config["qc"]["fail_on_invalid_bases"] else "",
        empty_flag="--fail-on-empty-sequences" if config["qc"]["fail_on_empty_sequences"] else ""
    shell:
        r"""
        mkdir -p "$(dirname {output.fasta:q})" "$(dirname {log:q})"
        {params.python:q} {input.script:q} \
            --input {input.data:q} \
            --output {output.fasta:q} \
            --mapping {output.mapping:q} \
            --qc {output.qc:q} \
            --allowed-symbols {params.symbols:q} \
            --line-width {params.width} \
            {params.invalid_flag} {params.empty_flag} > {log:q} 2>&1
        """
