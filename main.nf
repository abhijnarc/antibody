nextflow.enable.dsl = 2

params.data_dir = './humanize'
params.scripts_dir = '/data/Abhijna/scripts'
params.input_file = ''
params.prefix = ''

process IgFoldFull {
    input:
    path fasta_file from file(params.input_file)
    output:
    path "/data/Abhijna/humanize/${params.prefix}_igf.pdb"

    script:
    """
    python3 ${params.scripts_dir}/igf_full.py ${fasta_file.baseName}
    """
}

process CleanAntibody {
    input:
    path pdb_file from IgFoldFull.out
    output:
    path "/data/Abhijna/humanize/${params.prefix}_igf_clean.pdb"

    script:
    """
    bash ${params.scripts_dir}/clean_ab.sh $pdb_file > ${params.prefix}_igf_clean.pdb
    """
}

process ParseCDR {
    input:
    path cleaned_pdb from CleanAntibody.out
    output:
    path "/data/Abhijna/humanize/${params.prefix}_coord.csv"

    script:
    """
    python3 ${params.scripts_dir}/cdr_parse.py $cleaned_pdb > ${params.prefix}_coord.csv
    """
}

process AntibodyRestraint {
    input:
    path cdr_coord from AntibodyRestraint.out
    output:
    path "/data/Abhijna/humanize/${params.prefix}_res.txt"

    script:
    """
    python3 ${params.scripts_dir}/ab_res.py $cdr_coord $fasta_file $ > ${params.prefix}_.res.txt
    """
}

process GenerateAmbig {
    input:
    path ab_res from AntibodyRestraint.out
    output:
    path "${params.prefix}"

    script:
    """
    bash ${params.scripts_dir}/ambig.sh $ab_res > ambig${params.prefix}.tbl
    """
}
process GenerateUnambig {
    input:
    path cleaned_pdb from CleanAntibody.out
    output:
    path "${params.prefix}"

    script:
    """
    bash ${params.scripts_dir}/ambig.sh $cleaned_pdb > unambig${params.prefix}.tbl
    """
}
process GenerateConfig {
    input:
    path prefix from file(params.input_file)
    output:
    path "${params.prefix}.cfg"

    script:
    """
    bash ${params.scripts_dir}/generate_config.sh $prefix > ${params.prefix}.cfg
    """
}

workflow {
    if (!params.input_file.endsWith('.fa')) {
        exit 1, 'Error: Please provide an input file with .fa extension'
    }
    IgFoldFull()
    CleanAntibody()
    ParseCDR()
    AntibodyRestraint()
    GenerateAmbig()
    GenerateUnambig()
    GenerateConfig()
}
