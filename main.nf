nextflow.enable.dsl = 2

params.data_dir = './humanize'
params.scripts_dir = '/data/Abhijna/scripts'
params.input_file = ''
params.prefix = ''

process IgFoldFull {
    input:
    path fasta_file from file(params.input_file)
    output:
    path "${params.prefix}_igf.pdb"

    script:
    """
    python3 ${params.scripts_dir}/igf_full.py ${fasta_file.baseName}
    """
}

process CleanAntibody {
    input:
    path pdb_file from IgFoldFull.out
    output:
    path "${params.prefix}_cleaned_pdb.pdb"

    script:
    """
    bash ${params.scripts_dir}/clean_ab.sh $pdb_file > ${params.prefix}_cleaned_pdb.pdb
    """
}

process ParseCDR {
    input:
    path cleaned_pdb from CleanAntibody.out
    output:
    path "${params.prefix}_cdr_regions.txt"

    script:
    """
    python3 ${params.scripts_dir}/cdr_parse.py $cleaned_pdb > ${params.prefix}_cdr_regions.txt
    """
}

process ResolveAmbiguities {
    input:
    path cdr_file from ParseCDR.out
    output:
    path "${params.prefix}_resolved_cdr.txt"

    script:
    """
    bash ${params.scripts_dir}/ambig.sh $cdr_file > ${params.prefix}_resolved_cdr.txt
    """
}

process GenerateConfig {
    input:
    path resolved_cdr from ResolveAmbiguities.out
    output:
    path "${params.prefix}_config.txt"

    script:
    """
    bash ${params.scripts_dir}/generate_config.sh $resolved_cdr > ${params.prefix}_config.txt
    """
}

workflow {
    if (!params.input_file.endsWith('.fa')) {
        exit 1, 'Error: Please provide an input file with .fa extension'
    }
    IgFoldFull()
    CleanAntibody()
    ParseCDR()
    ResolveAmbiguities()
    GenerateConfig()
}
