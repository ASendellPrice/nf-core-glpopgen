process PCANGSD_PCA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pcangsd:1.36.4--py313h5d164f8_1':
        'quay.io/biocontainers/pcangsd:1.36.4--py313h5d164f8_1' }"

    input:
    tuple val(meta), path(beagle_file)

    output:
    tuple val(meta), path("*.cov"), emit: covariance
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('pcangsd'), eval("pcangsd --version"), topic: versions, emit: versions_pcangsd

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pcangsd \\
        --threads ${task.cpus} \\
        --beagle ${beagle_file} \\
        --out ${prefix} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.cov
    touch ${prefix}.log
    """
}
