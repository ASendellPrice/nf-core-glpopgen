process BEAGLE_SUBSET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(beagle_file), val(all_samples), val(keep_samples)

    output:
    tuple val(meta), path("*.beagle.gz")  , emit: beagle
    tuple val(meta), path("*.samples.txt"), emit: samples
    tuple val("${task.process}"), val('gawk'), eval("gawk --version | sed '1!d;s/GNU Awk //;s/,.*//'"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!keep_samples) {
        error("BEAGLE_SUBSET: no samples to keep for '${meta.id}'.")
    }
    def missing = keep_samples.findAll { sample -> !(sample in all_samples) }
    if (missing) {
        error("BEAGLE_SUBSET: samples ${missing} for '${meta.id}' are not in the Beagle sample list ${all_samples}.")
    }
    if (keep_samples.unique(false).size() != keep_samples.size()) {
        error("BEAGLE_SUBSET: duplicate samples in the list to keep for '${meta.id}': ${keep_samples}.")
    }
    if ("${prefix}.beagle.gz" == beagle_file.name) {
        error("BEAGLE_SUBSET: output name '${prefix}.beagle.gz' is the same as the input. Set a different ext.prefix.")
    }
    // Beagle columns: marker, allele1, allele2, then 3 likelihood columns per sample (IndN = Nth sample in all_samples)
    def indices = keep_samples.collect { sample -> all_samples.indexOf(sample) }
    def names   = keep_samples.collect { sample -> "'${sample}'" }.join(' ')
    // Keep the marker/allele columns (1-3) plus the three likelihood columns of each kept sample.
    // Kept samples are renumbered Ind0..Ind(n-1) in the header, in keep_samples order,
    // and the matching sample IDs are written to <prefix>.samples.txt (one per line, Ind0 first).
    """
    printf '%s\\n' ${names} > ${prefix}.samples.txt

    zcat ${beagle_file} \\
        | gawk -v idx="${indices.join(' ')}" -v total=${all_samples.size()} '
            BEGIN { FS = OFS = "\\t"; n = split(idx, sel, " ") }
            NR == 1 && (NF - 3) / 3 != total {
                printf("Beagle file has %d samples but the sample list has %d\\n", (NF - 3) / 3, total) > "/dev/stderr"
                exit 1
            }
            {
                out = \$1 OFS \$2 OFS \$3
                for (i = 1; i <= n; i++) {
                    col = 3 + 3 * sel[i]
                    if (NR == 1) {
                        name = "Ind" (i - 1)
                        out = out OFS name OFS name OFS name
                    } else {
                        out = out OFS \$(col + 1) OFS \$(col + 2) OFS \$(col + 3)
                    }
                }
                print out
            }' \\
        | gzip -c > ${prefix}.beagle.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.beagle.gz
    touch ${prefix}.samples.txt
    """
}
