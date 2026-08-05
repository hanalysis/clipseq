process GET_INIT_ALIGNED_CROSSLINKS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*init_xl_coord.bam")     , emit: bam
    path  "versions.yml"                            ,emit: versions


    script:
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    bedtools bamtobed -i $bam > dedup.bed
    bedtools shift -m 1 -p -1 -i dedup.bed -g $fai > shifted.bed
    bedtools bedtobam -g $fai > ${prefix}_init_xl_coord.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linux: NOVERSION
        bedtools: `bedtools --version | head -n 1`
    END_VERSIONS
    """
}
