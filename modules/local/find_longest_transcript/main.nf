process FIND_LONGEST_TRANSCRIPT {
    tag "$gtf"
    label "process_single"

    conda "bioconda::pyranges=0.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1':
        'quay.io/biocontainers/pyranges:0.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.txt")                 ,emit: longest_transcript
    tuple val(meta), path("*.fai")                 ,emit: longest_transcript_fai
    tuple val(meta), path("*.gtf")                 ,emit: longest_transcript_gtf
    tuple val(meta), path("*filtered.gtf")         ,emit: filtered_gtf
    path  "versions.yml"                           ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    output       = task.ext.output ?: "longest_transcript"
    template 'find_longest_transcript.py'
}
