process TELESCOPE_ASSIGN {
    tag "$meta_bam.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ea/eaa68e02f66957e51ecba656dfeaa8576bada172780d95592105bc316a01a65c/data':
        'community.wave.seqera.io/library/telescope:1.0.3_fix--d176f12022b914cf' }"

    input:
    tuple val(meta_bam), path(bam)
    tuple val(meta_gtf), path(gtf)

    output:
    tuple val(meta_bam), path("*updated.bam"), emit: updated_bam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*other.bam"), emit: other_bam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*updated.sam"), emit: updated_sam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*other.sam"), emit: other_sam, optional: true // only for --updated_sam
    tuple val(meta_bam), path("*.tsv"), emit: tsv, optional: true // for when there's no alignments
    tuple val(meta_bam), path("*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('telescope'), eval("telescope --version | sed '1!d;s/.* //'"), emit: versions_telescope, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_bam.id}"

    """

    echo -n "" > "$prefix"_telescope.log

    telescope \\
    assign \\
    $bam \\
    $gtf \\
    --project $prefix \\
    $args \\
    > "$prefix"_telescope.log 2>&1

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_bam.id}"

    """
    echo $args

    touch ${prefix}-updated.bam
    touch ${prefix}-other.bam
    touch ${prefix}-updated.sam
    touch ${prefix}-other.sam
    touch ${prefix}-telescope_report.tsv

    """
}
