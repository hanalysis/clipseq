process MULTIMAP_CLASS_BINNING {
    label "process_low"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96dde1efad90c922a0198cae64c642be95605c23cf1e53e3c35491817bf6c48b/data':
    'community.wave.seqera.io/library/bedtools_pybedtools_pysam_matplotlib_pruned:79786472f5bd377f' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta_gtf), path(gtf)

    output:
    path "discarded_ambiguous_reads.tsv", emit: discarded_reads
    path "repeat_ncRNA_read_counts.tsv", emit: counted_reads
    tuple val(meta), path("*collapsed.bed"), emit: collapsed_bed
    tuple val(meta), path("*collapsed.bam"), emit: collapsed_bam
    tuple val(meta), path("*collapsed.bam.bai"), emit: collapsed_bai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    python ${projectDir}/modules/local/multimap_class_binning/templates/multimap_class_binning.py \
    --bamdir . \
    --gtf ${gtf} \
    -i ${bam} \
    -o . \
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pandas: \$(python -c "import pybedtools; print(pybedtools.__version__)")
        pandas: \$(python -c "import pysam; print(pysam.__version__)")
        pandas: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pandas: \$(python -c "import pybedtools; print(pybedtools.__version__)")
        pandas: \$(python -c "import pysam; print(pysam.__version__)")
        pandas: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

}
