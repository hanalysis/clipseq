process LINUX_COMMAND {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta) , path(input)
    path input2
    val copy_input

    output:
    tuple val(meta), path("*.cmd.*"), emit: file
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def ext    = task.ext.ext ?: 'txt'
    def cmd1   = task.ext.cmd1 ?: 'echo "NO-ARGS"'
    def cmd2   = task.ext.cmd2 ? "CMD2=`cat $input2 | ${task.ext.cmd2}`" : ''
    if(copy_input) {
        cmd2 = task.ext.cmd2 ? "CMD2=`cat $input | ${task.ext.cmd2}`" : ''
    }

    """
    $cmd2
    cat $input | $cmd1 > ${prefix}.cmd.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linux: NOVERSION
    END_VERSIONS
    """
}
