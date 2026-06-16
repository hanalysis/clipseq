process ICOUNTMINI_SUMMARY_MQC{
     tag "$meta.id"
     label 'process_single'

     container 'community.wave.seqera.io/library/r-tidyverse:2.0.0--386074aa810af1c3'

     input:
     path ch_icountmini_summary_qc

     output:
     path "summary_subtype_qc.tsv", emit: subtype
     path "summary_type_qc.tsv", emit: type

     when:
     task.ext.when == null || task.ext.when

     script:
     def args = task.ext.args ?: ''

     """
     RScript $projectDir/modules/local/icountmini_summary_mqc/templates/icount_mini_summary_reformat.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """
}
