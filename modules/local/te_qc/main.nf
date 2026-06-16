process TE_QC{
    label "process_single"

    container 'community.wave.seqera.io/library/pandas_python:75380848c527810f'

    input:
    path summary_types
    path summary_subtypes

    output:
    path "*telescope_qc.tsv", optional: true, emit: tele_qc
    path "*tetranscripts_qc.tsv", optional:true, emit: tetr_qc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    python ${projectDir}/modules/local/te_qc/templates/telescope_qc.py
    python ${projectDir}/modules/local/te_qc/templates/tetranscripts_qc.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

}
