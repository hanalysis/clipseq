//
// Take consensus peaks and all crosslinks, return consensus map table
//

//
// MODULES
//
include { BEDTOOLS_MAP as CONSENSUS_MAP } from '../../modules/nf-core/bedtools/map/main'
include { BEDTOOLS_SORT as CONSENSUS_SORT       } from '../../modules/nf-core/bedtools/sort/main'
include { DESEQ2_QC              } from '../../modules/local/deseq2_qc'
//
// SUBWORKFLOWS
//
include { GET_CONSENSUS_COUNTS          } from '../../modules/local/consensus_count_table'


workflow CONSENSUS_PEAK_TABLE {
    take:
    all_crosslinks           // channel: [ val(meta), [ xl.bed ] ]
    consensus_peaks          // channel: [ val(meta), [ peaks.bed ] ]
    genome_fasta
    genome_fai               // channel: [ val(meta), [ index ] ]
    gtf                      // channel: [ val(meta), [ gtf ] ]
    table_name               // val "Clippy_Consensus_Counts.tsv" must end in tsv
    ch_deseq2_pca_header_multiqc
    ch_deseq2_clustering_header_multiqc
    skip_deseq2_qc   

    main:
    ch_versions = Channel.empty()


    CONSENSUS_SORT ( all_crosslinks , [])
    ch_consensus_sorted = CONSENSUS_SORT.out.sorted


    consensus_peaks
        .combine(ch_consensus_sorted)
        .map{ meta1, consensuspeaks, meta2, crosslink -> [meta2, consensuspeaks, crosslink] }
        .set { ch_consensus_map }


    //ch_consensus_map.view { item -> "consensus for bedtools map: $item" }

    ch_consensus_map.view()

    CONSENSUS_MAP (
        ch_consensus_map,
        genome_fai,
    )

    CONSENSUS_MAP.out.mapped
        .collect { it[1] }
        .map { crosslinks -> [[id:"allXL"], crosslinks]}
        .set { ch_mapped_xls }

    GET_CONSENSUS_COUNTS (
        ch_mapped_xls,
        gtf,
        table_name
    )

    //
    // Generate QC plots with DESeq2
    //
    ch_deseq2_qc_pdf           = Channel.empty()
    ch_deseq2_qc_rdata         = Channel.empty()
    ch_deseq2_qc_rds           = Channel.empty()
    ch_deseq2_qc_pca_txt       = Channel.empty()
    ch_deseq2_qc_pca_multiqc   = Channel.empty()
    ch_deseq2_qc_dists_txt     = Channel.empty()
    ch_deseq2_qc_dists_multiqc = Channel.empty()
    ch_deseq2_qc_log           = Channel.empty()
    ch_deseq2_qc_size_factors  = Channel.empty()
    if (!skip_deseq2_qc) {
        DESEQ2_QC (
            GET_CONSENSUS_COUNTS.out.tsv,
            ch_deseq2_pca_header_multiqc,
            ch_deseq2_clustering_header_multiqc
        )
        ch_deseq2_qc_pdf           = DESEQ2_QC.out.pdf
        ch_deseq2_qc_rdata         = DESEQ2_QC.out.rdata
        ch_deseq2_qc_rds           = DESEQ2_QC.out.rds
        ch_deseq2_qc_pca_txt       = DESEQ2_QC.out.pca_txt
        ch_deseq2_qc_pca_multiqc   = DESEQ2_QC.out.pca_multiqc
        ch_deseq2_qc_dists_txt     = DESEQ2_QC.out.dists_txt
        ch_deseq2_qc_dists_multiqc = DESEQ2_QC.out.dists_multiqc
        ch_deseq2_qc_log           = DESEQ2_QC.out.log
        ch_deseq2_qc_size_factors  = DESEQ2_QC.out.size_factors
        ch_versions = ch_versions.mix(DESEQ2_QC.out.versions)
    }

    emit:
    mapped_table            = GET_CONSENSUS_COUNTS.out.tsv    // channel: [ val(meta), [ tsv ] ]


    deseq2_qc_pdf           = ch_deseq2_qc_pdf                  // channel: [ pdf ]
    deseq2_qc_rdata         = ch_deseq2_qc_rdata                // channel: [ rdata ]
    deseq2_qc_rds           = ch_deseq2_qc_rds                  // channel: [ rds ]
    deseq2_qc_pca_txt       = ch_deseq2_qc_pca_txt              // channel: [ txt ]
    deseq2_qc_pca_multiqc   = ch_deseq2_qc_pca_multiqc          // channel: [ txt ]
    deseq2_qc_dists_txt     = ch_deseq2_qc_dists_txt            // channel: [ txt ]
    deseq2_qc_dists_multiqc = ch_deseq2_qc_dists_multiqc        // channel: [ txt ]
    deseq2_qc_log           = ch_deseq2_qc_log                  // channel: [ txt ]
    deseq2_qc_size_factors  = ch_deseq2_qc_size_factors         // channel: [ txt ]
    versions         = ch_versions                   // channel: [ versions.yml ]
}
