//
// Take consensus peaks and all crosslinks, return consensus map table
//

//
// MODULES
//
include { LINUX_COMMAND as CONSENSUS_PEAKS_SORT  } from '../../modules/local/linux_command'
include { LINUX_COMMAND as CROSSLINKS_SORT       } from '../../modules/local/linux_command'
include { BEDTOOLS_MAP  as CONSENSUS_MAP         } from '../../modules/nf-core/bedtools/map/main'
//
// SUBWORKFLOWS
//
include { GET_CONSENSUS_COUNTS          } from '../../modules/local/consensus_count_table'


workflow CONSENSUS_PEAK_TABLE {
    take:
    all_crosslinks           // channel: [ val(meta), [ xl.bed ] ]
    consensus_peaks          // channel: [ val(meta), [ peaks.bed ] ]
    genome_fai               // channel: [ val(meta), [ index ] ]
    gtf                      // channel: [ val(meta), [ gtf ] ]
    table_name               // val "Clippy_Consensus_Counts.tsv" must end in tsv

    main:
    ch_versions = Channel.empty()

    // Sort consensus peaks according to genome file order
    CONSENSUS_PEAKS_SORT (
        consensus_peaks,
        [],
        false
    )

    // Sort crosslinks according to genome file order
    CROSSLINKS_SORT (
        all_crosslinks,
        [],
        false
    )

    // Combine sorted peaks with sorted crosslinks
    CONSENSUS_PEAKS_SORT.out.file
        .combine(CROSSLINKS_SORT.out.file)
        .map{ meta1, consensuspeaks, meta2, crosslink -> [meta2, consensuspeaks, crosslink] }
        .set { ch_consensus_map }

    CONSENSUS_MAP (
        ch_consensus_map,
        genome_fai
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

    emit:
    mapped_table     = GET_CONSENSUS_COUNTS.out.tsv    // channel: [ val(meta), [ tsv ] ]
    versions         = ch_versions                   // channel: [ versions.yml ]
}
