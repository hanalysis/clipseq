include { SAMTOOLS_SIMPLE_VIEW as FILTER_TRANSCRIPTS       } from '../../modules/local/samtools_simple_view'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as TRANS_DEDUP    } from './bam_dedup_samtools_umicollapse'
include { CALC_CROSSLINKS as CALC_TRANSCRIPT_CROSSLINKS    } from './calc_crosslinks'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILT_TRANSCRIPT   } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILT_TRANSCRIPT } from '../../modules/nf-core/samtools/index/main'
include { CLIPPY as CLIPPY_TRANSCRIPTOME                   } from "../../modules/nf-core/clippy/main"
include { PARACLU as PARACLU_TRANSCRIPTOME                 } from "../../modules/nf-core/paraclu/main"



workflow TRANSCRIPTOME_PROCESSING {
    take:
    ch_transcript_bam          // channel: [ val(meta), [ bam ] ]
    ch_transcript_bai          // channel: [ val(meta), [ bai ] ]
    ch_longest_transcript      // channel: []
    ch_longest_transcript_gtf  // channel: []
    ch_longest_transcript_fai  // channel: []
    ch_fasta                   // channel: []
    callers
    ch_paraclu_mincluster

    main:
    ch_versions = Channel.empty()
    if(params.source == "fastq" & params.run_filtering) {
        //
        // CHANNEL: Combine bam and bai files on id
        //
        ch_transcript_bam_bai = ch_transcript_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }
        //ch_transcript_bam_bai | view

        //
        // MODULE: Filter transcriptome bam on longest transcripts
        //
        FILTER_TRANSCRIPTS (
            ch_transcript_bam_bai,
            [],
            ch_longest_transcript
        )
        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS.out.versions)

        //
        // Sort, index filtered bam
        //
        SAMTOOLS_SORT_FILT_TRANSCRIPT ( FILTER_TRANSCRIPTS.out.bam )
        ch_versions       = ch_versions.mix(SAMTOOLS_SORT_FILT_TRANSCRIPT.out.versions)
        ch_transcript_bam = SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam

        SAMTOOLS_INDEX_FILT_TRANSCRIPT ( SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam )
        ch_versions       = ch_versions.mix(SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.versions)
        ch_transcript_bai = SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.bai
    }

    // UMI DEDUPLICATION
    ch_trans_umi_log  = Channel.empty()
    if(params.source == "fastq" & params.run_dedup) {
        ch_transcript_bam_bai = ch_transcript_bam
        .map { row -> [row[0].id, row ].flatten()}
        .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4]] }

        //
        // SUBWORKFLOW: Run umi deduplication on transcript-level alignments
        //
        TRANS_DEDUP (
            ch_transcript_bam_bai,
            ch_fasta
        )
        ch_versions        = ch_versions.mix(TRANS_DEDUP.out.versions)
        ch_transcript_bam  = TRANS_DEDUP.out.bam
        ch_transcript_bai  = TRANS_DEDUP.out.bai
        ch_trans_umi_log   = TRANS_DEDUP.out.umi_log
    }

    ch_trans_crosslink_bed            = Channel.empty()
    ch_trans_crosslink_coverage       = Channel.empty()
    ch_trans_crosslink_coverage_norm  = Channel.empty()

    if(params.run_crosslinking) {
        //
        // SUBWORKFLOW: Run crosslink calculation for transcripts
        //
        CALC_TRANSCRIPT_CROSSLINKS (
            ch_transcript_bam,
            ch_longest_transcript_fai.map{ [[id:it.baseName], it] }
        )
        ch_versions                      = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS.out.versions)
        ch_trans_crosslink_bed           = CALC_TRANSCRIPT_CROSSLINKS.out.bed
        ch_trans_crosslink_coverage      = CALC_TRANSCRIPT_CROSSLINKS.out.coverage
        ch_trans_crosslink_coverage_norm = CALC_TRANSCRIPT_CROSSLINKS.out.coverage_norm
    }


    ch_clippy_transcriptome_peaks   = Channel.empty()
    ch_paraclu_transcriptome_peaks  = Channel.empty()
    if(params.run_peakcalling) {
        if('clippy' in callers) {
            CLIPPY_TRANSCRIPTOME (
                ch_trans_crosslink_bed,
                ch_longest_transcript_gtf,
                ch_longest_transcript_fai
            )
            
            ch_clippy_transcriptome_peaks    = CLIPPY_TRANSCRIPTOME.out.peaks
            ch_versions                      = ch_versions.mix(CLIPPY_TRANSCRIPTOME.out.versions)
        }

        if('paraclu' in callers) {
            PARACLU_TRANSCRIPTOME (
                ch_trans_crosslink_bed,
                ch_paraclu_mincluster
            )
            ch_paraclu_transcriptome_peaks          = PARACLU_TRANSCRIPTOME.out.bed
            ch_versions                             = ch_versions.mix(PARACLU_TRANSCRIPTOME.out.versions)
        }
    }


    emit:
    transcript_dedupe_bam   = ch_transcript_bam                  // channel: [ val(meta), [ bam ] ]
    transcript_dedupe_bai   = ch_transcript_bai                  // channel: [ val(meta), [ bai ] ]

    crosslink_bed           = ch_trans_crosslink_bed             // channel: [ val(meta), [ bed ] ] 
    crosslink_coverage      = ch_trans_crosslink_coverage        // channel: [ val(meta), [ bedGraph ] ] 
    crosslink_coverage_norm = ch_trans_crosslink_coverage_norm   // channel: [ val(meta), [ bedGraph ] ] 

    clippy_peaks            = ch_clippy_transcriptome_peaks      // channel: [ val(meta), [ bed ] ]
    paraclu_peaks           = ch_paraclu_transcriptome_peaks     // channel: [ val(meta), [ bed ] ]
 
    versions                = ch_versions                        // channel: [ versions.yml ]
}