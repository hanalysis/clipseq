//
// Align to ncrna and primary genomes before indexing and sorting
//

//
// MODULES
//
include { BOWTIE_ALIGN                           } from '../../modules/nf-core/bowtie/align/main.nf'
include { STAR_ALIGN                             } from '../../modules/nf-core/star/align/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_COORD } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRANS   } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANS } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NCRNA  } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_NCRNA } from '../../modules/nf-core/samtools/index/main'

//
// SUBWORKFLOWS
//
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_NCRNA } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_STATS_SAMTOOLS                                       } from '../../subworkflows/nf-core/bam_stats_samtools/main'

workflow RNA_ALIGN {
    take:
    fastq      // channel: [ val(meta), [ fastq ] ]
    bt2_index  // channel: [ val(meta), [ index ] ]
    star_index // channel: [ val(meta), [ index ] ]
    gtf        // channel: [ val(meta), [ gtf ] ]
    fasta      // channel: [ val(meta), [ fasta/fa ]

    main:
    ch_versions = Channel.empty()
    ch_gtf = gtf.first() //for some reason it comes in as a queue channel so need to convert to value channel
    //
    // MODULE: Align reads to ncrna genome
    //
    BOWTIE_ALIGN (
        fastq,
        bt2_index.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_SORT_NCRNA( BOWTIE_ALIGN.out.bam )
    SAMTOOLS_INDEX_NCRNA( SAMTOOLS_SORT_NCRNA.out.bam )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT_NCRNA.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_NCRNA.out.versions)

    //
    // MODULE: Align reads that did not align to the ncrna genome to the primary genome
    //

    STAR_ALIGN (
        BOWTIE_ALIGN.out.fastq,
        star_index,
        ch_gtf,
        false,
        '',
        ''
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //
    // MODULE: Index the coord reads
    //
    SAMTOOLS_INDEX_COORD( STAR_ALIGN.out.bam_sorted )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_COORD.out.versions)

    //
    // CHANNEL: Merge reads and index
    //
    ch_coord_bam_bai = STAR_ALIGN.out.bam_sorted
        .join(SAMTOOLS_INDEX_COORD.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX_COORD.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }

    //
    // MODULE: Sort and index transcript based BAM file
    //
    SAMTOOLS_SORT_TRANS( STAR_ALIGN.out.bam_transcript )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRANS.out.versions)
    SAMTOOLS_INDEX_TRANS( SAMTOOLS_SORT_TRANS.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRANS.out.versions)

    emit:
    ncrna_bam        = SAMTOOLS_SORT_NCRNA.out.bam      // channel: [ val(meta), [ bam ] ]
    ncrna_bai        = SAMTOOLS_INDEX_NCRNA.out.bai      // channel: [ val(meta), [ bai ] ]
    ncrna_log        = BOWTIE_ALIGN.out.log                       // channel: [ val(meta), [ txt ] ]

    genome_log       = STAR_ALIGN.out.log                         // channel: [ val(meta), [ txt ] ]
    genome_log_final = STAR_ALIGN.out.log_final                   // channel: [ val(meta), [ txt ] ]
    genome_bam       = STAR_ALIGN.out.bam_sorted                  // channel: [ val(meta), [ bam ] ]
    genome_bai       = SAMTOOLS_INDEX_COORD.out.bai               // channel: [ val(meta), [ bai ] ]

    transcript_bam   = SAMTOOLS_SORT_TRANS.out.bam                // channel: [ val(meta), [ bam ] ]
    transcript_bai   = SAMTOOLS_INDEX_TRANS.out.bai               // channel: [ val(meta), [ bai ] ]

    versions         = ch_versions                                // channel: [ versions.yml ]
}
