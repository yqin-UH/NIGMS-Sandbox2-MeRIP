/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome'
include { ALIGN_STAR                       } from '../subworkflows/local/align_star'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { PEAKS_CALL_MACS3_HOMER           } from '../subworkflows/local/peaks_call_macs3_homer'
include { BED_CONSENSUS_MACS3 } from '../subworkflows/local/bed_consensus_macs3'

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { DESEQ2_QC              } from '../modules/local/deseq2_qc'
include { PICARD_MERGESAMFILES   } from '../modules/nf-core/picard/mergesamfiles/main'
include { KHMER_UNIQUEKMERS      } from '../modules/nf-core/khmer/uniquekmers/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_meripseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Header files for MultiQC
ch_peak_count_header         = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header         = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header    = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header        = Channel.value(file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true))
ch_deseq2_clustering_header = Channel.value(file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true))

workflow MERIPSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // SUBWORKFLOW: Prepare genome files
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.star_index,
        params.rsem_index
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    ch_fasta    = PREPARE_GENOME.out.fasta
    ch_gtf      = PREPARE_GENOME.out.gtf

    //
    // SUBWORKFLOW:  Run RNA-seq FASTQ preprocessing subworkflow
    //

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
        ch_samplesheet,
        params.skip_fastqc || params.skip_qc,
        false,
        false,
        params.skip_trimming,
        0,
        10000
    )  
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions.first())

    //
    // MODULE: Run STAR alignment on trimmed reads
    //
    ch_genome_bam          = Channel.empty()
    ch_genome_bam_index    = Channel.empty()
    ch_star_log            = Channel.empty()
    ch_unaligned_sequences = Channel.empty()
    ch_transcriptome_bam   = Channel.empty()

    def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

    ALIGN_STAR (
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads,  // Trimmed FASTQ files from Trim Galore!
        PREPARE_GENOME.out.star_index.map { [ [:], it ] },  // Transform to [meta2, path]
        PREPARE_GENOME.out.gtf.map { [ [:], it ] },        // Transform to [meta3, path]
        params.star_ignore_sjdbgtf ?: false,
        '',
        params.seq_center ?: '',
        is_aws_igenome,
        PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )

    ch_genome_bam          = ALIGN_STAR.out.bam
    ch_genome_bam_index    = ALIGN_STAR.out.bai
    ch_transcriptome_bam   = ALIGN_STAR.out.bam_transcript
    ch_star_log            = ALIGN_STAR.out.log_final
    ch_unaligned_sequences = ALIGN_STAR.out.fastq
    ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})
    ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    
    ch_multiqc_files = ch_multiqc_files
                .mix(ALIGN_STAR.out.stats.collect{it[1]})
                .mix(ALIGN_STAR.out.flagstat.collect{it[1]})
                .mix(ALIGN_STAR.out.idxstats.collect{it[1]})

    //
    // MODULE: Merge resequenced BAM files
    //
    /*
    ch_genome_bam
        .map {
            meta, bam ->
                def meta_clone = meta.clone()
                //meta_clone.remove('read_group')
                meta_clone.id = meta_clone.id - ~/_T\d+$/
                [ meta_clone, bam ]
        }
        .groupTuple(by: [0])
        .map {
            meta, bam ->
                [ meta, bam.flatten() ]
        }
        .set { ch_sort_bam }

    PICARD_MERGESAMFILES (
        ch_sort_bam
    )
    ch_sort_bam = PICARD_MERGESAMFILES.out.bam
        .map { meta, bam -> [meta, file(bam).name.replace('_merge.bam', '.bam')] }
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first())
    */

    ch_genome_bam.join(ch_genome_bam_index, by: [0])
        .set {ch_genome_bam_bai }
    ch_genome_bam_bai
        .map {
            meta, bam, bai ->
                meta.control ? null : [ meta.id, [ bam ] , [ bai ] ]
        }
        .set { ch_control_bam_bai }
    ch_genome_bam_bai
        .map {
            meta, bam, bai ->
                meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null
        }
        .combine(ch_control_bam_bai, by: 0)
        .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
        .set { ch_ip_control_bam_bai }


    //
    // MODULE: Calculute genome size with khmer
    //
    ch_macs_gsize     = Channel.empty()
    ch_macs_gsize     = params.macs_gsize
    if (!params.macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
    }

    // Create channels: [ meta, ip_bam, control_bam ]
    ch_ip_control_bam_bai
        .map {
            meta, bams, bais ->
                [ meta , bams[0], bams[1] ]
        }
        .set { ch_ip_control_bam }

    //
    // SUBWORKFLOW: Call peaks with MACS3, annotate with HOMER and perform downstream QC
    //
    PEAKS_CALL_MACS3_HOMER (
        ch_ip_control_bam,
        ch_fasta,
        ch_gtf,
        ch_macs_gsize,
        "_peaks.annotatePeaks.txt",
        ch_peak_count_header,
        ch_frip_score_header,
        ch_peak_annotation_header,
        params.narrow_peak,
        params.skip_peak_annotation,
        params.skip_peak_qc
    )
    ch_versions = ch_versions.mix(PEAKS_CALL_MACS3_HOMER.out.versions)

    //
    //  Consensus peaks analysis
    
    ch_subreadfeaturecounts_multiqc   = Channel.empty()//
    ch_macs3_consensus_bed_lib   = Channel.empty()
    ch_macs3_consensus_txt_lib   = Channel.empty()
    ch_deseq2_pca_multiqc        = Channel.empty()
    ch_deseq2_clustering_multiqc = Channel.empty()
    if (!params.skip_consensus_peaks) {
        // Create channels: [group, [ ip_bams ] ]
        ch_ip_control_bam
            .map {
                meta, ip_bam, control_bam ->
                   [ meta.group, ip_bam ]
            }
            .groupTuple()
            .set { ch_group_bams }

        BED_CONSENSUS_MACS3 (
            PEAKS_CALL_MACS3_HOMER.out.peaks,
            ch_group_bams,
            ch_fasta,
            ch_gtf,
            ch_deseq2_pca_header,
            ch_deseq2_clustering_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_deseq2_qc
        )
        ch_macs3_consensus_bed_lib       = BED_CONSENSUS_MACS3.out.consensus_bed
        ch_macs3_consensus_txt_lib       = BED_CONSENSUS_MACS3.out.consensus_txt
        ch_subreadfeaturecounts_multiqc  = BED_CONSENSUS_MACS3.out.featurecounts_summary
        ch_deseq2_pca_multiqc            = BED_CONSENSUS_MACS3.out.deseq2_qc_pca_multiqc
        ch_deseq2_clustering_multiqc     = BED_CONSENSUS_MACS3.out.deseq2_qc_dists_multiqc
        ch_versions = ch_versions.mix(BED_CONSENSUS_MACS3.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'meripseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
