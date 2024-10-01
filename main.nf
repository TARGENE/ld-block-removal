#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Default parameters
params.TRAITS_DATASET = "You need to provide a Traits dataset."
params.UKB_ENCODING_FILE = "NO_UKB_ENCODING_FILE"
params.MAF_THRESHOLD = 0.01
params.COHORT = "UKB"
params.UKB_CONFIG = "${projectDir}/assets/ukbconfig.yaml"
params.UKB_WITHDRAWAL_LIST = "${projectDir}/assets/NO_WITHDRAWAL_LIST"
params.OUTDIR = "${launchDir}/results"
params.NB_PCS = 50
params.QC_FILE = "${projectDir}/assets/NO_QC_FILE"
params.FLASHPCA_EXCLUSION_REGIONS = "${projectDir}/assets/exclusion_regions_hg19.txt"

// Import subworkflows
include { ImportSNPs; ComputeLD } from './subworkflows/ldblocks.nf'
include { IIDGenotypes; FlashPCA; ScreePlot } from './subworkflows/confounders.nf'
include { ExtractTraits } from './subworkflows/extract_traits.nf'

// Define workflow
workflow {
    // Define Parameters
    ukb_encoding_file = params.UKB_ENCODING_FILE
    ukb_config = Channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

    qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect().toList()

    // Compute LD blocks
    ImportSNPs(bgen_files)

    ComputeLD(ImportSNPs.out)
    
    // Extract Traits
    ExtractTraits(
        traits_dataset,
        ukb_config,
        ukb_withdrawal_list,
        ukb_encoding_file,
    )
    
    // IID Genotypes
    IIDGenotypes(
        flashpca_excl_reg,
        ComputeLD.out.ld_blocks,
        bed_files,
        qc_file,
        ExtractTraits.out,
    )
    
    // PCA
    FlashPCA(IIDGenotypes.out)
    ScreePlot(FlashPCA.out.pve)
    
}

