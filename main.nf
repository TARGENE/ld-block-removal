#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Default parameters
params.OUTDIR = "${launchDir}/results"
params.DECRYPTED_DATASET = "NO_FILE"
params.COHORT = "UKBB"
params.TRAITS_CONFIG = "NO_UKB_TRAIT_CONFIG"
params.WITHDRAWAL_LIST = 'NO_WITHDRAWAL_LIST'
params.QC_FILE = "NO_QC_FILE"
params.MAF_THRESHOLD = 0.01

// Define functions
def longest_common_prefix(List<String> strs) {    
    def longest_common = null    
    for (s in strs) {    
        if (longest_common == null) {    
            longest_common = s    
        } else {    
            for (i=0; i<=longest_common.length() - 1; i++) {    
                if (longest_common.charAt(i) != s.charAt(i)) {    
                    longest_common = longest_common.substring(0, i)    
                    break    
                }    
            }    
        }    
    }    
    return longest_common    
}    

// Import processes
include { pull_ld; compile_ld_information } from './modules/extract.nf'
include { IIDGenotypes } from './modules/genotypes.nf'
include { FlashPCA; AdaptFlashPCA; ScreePlot } from './modules/confounders.nf'

workflow importSNPs {
    main:
        Channel
            .fromPath(params.INPUT_SNPS, checkIfExists: true)
            .splitCsv(header:true)
            .map {
                row -> tuple(row.RSID_LABEL, row.RSID, row.CHR, row.POS)
            }
            .set { snps_ch }

    emit:
        snps_ch
}

workflow configureBGEN {
    main:
        bgen_files = Channel.fromFilePairs("$params.BGEN_FILES", size: 3, checkIfExists: true)
        // these files/keys have a matching prefix that is generally associated with the cohort name
        bgen_prefix = bgen_files.map({k, v -> k}).collect().map({v -> longest_common_prefix(v)})

        bgen_files.map{ key, files -> files }
            .collect()
            .set { bgen_files_ch }

    emit:
        files = bgen_files_ch
        prefix = bgen_prefix
}

workflow computeLD {
    take:
        snps
        bgen_files
        bgen_prefix

    main:
        ld_ch = pull_ld(snps, bgen_files, bgen_prefix)
        // filter for just sqlite files and collect
        ld_ch.map{ chr, pos, snp, snp_label, sqlite_files -> sqlite_files }
            .collect()
            .set { sqlite_ch }

        // compile LD_block information
        csv_ch = Channel.fromPath(params.INPUT_SNPS, checkIfExists: true)
        script_ch = Channel.fromPath("$projectDir/py/convert_sqlite.py")
        qtls_ch = compile_ld_information(sqlite_ch, csv_ch, script_ch)
    emit:
        ld_blocks = qtls_ch.ld_pca
}

workflow extractTraits {
    traits_config = Channel.value(file("$params.TRAITS_CONFIG"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))
    if (params.DECRYPTED_DATASET == "NO_FILE") {
        encrypted_dataset = Channel.value(file("$params.ENCRYPTED_DATASET"))
        encoding_file = Channel.value(file("$params.ENCODING_FILE"))
        UKBFieldsList(traits_config)
        decrypted_dataset = UKBConv(UKBFieldsList.out, encrypted_dataset, encoding_file)
    }
    else {
        decrypted_dataset = Channel.value(file("$params.DECRYPTED_DATASET"))
    }

    TraitsFromUKB(decrypted_dataset, traits_config, withdrawal_list)

    emit:
        TraitsFromUKB.out
}

workflow generateIIDGenotypes {
    take:
        traits
        ld_blocks

    main:
        qc_file = Channel.value(file("$params.QC_FILE"))
        flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
        bed_files_ch = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

        IIDGenotypes(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file, traits)

    emit:
        IIDGenotypes.out
}

workflow geneticConfounders {
    take:
        iid_genotypes

    main:
        FlashPCA(iid_genotypes)
        AdaptFlashPCA(FlashPCA.out.pcs)
        ScreePlot(FlashPCA.out.pve)
}

// Define workflow
workflow {
    // compute regions of LD around each VOI
    snps_ch = importSNPs()

    bgen_ch = configureBGEN()
   
    ld_ch = computeLD(snps_ch, bgen_ch.files, bgen_ch.prefix)

    // run PCA on entire cohort, excluding regions of LD
    if (params.COHORT == "UKBB") {
        extractTraits()
        phenoInput = extractTraits.out
    } else {
        phenoInput = Channel.fromPath("$params.DECRYPTED_DATASET", checkIfExists: true)
    }

    // Generate IID Genotypes
    generateIIDGenotypes(phenoInput, ld_ch.ld_blocks)

    // Genetic confounders up to NB_PCS
    geneticConfounders(generateIIDGenotypes.out)
}


