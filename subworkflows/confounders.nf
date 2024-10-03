include { filterBED; thinByLD; mergeBEDS; SampleQCFilter; FlashPCA; ScreePlot } from '../modules/confounders.nf'
include { leave_chr_out } from '../modules/utils.nf'

workflow IIDGenotypes{
    take:
        flashpca_excl_reg
        ld_blocks
        bed_files
        qc_file
        traits

    main:
        bed_files_tuple = bed_files
            .combine(qc_file)
            .combine(ld_blocks)
            .combine(traits)
        filtered_bedfiles = filterBED(bed_files_tuple)
        ld_pruned = thinByLD(filtered_bedfiles.combine(flashpca_excl_reg))
        bedfiles_to_be_merged = ld_pruned.collect()
            .map{it -> ["all_genotypes", it]}
        mergeBEDS(bedfiles_to_be_merged)
        SampleQCFilter(mergeBEDS.out.collect())

    emit:
        SampleQCFilter.out
}
