include { pull_ld; compile_ld_information } from '../modules/ldblocks.nf'

workflow ImportSNPs {
    take:
        bgen_files
    main:
        Channel
            .fromPath(params.INPUT_SNPS, checkIfExists: true)
            .splitCsv(header:true)
            .map {
                row -> tuple(row.RSID, row.CHR, row.POS)
            }
            .combine(bgen_files)
            .set { snps_bgen_ch }

    emit:
        snps_bgen_ch
}

workflow ComputeLD {
    take:
        snps_bgen

    main:
        ld_ch = pull_ld(snps_bgen)
        // filter for just sqlite files and collect
        ld_ch.map{ snp, chr, pos, sqlite_files -> sqlite_files }
            .collect()
            .set { sqlite_ch }

        // compile LD_block information
        csv_ch = Channel.fromPath(params.INPUT_SNPS, checkIfExists: true)
        script_ch = Channel.fromPath("$projectDir/py/convert_sqlite.py")
        qtls_ch = compile_ld_information(sqlite_ch, csv_ch, script_ch)
    emit:
        ld_blocks = qtls_ch.ld_pca
}