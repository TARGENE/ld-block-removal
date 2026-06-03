include { pull_ld; compile_ld_information } from '../modules/ldblocks.nf'

workflow ImportSNPs {
    take:
        bgen_files
    main:
        channel
            .fromPath(params.INPUT_SNPS, checkIfExists: true)
            .splitCsv(header:true)
            .map {
                row -> tuple(row.CHR, row.RSID, row.POS)
            }
            .set {snps_ch}
        
        bgen_files
            .map { prefix, files ->
                def match = (prefix =~ /\.chr(\d+|X|Y)/)
                def chr = match.find() ? match.group(1) : null
                if (chr == null) {
                    println "Warning: Could not extract chromosome number from '${prefix}'"
                    chr = "null"
                }
                return [chr, prefix, files]
            }
            .set { bgen_ch }

        // add proper chromosome BGEN files into SNP channel
        snps_bgen_ch = snps_ch.combine(bgen_ch)
            .filter { it -> it[0] == it[3] }
            .map { chr, snp, pos, _chr2, prefix, files -> [snp, chr, pos, prefix, files]}
        
    emit:
        snps_bgen_ch
}

workflow ComputeLD {
    take:
        snps_bgen

    main:
        ld_ch = pull_ld(snps_bgen)

        // filter for just sqlite files and collect
        ld_ch.map{ _snp, _chr, _pos, sqlite_files -> sqlite_files }
            .collect()
            .set { sqlite_ch }

        // compile LD_block information
        csv_ch = channel.fromPath(params.INPUT_SNPS, checkIfExists: true)
        script_ch = channel.fromPath("$projectDir/py/convert_sqlite.py")
        qtls_ch = compile_ld_information(sqlite_ch, csv_ch, script_ch)
    emit:
        ld_blocks = qtls_ch.ld_blocks
}
