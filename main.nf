#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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

// Define input channels
// Import a CSV file with all SNPs
workflow import_snps {
    main:
    Channel
        .fromPath(params.INPUT_SNPS, checkIfExists: true)
        .splitCsv(header:true)
        .map {
            row -> tuple(row.CHR, row.POS, row.RSID)
        }
        .set { snps }

    emit:
    snps
}

// Define processes
process removePrefix {
    input:
    tuple val(prefix), path(files)
    each common_prefix

    output:
    tuple val(stripped_prefix), path(files)

    script:
    stripped_prefix = prefix.replaceAll(common_prefix,'')
    "echo Stripped prefix ${common_prefix} from filename"
}

// Define workflow
workflow {
    snps_ch = import_snps()

    bgen_files = Channel.fromFilePairs("$params.BGEN_FILES", size: 3, checkIfExists: true)

    // these files/keys have a matching prefix that is generally associated with the cohort name
    prefix = bgen_files.map({k, v -> k}).collect().map({v -> longest_common_prefix(v)})

    bgen_files_ch = removePrefix(bgen_files, prefix)
    snps_ch.join(bgen_files_ch)
        .set{ input_ch }

    input_ch.view()
}


