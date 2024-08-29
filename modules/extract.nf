process pull_ld {
    label 'qctool_image'

    input:
    tuple val(RSID), val(CHR), val(POS), val(PREFIX), path(BGEN_FILES)

    output:
    tuple val(RSID), val(CHR), val(POS), path("*.sqlite")
 
    script:
    """
    if [[ "${params.COHORT}" == "UKBB" ]]; then 
        CHR_FORMAT=\$(echo ${CHR} | xargs printf "%02d" )
    elif [[ "${params.COHORT}" == "GENOMICC" ]]; then 
        CHR_FORMAT=\$(echo chr${CHR} )
    fi

    BGEN_FILE="${PREFIX}${CHR}.bgen"
    SAMPLE_FILE="${PREFIX}${CHR}.sample"

    # Define lower and upper bounds to scan for LD
    LOWER=\$((${POS} - 10000000))
    UPPER=\$((${POS} + 10000000))

    # Account for the possibility that the scan range becomes negative.
    if [ \$LOWER -lt 0 ]
    then  
      LOWER="0"
    fi

    if [ \$UPPER -lt 0 ]
    then  
      UPPER="0"
    fi

    OUTNAME=\$( echo ${RSID} | sed 's/:/_/g' )

    # Use bgenix to create a new .bgen file with the info for just the SNP of interest + index file
    bgenix -g \$BGEN_FILE -incl-rsids ${RSID} > "\$OUTNAME.bgen"
    bgenix -g "\$OUTNAME.bgen" -index

    bgenix -g \$BGEN_FILE -incl-range "\$CHR_FORMAT:\$LOWER-\$UPPER" | qctool -g - -filetype bgen \
           -s \$SAMPLE_FILE -compute-ld-with "\$OUTNAME.bgen" \$SAMPLE_FILE -old "sqlite://\$OUTNAME.sqlite:LD"  \
           -min-r2 0.05
    """

}

process compile_ld_information {
    label 'bgen_python_image'
    publishDir("$params.OUTDIR/ld_blocks/", mode: "copy", pattern: "*LD_block*")

    input:
    path sqlite_files
    path csv
    path script

    output:
    path "LD_for_PCA.csv", emit: ld_pca
    path "*_with_LD_blocks.csv"
    path "LD_block_length_histogram.png"

    script:
    """
    python ${script} --input ${csv}
    """

}

