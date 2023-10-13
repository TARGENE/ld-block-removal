process pull_ld {
    input:
    tuple val(RSID_LABEL), val(RSID), val(CHR), val(POS)
    path BGEN_FILES
    each PREFIX

    output:
    tuple val(CHR), val(POS), val(RSID), val(RSID_LABEL), path("${RSID}.sqlite")
 
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

    # Use bgenix to create a new .bgen file with the info for just the SNP of interest + index file
    bgenix -g \$BGEN_FILE -incl-rsids ${RSID_LABEL} > "${RSID}.bgen"
    bgenix -g "${RSID}.bgen" -index

    bgenix -g \$BGEN_FILE -incl-range "\$CHR_FORMAT:\$LOWER-\$UPPER" | qctool -g - -filetype bgen \
           -s \$SAMPLE_FILE -compute-ld-with "${RSID}.bgen" \$SAMPLE_FILE -old "sqlite://${RSID}.sqlite:LD"  \
           -min-r2 0.05
    """

}

process compile_ld_information {
    label 'python'
    publishDir("$projectDir/output/", mode: "copy")

    input:
    path sqlite_files
    path csv
    path script

    output:
    path "*.csv"
    path "LD_block_length_histogram.png"

    script:
    """
    python ${script} --input ${csv}
    """

}
