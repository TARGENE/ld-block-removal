# ld-block-removal

## Data format
* Input qtls.csv file must be specified in the nextflow.config file and have the columns RSID, CHR and POS for each SNP you would like to examine
* RSID must match the format used in the BGEN files for the cohort-of-interest
* CHR must be just an integer values, without the prefix `chr`
* POS must also be an integer value and must match the corresponding BGEN files and genome assembly chosen

