# ld-block-removal

## Software requirements
To run this pipeline, you must have nextflow and singularity installed. This can be done using the conda environment `env.yaml`. First conda must be installed, which can be done following the instructions here: https://docs.anaconda.com/miniconda/miniconda-install/. Once installed, the environment can be created as follows:

```
conda env create --name env.yaml
```

Once installed, activate the environment in order to run the pipeline by:

```
conda activate nextflow
```

## Setup

### Configuration
Input data and parameters are specified in a run-specific configuration file created by the user. An example is found in `test/test.config`, and can include the following parameters:
* `COHORT` [ default = "UKB" ] : Genomic cohort you are working with. Current supported cohorts are UKB and GENOMICC. These are configured to format the chromosome values as they are represented in the BGEN files for a given cohort. If UKB or GENOMICC are not specified, the default will be tot ake the integer value for each chromosome, with no prefix.
* `INPUT_SNPS` [ required ] : Path to a CSV file, where each line represents a SNP-of-interest. Columns must be RSID, CHR, POS. RSID must be in the format as it appears in the BGEN files. 
* `BGEN_FILES` [ required ] : Path to BGEN files for the cohort-of-interest. Must be split by chromosome and specified as so: "/path/to/bgen/files/test_cohort_chr*.{bgen,sample,bgen.bgi}".
* `BED_FILES` [ required ] : Path to BED files to use for PCA after LD blocks have been removed. Must be split by chromosome and specified as so: "/path/to/bed/files/test_cohort_chr*.{bed,bim,fam}".
* `TRAITS_DATASET` [ required ] : see https://targene.github.io/targene-pipeline/stable/all_workflows_parameters/.
* `RUN_TYPE` [ default = "DISCOVERY" ] : Any non-test run.
* `NB_PCS` [ default = 50 ] : Number of PCs to compute usign FlashPCA2 after removing LD blocks. 
* `FLASHPCA_EXCLUSION_REGIONS` [ default = assets/exclusion_regions_hg19.txt ] : Regions to exclude from BED files before PCA. hg38 also available in assets/exclusion_regions_hg38.txt.
* `UKB_ENCODING_FILE` [ default = "NO_UKB_ENCODING_FILE" ] : see https://targene.github.io/targene-pipeline/stable/all_workflows_parameters/.
* `UKB_WITHDRAWAL_LIST` [ default = "assets/NO_WITHDRAWAL_LIST" ] : see https://targene.github.io/targene-pipeline/stable/all_workflows_parameters/.
* `UKB_CONFIG` [ default = "assets/ukbconfig.yaml" ] : see https://targene.github.io/targene-pipeline/stable/all_workflows_parameters/.
* `QC_FILE` [ default = "assets/NO_QC_FILE" ] : see https://targene.github.io/targene-pipeline/stable/all_workflows_parameters/.
* `MAF_THRESHOLD` [ default = 0.01 ] : minor allele frequency threshold when filtering BED files during LD block removal.
* `OUTDIR` [ default = "results/" ] : Output directory for final results. 

## Run
Once your configuration file has been set up, you can run the pipeline at the command-line using nextflow.

This pipeline can be run with different profiles depending on what platform the pipeline is being run on. Current supported platforms/profiles are: eddie, ultra2, local. 

```
nextflow run main.nf -c your_config.config -profile your_profile
```

## Results
Results will be deposited in `results/` by default, unless a different OUTDIR is specified in your configuration.
