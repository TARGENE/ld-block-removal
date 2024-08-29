process FlashPCA {
    label "multithreaded"
    label 'flashpca_image'

    input:
        path bedfiles
    
    output:
        path "pcs.txt", emit: pcs
        path "pve.txt", emit: pve
    
    script:
        prefix = bedfiles[0].toString().minus('.bed')
        "/home/flashpca-user/flashpca/flashpca --bfile $prefix --ndim $params.NB_PCS --numthreads $task.cpus"
}

process AdaptFlashPCA {
    label 'tl_core_image'
    publishDir "$params.OUTDIR/covariates/", mode: 'symlink'
    label 'bigmem'
    
    input:
        path flashpca_out
    
    output:
        path "pcs.csv"
    
    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input $flashpca_out --output pcs.csv adapt
        """
}

process ScreePlot {
    label 'bgen_python_image'
    publishDir "$params.OUTDIR/pca/", mode: 'symlink'

    input:
        path pve

    output:
        path "ScreePlot.pdf"

    script:
        """
        #!/usr/bin/env Rscript
        library(ggplot2)

        pve <- read.table("$pve")
        pve[['PC']] <- 1:nrow(pve)
        colnames(pve)[1] <- "pve"

        pdf("ScreePlot.pdf")
        ggplot(data = pve, mapping = aes(x = PC, y = pve)) +
          geom_point() +
          geom_line() +
          ylab("Proportion of Variance Explained") +
          theme_bw()
        dev.off()
        """

}

