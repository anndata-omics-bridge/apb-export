# The R tools that read apb-export's files: QFeatures and MSstats, on the Bioconductor release
# their export rules' target_version_pattern names. MuData comes from GitHub at the merge of
# MuData#11, which reads the nullable strings and None values anndata writes; no Bioconductor
# release carries it yet.
FROM bioconductor/bioconductor_docker:RELEASE_3_23

RUN Rscript -e 'BiocManager::install(c("remotes", "QFeatures", "MSstats"), ask = FALSE, update = FALSE)' \
 && Rscript -e 'BiocManager::install(c("ilia-kats/MuData@9314788d61b4c04a38a4239a5201f60aab324ed4"), ask = FALSE, update = FALSE)' \
 && Rscript -e 'stopifnot(packageVersion("QFeatures") >= "1.22", packageVersion("MSstats") >= "4.20")'

COPY check_qfeatures.R check_msstats.R /checks/
