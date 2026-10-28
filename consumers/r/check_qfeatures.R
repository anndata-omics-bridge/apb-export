# QFeatures links apb-export's levels and aggregates precursors to proteins.
# Inside the R consumer image: Rscript /checks/check_qfeatures.R FILE...

suppressPackageStartupMessages(library(QFeatures))

# Each modality links to the next coarser one present, joining its rowData on that one's key.
COARSER <- list(precursors = c("peptides", "proteins"), peptides = "proteins")
KEY <- c(peptides = "sequence", proteins = "proteins")

# readH5MU returns each modality as a SingleCellExperiment whose X is an unnamed assay;
# coercing it with as(x, "SummarizedExperiment") drops the rowData, so build the assay here.
# Categorical var columns arrive as factors; the join compares them as strings.
as_assay <- function(x) {
  features <- rowData(x)
  for (column in names(features)) {
    if (is.factor(features[[column]])) features[[column]] <- as.character(features[[column]])
  }
  SummarizedExperiment(assays = list(intensity = assay(x)), rowData = features)
}

check <- function(path) {
  mae <- MuData::readH5MU(path)
  assays <- lapply(MultiAssayExperiment::experiments(mae), as_assay)
  # The annotations sit in each modality's obs and the MuData obs is empty; one modality's
  # colData serves all, so aggregateFeatures meets no conflicting colData.
  qf <- QFeatures(assays, colData = colData(mae[[1]]))
  links <- character()
  for (from in intersect(names(COARSER), names(qf))) {
    to <- intersect(COARSER[[from]], names(qf))[1]
    if (is.na(to)) next
    qf <- addAssayLink(qf, from = from, to = to, varFrom = KEY[[to]], varTo = KEY[[to]])
    # QFeatures keeps the link on the coarser assay; its hits run from the finer one.
    linked <- length(unique(S4Vectors::queryHits(assayLink(qf, to)@hits)))
    if (linked != nrow(qf[[from]])) {
      stop(path, ": ", nrow(qf[[from]]) - linked, " ", from, " link to no ", to)
    }
    links <- c(links, sprintf("%s>%s %d", from, to, linked))
  }
  qf <- aggregateFeatures(
    qf, "precursors", fcol = "proteins", name = "aggregated", fun = colSums, na.rm = TRUE
  )
  proteins <- nrow(qf[["aggregated"]])
  if (proteins == 0) stop(path, ": aggregateFeatures returned no proteins")
  cat(basename(path), ": ", paste(links, collapse = ", "), "; ", proteins,
      " proteins from precursors\n", sep = "")
}

for (path in commandArgs(trailingOnly = TRUE)) check(path)
