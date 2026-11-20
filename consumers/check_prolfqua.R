# prolfquapp rebuilds LFQData from apb-export's prolfqua file, hierarchy and factors intact.
# Inside the consumer image: Rscript check_prolfqua.R EXAMPLES

args <- commandArgs(trailingOnly = TRUE)
adata <- anndataR::read_h5ad(file.path(args[[1]], "prolfqua.h5ad"))
restored <- prolfquapp::LFQData_from_anndata(adata)
config <- restored$lfqdata$get_config()

# HDF5 may return uns mappings in name order; the hierarchy must keep its own order.
hierarchy <- config$hierarchy_keys()
if (!identical(hierarchy, c("protein_Id", "peptide_Id", "precursor_Id"))) {
  stop("hierarchy came back as ", paste(hierarchy, collapse = ", "))
}
if (!identical(config$factor_keys(), "condition")) {
  stop("factors came back as ", paste(config$factor_keys(), collapse = ", "))
}
cat(
  "prolfquapp rebuilds LFQData:", nrow(restored$lfqdata$data), "rows,",
  paste(hierarchy, collapse = " > "), "\n"
)

# prolfquapp's APB reader filters by the q-value layers and sums precursors to peptides.
samples <- utils::read.delim(file.path(args[[1]], "samples.tsv"))
annotation <- prolfquapp::read_annotation(
  data.frame(file = samples$raw_file, group = samples$condition),
  QC = TRUE
)
xd <- prolfquapp::preprocess_APB(file.path(args[[1]], "prolfqua.h5ad"), character(), annotation)
counts <- xd$lfqdata$hierarchy_counts()
if (!identical(xd$lfqdata$hierarchy_keys(), c("protein_Id", "peptide_Id")) || counts$protein_Id == 0) {
  stop("APB reader returned ", paste(xd$lfqdata$hierarchy_keys(), collapse = ", "), " and ", counts$protein_Id, " proteins")
}
cat("prolfquapp APB reader:", counts$protein_Id, "proteins,", counts$peptide_Id, "peptides\n")

# The APB reader gives the same peptide values, and medpolish the same protein values, as
# prolfquapp's own DIA-NN reader on report.tsv.
reference <- prolfquapp::preprocess_DIANN(file.path(args[[1]], "report.tsv"), character(), annotation)
values <- function(lfq) {
  keys <- c(lfq$sample_name(), lfq$hierarchy_keys())
  data <- lfq$data_long()
  data <- data[do.call(order, data[keys]), c(keys, lfq$response())]
  names(data)[ncol(data)] <- "value"
  as.data.frame(data)
}
protein <- function(lfq) {
  lfq$get_Transformer()$log2()$lfq$get_Aggregator("medpolish")$aggregate()
}
for (level in c("peptide", "protein")) {
  ours <- xd$lfqdata
  theirs <- reference$lfqdata
  if (level == "protein") {
    ours <- protein(ours)
    theirs <- protein(theirs)
  }
  same <- all.equal(values(ours), values(theirs), check.attributes = FALSE, tolerance = 1e-12)
  if (!isTRUE(same)) {
    stop(level, " values differ from preprocess_DIANN: ", paste(same, collapse = "; "))
  }
  cat("prolfquapp", level, "values:", nrow(values(theirs)), "rows, equal to preprocess_DIANN\n")
}
