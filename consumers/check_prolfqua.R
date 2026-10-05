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
