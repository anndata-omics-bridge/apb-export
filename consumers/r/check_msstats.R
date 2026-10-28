# MSstats summarizes apb-export's long table to proteins and compares the two conditions.
# Inside the R consumer image: Rscript /checks/check_msstats.R FILE

suppressPackageStartupMessages(library(MSstats))

path <- commandArgs(trailingOnly = TRUE)[[1]]
input <- utils::read.csv(path)
proteins <- length(unique(input$ProteinName))

# dataProcess censors the lowest intensities by default, so a protein measured only there,
# such as the synthetic P1, may get no summary; the check needs some, not all.
processed <- dataProcess(input, use_log_file = FALSE, verbose = FALSE)
summarized <- length(unique(processed$ProteinLevelData$Protein))
if (summarized == 0) stop(path, ": dataProcess summarized none of ", proteins, " proteins")

conditions <- levels(processed$ProteinLevelData$GROUP)
contrast <- matrix(c(-1, 1), nrow = 1, dimnames = list(paste(rev(conditions), collapse = "-"), conditions))
compared <- groupComparison(contrast, processed, use_log_file = FALSE, verbose = FALSE)
rows <- nrow(compared$ComparisonResult)
if (rows != summarized) stop(path, ": groupComparison returned ", rows, " of ", summarized, " proteins")

cat(basename(path), ": ", nrow(input), " rows; ", summarized, " of ", proteins,
    " proteins summarized and compared (", rownames(contrast), ")\n", sep = "")
