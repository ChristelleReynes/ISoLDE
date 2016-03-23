library(ISoLDE)

targetfile <- system.file("extdata", "target_file.txt", package = "ISoLDE")
rawfile <- system.file("extdata", "rawASRcounts_file.txt", package = "ISoLDE")
normfile <- system.file("extdata", "normASRcounts_file.txt", package = "ISoLDE")

rawASRcounts <- readRawInput(raw_file = rawfile, del = "tab", rownames = TRUE, colnames = TRUE)
normASRcounts <- readNormInput(norm_file = normfile, del = "tab", rownames = TRUE, colnames = TRUE, dec = ".")
target <- readTarget(target_file = targetfile, asr_counts = rawASRcounts, del = "tab")

res_filterT <- filterT(rawASRcounts = rawASRcounts, normASRcounts = normASRcounts, target = target, bias="parental")
filteredASRcounts <- res_filterT$filteredASRcounts
resiso<-isolde_test(bias = "parental", asr_counts = filteredASRcounts, target = target, graph = TRUE, text = FALSE, ext = "png", nboot = 3000)
res_filterT <- filterT(rawASRcounts = rawASRcounts, normASRcounts = normASRcounts, target = target, bias="strain")
filteredASRcounts <- res_filterT$filteredASRcounts
resiso<-isolde_test(bias = "strain", asr_counts = filteredASRcounts, target = target, graph = TRUE, text = FALSE, ext = "png", nboot = 3000, method = "threshold")
