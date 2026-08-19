.libPaths(c("F:/R/win-library/4.6", .libPaths()))
library(urca)
writeLines(deparse(ca.po), "F:/Cpp_Hub/tests/fixtures/timeseries/ca_po_source.txt")
cat("done\n")
