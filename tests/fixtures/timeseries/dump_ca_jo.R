.libPaths(c("F:/R/win-library/4.6", .libPaths()))
library(urca)
writeLines(deparse(ca.jo), "F:/Cpp_Hub/tests/fixtures/timeseries/ca_jo_source.txt")
cat("done\n")
