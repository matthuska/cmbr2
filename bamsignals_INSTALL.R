# run it using "R(the R version you are using) --slave < bamsignals_INSTALL.R"
# Make sure you are in the same directory where bamsignals_INSTALL.R
# and the directory bamsignals/ are.
# After that the package is installed, to use it just do
# from within R "library(bamsignals)".
system("R CMD build bamsignals")
install.packages("bamsignals_1.1.tar.gz")
system("rm bamsignals_1.1.tar.gz")
