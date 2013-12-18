# run it using "R(the R version you are using) --vanilla < installBamsignals.R"
# Make sure you are in the same directory where installBamsignals.R
# and the directory bamsignals/ are.
# After that the package is installed, to use it just do
# from within R "library(bamsignals)".
system("R CMD build bamsignals")
install.packages("bamsignals_1.0.tar.gz")
system("rm bamsignals_1.0.tar.gz")
