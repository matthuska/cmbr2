# From Jonathan Goeke
# Calculates normalized CpG frequencies for a set of sequences
norm_cpg <- function(seqs) {
    require(BSgenome)
    if (class(seqs) == "character") {
        seqs <- DNAStringSet(seqs)
    }

    nf <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
    dnf <- dinucleotideFrequency(seqs, as.prob = TRUE)
    ncpg <- dnf[,"CG"] / pmax(0.00000000001, (((nf[,"G"] + nf[,"C"]) / 2)^2))
    return(ncpg)
}
