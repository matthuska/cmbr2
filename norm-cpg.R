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

#' Computes GC content and CpGodd for a given GRanges object
#'
#' As a threshold between HCP and LCP a CpG odd ratio of .4 is a good approximate.
#'
#' You can construct a 3kb-flanking-TSS representation for a given set of granges with:
#'
#'    require(BSgenome)
#'    require(BSgenome.Hsapiens.UCSC.hg19)
#'    require(GenomicFeatures)
#'    granges.tss = flank(granges, 1500, start=T)
#'    
#'    computeCpGOddInGRanges( Hsapiens, granges, tss)
#'
computeCpGOddInGRanges <- function(Hsapiens=Hsapiens, granges=granges) {

}

