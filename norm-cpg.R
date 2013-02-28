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

# From Matthias. Calculates some stats about CpG frequencies given a GR object
# specifying the regions and a BSgenome object specifying the genome
cpg_GR <- function(genome, regions) {
  stopifnot(class(regions) == "GRanges")
  chrom = as.character(GenomicRanges::seqnames(regions))
  seq = getSeq(genome, chrom, start(regions), end(regions))

  seq = DNAStringSet(seq)

  # for each seq compute dinucleotide frequencies
  cpg = dinucleotideFrequency(seq)
  cpg = cpg[,"CG"] / apply(cpg, 1, sum)

  # compute the expected cpg content as freq of c*g from a small
  # chromosome (20)
  freq = alphabetFrequency(genome[[20]])

  expected = freq["C"] / sum(freq) * freq["G"] / sum(freq)
  ratio = cpg / expected

  cpg.status = c("lcpg", "hcpg")[as.numeric(ratio > 0.5) + 1]

  meta = elementMetadata(regions)
  if ("ID" %in% colnames(meta)) {
    id = meta[,"ID"]
  } else {
    id = 1:length(regions)
  }

  cpg.tab = data.frame(id=id, cpg, expected=as.numeric(expected), ratio, cpg.status, seq=as.character(seq), stringsAsFactors=F)

  return(cpg.tab)
}
