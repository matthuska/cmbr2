#' Get conservation scores for a given set of genomic ranges
#'
#' The bigWig file of conservation scores can be generated from data
#' downloaded from UCSC. It also requires the command line tools
#' "fetchChromSizes" and "wigToBigWig" from UCSC. These are the
#' commands I used:
#'
#' $ ./fetchChromSizes hg19 > chrom.sizes
#' $ zcat chr{1..22}.phastCons46way.wigFix.gz \
#'     chr{X,Y}.phastCons46way.wigFix.gz | \
#'     ./wigToBigWig -clip stdin chrom.sizes all.bw
#'
#' Download the wigFix phastCons files for hg19 created using 46
#' vertebrate species here:
#'
#' http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate/
#'
#' And download the UCSC command line tools for Linux here:
#'
#' http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
#'
#' FIXME: some regions have no conservation score at all (not even a 0
#' is returned). If the whole range has no score then it is set to
#' NULL, but if a subset of the seqeunce is not available it is
#' sometimes just left out completely. This will mess up aggregate
#' calculations and cause a few other headaches.
#'
#' @param gr GRanges object with the ranges you want conservation data
#' @param cons_file a single bigWig file with conservation data for
#' all chromosomes you are interested in
#' @return a list with two elements:
#' cons = a list of the same length as the GRanges
#' object passed in, where each element is a list of conservation
#' scores associated with the matching region in the supplied GRanges
#' object
#'
#' incomplete = a binary vector that is TRUE if part of a range had no
#' information in the bigWig file. In my eperience this tends to
#' happen at the ends of chromosomes, and in various other places.
#' These missing values can cause errors in aggregate calculations for
#' example, so this vector can be used to subset the conservation list
#' as well as the original GRanges object before further analysis.
#'
getConservation <- function(gr, cons_file) {
  require(rtracklayer)
  cons <- import.bw(cons_file, format = "bw", which = gr, asRangedData = TRUE)
  over <- findOverlaps(gr, cons)
  cons_list <- split(score(cons), queryHits(over))

  # There are some ranges with no conservation data (import.bw() does
  # not return 0 for these regions and there's no obvious way to force
  # that behaviour). We need to add these elements back manually, so
  # that the returned list matches up with the original GRanges
  # object.
  cons_all <- vector("list", length(gr))
  names(cons_all) <- as.character(1:length(gr))
  cons_all[names(cons_list)] <- cons_list

  incomplete <- width(gr) != lapply(cons_all, length)
  stopifnot(length(incomplete) == length(gr))
  return(list(cons=cons_all, notfull=incomplete))
}
