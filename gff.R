
extract <- function (pattern, string, perl = TRUE) {
  r <- paste(".*", pattern, ".*", sep = "")
  matched <- grep(r, string, perl = perl)
  result <- rep(NA, length(string))
  result[matched] <- sub(r, "\\1", string[matched], perl = perl)
  return(result)
}

gff2GR <- function(filename, gffAttrNames=NULL) {
  # read gff into genomic ranges
  require(GenomicRanges)
  regions = scan(filename, what=list(character(), character(), character(), numeric(), numeric(), character(), character(), character(), character()), comment.char="#", sep="\t")

  strand = regions[[7]]
  strand[strand == "."] = "*"
  strand[strand == "1"] = "+"
  strand[strand == "-1"] = "-"

  gr = GRanges(seqnames=regions[[1]],
    ranges=IRanges(start=regions[[4]], end=regions[[5]]),
    strand=strand)

  src = regions[[2]]
  type = regions[[3]]
  score = regions[[6]]
  df = DataFrame(src, type, score)

  if (!is.null(gffAttrNames)) {
    df = cbind(df, DataFrame(sapply(gffAttrNames, function(n)
      extract(paste(n , "=(.+?)(;|$)", sep=""), regions[[9]]))))
  }
  elementMetadata(gr) = df

  return(gr)
}



GR2gff <- function(regions, filename, feature.type="experimental_feature", src="GenomicRanges", score=".", phase=".") {
  require(GenomicRanges)

  strnd = as.character(strand(regions))
  strnd[strnd == "*"] = "."

  tab = data.frame(as.character(seqnames(regions)), src, feature.type, as.numeric(start(regions)), as.numeric(end(regions)), score, strnd, phase, makeGffAttributes(as.data.frame(elementMetadata(regions))), stringsAsFactors=F)

  write.table(tab, file=filename, sep="\t", quote=F, row.names=F, col.names=F)
}

makeGtfAttributes <- function(df, cols=NULL) {
  if (is.null(cols))
    cols = colnames(df)
  # make sure that gene_id and transcript_id are the first two columns
  mandatory = c("gene_id", "transcript_id")
  o = match(c(mandatory, setdiff(cols, mandatory)), cols)
  if (any(is.na(o[1:length(mandatory)]))) {
    warning("mandatory gtf attributes gene_id or transcript_id missing")
    o = o[!is.na(o)]
  }
  cols = cols[o]
  return(paste(apply(sapply(cols, function(s) {
    content = df[,s]
    if (is.character(content) | is.factor(content)) {
      content = paste('"', content, '"', sep="")
    }
    paste(gsub(".", "_", s, fixed=T), content, sep=" ")
  }), 1, paste, collapse="; "), ";", sep=""))
}


gtf2GR <- function(filename, gtfAttrNames=NULL) {
  # read gff into genomic ranges
  require(GenomicRanges)
  regions = scan(filename, what=list(character(), character(), character(), numeric(), numeric(), character(), character(), character(), character()), comment.char="#", sep="\t")

  strand = regions[[7]]
  strand[strand == "."] = "*"

  gr = GRanges(seqnames=regions[[1]],
    ranges=IRanges(start=regions[[4]], end=regions[[5]]),
    strand=strand)

  src = regions[[2]]
  type = regions[[3]]

  df = data.frame(src, type, stringsAsFactors=F)
  if (!is.null(gtfAttrNames)) {
    df = data.frame(df, sapply(gtfAttrNames, function(n)
      extract(paste(n , " (.+?)(;|$)", sep=""), regions[[9]])),
      stringsAsFactors=F)
  }
  elementMetadata(gr) = df
  return(gr)
}

GR2gtf <- function(regions, filename, feature.type="experimental_feature", src="GenomicRanges", score=".", phase=".", attributes=NULL, ...) {
  require(GenomicRanges)

  strnd = as.character(strand(regions))
  strnd[strnd == "*"] = "."

  tab = data.frame(as.character(seqnames(regions)), src, feature.type, as.numeric(start(regions)), as.numeric(end(regions)), score, strnd, phase, makeGtfAttributes(as.data.frame(elementMetadata(regions)), cols=attributes), stringsAsFactors=F)

  write.table(tab, file=filename, sep="\t", quote=F, row.names=F, col.names=F, ...)
}

parseProperBEDLine <- function(bedline) {
  lineArgs <- strsplit(bedline, "\t")[[1]]
  argTypes <- lapply(lineArgs, function(arg) type.convert(arg, as.is=TRUE))
  if (length(lineArgs)>=3){
    if (is.integer(argTypes[[2]])  &&  is.integer(argTypes[[3]])){
      if (length(lineArgs)>=5){
        if (is.numeric(argTypes[[5]])){
          if (all(length(lineArgs)>=6, !(lineArgs[6] %in% c("+", "-","*",".")))){
            return (NA)
          }
          else{
            return(argTypes)
          }
        }
        else{
          return (NA)
        }
      }
      else {
        return(argTypes)
      }
    }
  }
  return (NA)
}

bed2GR2 <- function(filename, parseMetadata=TRUE) {
  #tries to guess the right number of columns and header lines to be skipped
  #should adapt the function bed2GR so that it can parse the metaData
  nLines <- 10
  firstLines <- scan(filename, what=character(), nlines = nLines, sep="\n")
  what <- NA
  lIndex <- 0
  while (all(is.na(what),lIndex < nLines) ){
    lIndex <- lIndex + 1
    what <- parseProperBEDLine(firstLines[lIndex])
  }

  if (lIndex >= nLines){
    stop("unable to find a proper bed line in the file")
  }
  
  if (!parseMetadata) what <- what[1:min(length(what), 6)]
  #todo: try to parse possible headers to get the column names right
  if (length(what)>6) names(what)[7:length(what)] <- paste("metadata",1:(length(what)-6), sep="_")
  bed2GR(filename, what=what, skip=lIndex-1)

}

# write a GR object into a bed file and all meta data as additional columns
GR2bed <- function(regions, filename, header=FALSE, writeMetadata=TRUE) {
  require(GenomicRanges)
  tab = data.frame(chrom=as.character(seqnames(regions)), start=start(regions)-1, end=end(regions))
  
  fieldNum = 3
  extraColNames <- names(elementMetadata(regions))[!names(elementMetadata(regions))%in%c("name","score")]
  if (length(extraColNames)>0 && writeMetadata) fieldNum = 7 #in this case the number of fields can be also above 7, not necessarily 7
  else if ((!is.null(strand(regions))) && any(strand(regions)!="*")) fieldNum = 6
  else if (!is.null(score(regions))) fieldNum = 5
  else if (!is.null(names(regions))) fieldNum = 4
  
  
  
  
  if (fieldNum > 3){
    if (!is.null(names(regions))) tab$name=names(regions)
    else tab$name=rep("*", length(regions))
  
    if (fieldNum > 4){
      if (!is.null(score(regions))) tab$score=score(regions)
      else tab$score=rep("*", length(regions))
    
      if (fieldNum > 5){
        strnd <- rep("*", length(strand(regions)))
        strnd[as.logical(strand(regions)=="+")] <- "+"
        strnd[as.logical(strand(regions)=="-")] <- "-"
        tab$strand=strnd
      }
      
      if (fieldNum > 6){
        tab <- data.frame(tab, as(elementMetadata(regions)[extraColNames], "data.frame"))
      }
    }
  }
  
  cnames <- F
  if (header) {cnames <- names(tab); cnames[1] <- paste("#", cnames[1], sep="")}
  write.table(tab, file=filename, sep="\t", quote=F, row.names=F, col.names=cnames)
}


bed2GR <- function(filename, nfields=6, skip=0, what=NA) {
  stopifnot(nfields >= 3)
  # read BED into genomic ranges
  require(GenomicRanges)
  if (!is.list(what) && is.na(what)){
    what = list(character(), numeric(), numeric(), character(), numeric(), character())[1:nfields]
  }
  else {
    nfields = length(what)
  }
  
  regions = scan(filename, what=what, sep="\t", skip=skip)

  if (nfields >= 6) {
    strand = regions[[6]]
    strand[strand == "."] = "*"
  } else {
    strand = "*"
  }

  
  # GRanges are 1-indexed and closed, while BED intervals are 0-indexed and half-open
  gr = GRanges(seqnames=regions[[1]],
  ranges=IRanges(start=regions[[2]]+1, end=regions[[3]]), strand=strand)
  

  if (nfields >= 4 && any(regions[[4]]!="*")) {
    names(gr) = regions[[4]]
  }
  
  if (nfields > 6) {
    elementMetadata(gr) = DataFrame(score=regions[[5]], regions[7:length(regions)])
  }
  else if (nfields >= 5) {
    elementMetadata(gr) = DataFrame(score=regions[[5]])
  }
  
  
  
  return(gr)
}

makeGffAttributes <- function(df, cols=NULL) {
  if (ncol(df) == 0)
    return(rep("", nrow(df)))
  if (is.null(cols))
    cols = colnames(df)
  return(apply(sapply(cols, function(s) paste(gsub(".", "_", s, fixed=T), df[,s], sep="=")), 1, paste, collapse=";"))
}

countBamInGRanges <- function(bam.file, granges, min.mapq=NULL, read.width=1) {
  require(GenomicRanges)
  require(Rsamtools)

  rds.counts <- numeric(length(granges))
  seq.names <- as.character(unique(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      print( paste("[", Sys.time(),"] Started processing count on chromosome", seq.name, "of file", bam.file) )
      granges.subset <- granges[seqnames(granges)==seq.name]
      strand(granges.subset) <- "*"
      rds <- scanBam(bam.file,param=ScanBamParam(what=c("pos","mapq"),which=range(granges.subset)))
      if (!is.null(min.mapq)) {
        mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      } else {
        mapq.test = rep(T, length(rds[[1]]$mapq))
      }
      if (sum(mapq.test) > 0) {
        rds.ranges <- GRanges(seq.name,IRanges(start=rds[[1]]$pos[mapq.test],width=read.width))
        rds.counts.seq.name <- countOverlaps(granges.subset,rds.ranges)
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- rds.counts.seq.name
      } else {
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
      }
    } else {
      rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
    }
    print( paste("[", Sys.time(),"] Finished processing count on chromosome", seq.name, "of file", bam.file) )
  }
  rds.counts
}

# Using countBam() is fast but a little messy: if you have overlapping ranges
# then you'll get duplicate counts (apparently). However if your ranges are
# reasonably distant from each other and you're okay with counting all reads
# that partially overlap the range, then this method is very fast. (~ 30x faster
# than countBamInGranges())
#
# Also, we do not filter by mapping quality (min.mapq) or allow the modification
# of the read width.
countBamInGRangesFast <- function(bam.file, granges) {
  require(GenomicRanges)
  require(Rsamtools)

  rds.counts <- numeric(length(granges))
  seq.names <- as.character(unique(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
  fields <- c("pos")
  strand(granges) <- "*"
  cnts <- countBam(bam.file,param=ScanBamParam(what=fields,which=granges))
  cnts$records
}

getBins <- function(chr=NULL, n=NULL, bin.size=NULL, genome=Rnorvegicus, offset=0) {
  stopifnot(!all(c(is.null(n), is.null(bin.size)), "specify either bin size or number of bins"))
  if (is.null(chr)) {
    chr = seqnames(genome)
  }
  if (!is.null(n)) {
    bin.size = floor((seqlengths(genome)[chr] - offset) / n)
    names(bin.size) = chr
    n = rep(n, length(chr))
    names(n) = chr
  } else {
    n = floor((seqlengths(genome)[chr] - offset) / bin.size)
    names(n) = chr
    bin.size = rep(bin.size, length(chr))
    names(bin.size) = chr
  }

  g = GRanges()
  for (ch in chr) {
    g = c(g, GRanges(seqnames=ch, IRanges(start=0:(n[ch] - 1) * bin.size[ch] + 1 + offset, width=bin.size)))
  }
  return(g)
}

# From Matthias 2013-01-15
# Requires all GRanges to have the same width
# Suggestions to make faster: use reduce on the granges.subset and ranges
# instead of range
coverageBamInGRanges <- function(bam.file, granges, min.mapq, reads.collapsed=FALSE, width=NULL) {
  require(GenomicRanges)
  require(Rsamtools)

  # first check that all granges have the same width
  w = width(granges[1])
  stopifnot(all(width(granges) == w))

  seq.names <- as.character(unique(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)

  grange.coverage = matrix(0, nrow=length(granges), ncol=w)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      print( paste("[", Sys.time(),"] Started processing coverage on chromosome", seq.name, "of file", bam.file) )
      granges.subset <- granges[seqnames(granges)==seq.name]
      strand(granges.subset) <- "*"
      what = c("pos", "mapq", "qwidth")
      if (reads.collapsed) {
        what = c(what, "qname")
      }
      rds <- scanBam(bam.file,param=ScanBamParam(what=what, which=range(granges.subset)))
      if (missing(min.mapq)) {
        mapq.test = rep(T, length(rds[[1]]$mapq))
      } else {
        mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      }
      if (sum(mapq.test) > 0) {
        if (is.null(width)) {
          width = rds[[1]]$qwidth[mapq.test]
        }
        rds.ranges <- GRanges(seq.name, IRanges(start=rds[[1]]$pos[mapq.test], width=width))
        if (reads.collapsed) {
          multiply = as.numeric(sapply(strsplit(rds[[1]]$qname[mapq.test], "_x"), "[", 2))
          select = unlist(lapply(1:length(multiply), function(i) rep(i, multiply[i])))
          rds.ranges = rds.ranges[select]
        }
        # set the seqlength, so the coverage Rle gets the right length
        len = seqlengths(granges)[seq.name]
        if (is.na(len)) {
          len = max(c(end(granges), end(rds.ranges)))
        }
        seqlengths(rds.ranges)[seq.name] = len

        coverage.seq.name <- coverage(rds.ranges)[[1]]
        v = Views(coverage.seq.name, start=start(granges.subset), end=end(granges.subset))
        cvg = t(sapply(v, as.numeric))
        grange.coverage[as.logical(seqnames(granges)==seq.name),] <- cvg
      }
      print( paste("[", Sys.time(),"] Finished processing coverage on chromosome", seq.name, "of file", bam.file) )
    }
  }
  # reverse the ones on the minus strand
  minus = as.logical(strand(granges) == "-")
  if (any(minus)) {
    grange.coverage[minus,] = t(apply(grange.coverage[minus,], 1, rev))
  }
  return(grange.coverage)
}

##' A fast(er) function to calculate coverage across a set of GRanges.
##'
##' An alternative version of the coverageBamInGRanges() function. Unlike that
##' function you can not filter the reads by mapping quality or collapse reads
##' (what does that even do?). Also unlike that function this one is pretty
##' fast. Calculating coverage across 6 bam files for 33,000 ranges, each of
##' which was 3kb in length took about 4.5 minutes with this function and
##' approximately an hour with the old function.
##'
##' TODO: add an option to return positive and negative strand coverage
##' separately.
##'
##' TODO: allow filtering by min.mapq
##'
##' TODO: consider the strand of the reads. Right now the start of the read is
##' always the "left-most" or "lowest" position (as returned by scanBam)
##'
##' @param bam.file the full path to a bam file. It should have an associated
##' index with the same name and .bai at the end.
##' @param granges a GRanges object where all ranges have the same width
##' @param frag.width an optional parameter to force the fragment width to
##' something other than the read width contained in the bam file. By default
##' the qwidth contained in the bam file is used.
##' @return a length(granges) x width(granges) dimension matrix where each row
##' is a grange and each column is a base pair relative to the start of the
##' GRange.
##' @author Matthew Huska
coverageBamInGRangesFast <- function(bam.file, granges, frag.width=NULL) {
  require(GenomicRanges, quietly=TRUE)
  require(Rsamtools, quietly=TRUE)
  # first check that all granges have the same width
  w <- width(granges[1])
  stopifnot(all(width(granges) == w))
  what <- c("pos", "mapq", "qwidth")
  if (is.null(frag.width)) {
    rds <- scanBam(bam.file, param=ScanBamParam(what=what, which=granges))
  } else {
    # Extend the width of our granges by "frag.width" in order to catch overlaps
    # that are further away from the ends of our ranges than the actual read
    # width (qwidth). This only matters when frag.width > qwidth.
    rds <- scanBam(bam.file, param=ScanBamParam(what=what, which=resize(granges, width(granges) + 2*frag.width, fix="center")))
  }
  read_pos <- lapply(rds, "[[", "pos")
  # Position of the read relative to the start of the grange
  relative_pos <- mapply("-", read_pos, start(granges) - 1)
  starts <- lapply(relative_pos, function(s) { s[s < 1] <- 1; s[s > w] <- w; s })
  if (is.null(frag.width)) {
    widths <- lapply(rds, "[[", "qwidth")
    ends <- mapply(function(x,y,w) {e <- x + y; e[e > w] <- w; e[e < 1] <- 1; e}, relative_pos, widths, w)
  } else {
    ends <- mapply(function(x,y,w) {e <- x + y; e[e > w] <- w; e[e < 1] <- 1; e}, relative_pos, frag.width, w)
  }
  start_counts <- lapply(starts, tabulate, nbins=w)
  end_counts <- lapply(ends, tabulate, nbins=w)
  start_sums <- lapply(start_counts, cumsum)
  end_sums <- lapply(end_counts, cumsum)
  grange.coverage <- do.call(rbind, start_sums) - do.call(rbind, end_sums)

  # reverse the ones on the minus strand
  minus <- as.logical(strand(granges) == "-")
  if (any(minus)) {
    grange.coverage[minus,] <- t(apply(grange.coverage[minus,], 1, rev))
  }

  invisible(grange.coverage)
}

# Window based coverage counting
#
# your highness 2013-01-25 helmuth@molgen.mpg.de
#
# Equally sized granges necessary
#
# TODO: Maybe it is better to use countBamInGRanges for this purpose to count for bins
#
# Returns a list of data.frames with window counts for each grange
coverageBamInGRangesWindows  <- function( bam.file, granges, window.width=300, sliding.window=F, FUN=coverageBamInGRanges, ...) {
  require( GenomicRanges )

  coverage = FUN( bam.file=bam.file, granges=granges, ... )

  if ( sliding.window ) {
    half.win.size = window.width / 2
    window.count = width( granges )[1] / half.win.size - 1
    window.coverage = do.call( "cbind", lapply(1:window.count, function(i) { rowSums(coverage[ , ((i-1) * half.win.size + 1):((i-1) * half.win.size + window.width)]) }) )
  } else {
    window.count = width( granges )[1] / window.width
    window.coverage = do.call( "cbind", lapply(1:window.count, function(i) { rowSums(coverage[ , ((i-1) * window.width + 1):(i * window.width)]) }) )
  }

  return (window.coverage)
}

# Parallel processing for a list of bams (Thanks to Mike Love)
#
# your highness 2013-01-25 helmuth@molgen.mpg.de
#
# Takes in a list of bamfiles and does counting on multiple processors
# The number of processors in determined by the length of the list of bamfiles or directly specified with mc.cores argument
# The function is provided via FUN parameter. Default is countBamInGRangesFast
# Extra parameters are given to FUN
#
# Returns a list of vectors or data.frames (depending on the plugged in function)
processListOfBamsInGRanges  <- function( bam.files, granges=granges, mc.cores=NA, FUN=countBamInGRangesFast, ... ) {
  require(multicore)

  if ( is.na(mc.cores) ) {
    mc.cores = length( bam.files )
  }
  counts = mclapply( 1:length(bam.files), function( i ) {
    print( paste("[", Sys.time(), "] Processor", i ,": Retrieving tag count for", bam.files[i], ".") )
    FUN( bam.file=bam.files[i], granges=granges, ... )
  }, mc.cores=mc.cores )

  names(counts) = sapply( bam.files, function( file ) { tail( unlist(strsplit( file, "/", fixed=T)), 1 ) } )

  return (counts)
}
