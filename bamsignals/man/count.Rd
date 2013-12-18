\name{count}
\alias{count}
\title{Count reads from a bam file.}
\usage{
  count(gr, bampath, mapqual = 0, shift = 0, ss = F)
}
\arguments{
  \item{gr}{GenomicRanges object used to specify the
  regions}

  \item{bampath}{path to the bam file storing the read. The
  file must be indexed.}

  \item{mapqual}{discard reads with mapping quality
  strictly lower than this parameter. The value 0 ensures
  that no read will be discarded, the value 254 that only
  reads with the highest possible mapping quality will be
  considered.}

  \item{shift}{shift the read position by a user-defined
  number of basepairs. This can be handy in the analysis of
  chip-seq data.}

  \item{ss}{produce a strand-specific count or ignore the
  strand of the read. Strand-specific counts will be
  returned in a 2*length(gr) matrix.}
}
\value{
  a vector or a matrix with the counts
}
\description{
  Count reads in the bins specified by a GenomicRanges
  object. A read position is always specified by its 5'
  end, so a read mapping to the reference strand is
  positioned at its leftmost coordinate, a read mapping to
  the alternative strand is positioned at its rightmost
  coordinate. To change that use the \code{shift}
  parameter.
}

