#include <Rcpp.h>
#include <stdio.h>  
#include <algorithm>
#include "sam.h"  

//in C: true==1(or something different than 0), false==0
//coding conventions:
//iterator for type t and variable v: It, variable: i_v,
//end of the iterator, e_v.
//use classes, not structs. Name of a class starts with upper-case letter
using namespace Rcpp;

inline bool isNegStrand(const bam1_t *b){
	return ((b->core).flag & BAM_FREVERSE) != 0;
}

inline int fiveprimepos(const bam1_t *b, int strand){
	if (strand){//read on the negative strand
		return bam_calend(&(b->core), bam1_cigar(b))  - 1;
	} else{//read on the reference strand
		return (b->core).pos;
	}
}



typedef NumericVector::iterator Idouble;
typedef IntegerVector::iterator Iint;


inline int getRefId(samfile_t* in, const std::string& refname){
	return bam_get_tid(in->header, refname.c_str());
}

//genomic array: represents a correspondence between a genomic range and
//a memory area.
class GArray {
	//initialized when parsing the GRanges
	public:
	
	int rid;//reference id as in the bamfile
	int loc;//0-indexed
	int len;
	int strand;//-1,0,1 represent alternative, unspecified, reference strand respectively
	inline int end(){ return loc + len; }
	//initialized when allocating memory
	int* array;
	int alen;
	GArray(int _rid, int _loc, int _len, int _strand){
		rid = _rid;
		loc = _loc;
		len = _len;
		strand = _strand;
	}
};

//this object won't be valid anymore when the RObject used for construction
//will be destroyed
class RleIter {
	IntegerVector rlens;
	IntegerVector values;
	CharacterVector names;
	
	public:
	
		int run;
		int rlen;
		int rpos;
		bool valid;

		RleIter(RObject& rle):
			rlens(as<IntegerVector>(rle.slot("lengths"))),
			values(as<IntegerVector>(rle.slot("values"))),
			names(as<CharacterVector>(values.attr("levels"))),
			run(0), rpos(-1)
		{
			next();
		}
		
		bool next(){
			++rpos;
			if (rpos == rlens[run]){ //end of the run, go to the next
				++run; rpos = 0;
				if (run == rlens.length())
					valid = false;
					return valid;
			}
			valid = true;
			return valid;
		}
		
		String getValue(){
			return names[values[run]-1];
		}
};

//parses the GR object.
void parseRegions(std::vector<GArray>& container, RObject& gr, samfile_t* in){
	if (not gr.inherits("GRanges"))
		stop("must provide a GRanges object");
	
	IntegerVector starts = as<IntegerVector>(as<RObject>(gr.slot("ranges")).slot("start"));
	IntegerVector lens =   as<IntegerVector>(as<RObject>(gr.slot("ranges")).slot("width"));
	
	RObject chrsRle = as<RObject>(gr.slot("seqnames"));
	RObject strandsRle = as<RObject>(gr.slot("strand"));
	RleIter chrs(chrsRle);
	RleIter strands(strandsRle);
	container.reserve(container.size() + starts.length());
	Iint e_starts = starts.end(); Iint i_starts = starts.begin(); Iint i_lens = lens.begin();
	
	int lastStrandRun = -1;
	int strand = -1;
	
	int lastChrsRun = -1;
	int rid = -1;
	 
	for (; i_starts < e_starts; ++i_starts, ++i_lens, chrs.next(), strands.next()){
		//if new run, update chromosome
		if (lastChrsRun != chrs.run){
			lastChrsRun = chrs.run;
			rid = getRefId(in, chrs.getValue());
			if (rid == -1)
				stop("chromosome " + (std::string)chrs.getValue() + " not present in the bam file");
		}
		
		//if new run, update strand 
		if (lastStrandRun != strands.run){
			lastStrandRun = strands.run;
			const std::string& s = strands.getValue();
			if (s == "-"){ strand = -1; }
			else if (s == "+"){ strand = +1; }
			else { strand = 0; }
		}
		
		container.push_back(GArray(rid, *i_starts - 1, *i_lens, strand));
	}
}

static List allocateAndDistributeMemory(std::vector<GArray>& ranges, int binsize, bool ss){
	int rnum = ranges.size();//number of ranges
	int mult = ss?2:1;
	int pos = 0;
	
	
	IntegerVector starts(rnum);//start of each region in the returned count vector
	IntegerVector ends(rnum);//end of each region in the returned count vector
	
	typedef std::vector<GArray>::iterator IGArray;
	IGArray i_ranges = ranges.begin(); IGArray e_ranges = ranges.end();
	Iint i_starts = starts.begin(); Iint i_ends = ends.begin();
	for (; i_ranges < e_ranges; ++i_ranges, ++i_ends, ++i_starts){
		*i_starts = pos+1;//in R indices are 1-indexed... 
		pos += mult*ceil((i_ranges->len)/((float)binsize));
		*i_ends = pos;//and intervals are right-closed
	}
	
	IntegerVector counts(pos);
	int* countsptr = counts.begin();
	
	i_ranges = ranges.begin(); i_starts = starts.begin(); i_ends = ends.begin();
	for (;i_ranges < e_ranges; ++i_ranges, ++i_ends, ++i_starts){
		i_ranges->array = countsptr + *i_starts - 1; 
		i_ranges->alen = *i_ends - *i_starts + 1;
	}
	
	return List::create(_("counts")=counts, _("starts")=starts, _("ends")=ends);
}

//non-standard c++ semantics. Close should be the destructor and
//you should implement copy constructor and copy assignment operator
//according to the rule of three. This way you have automatic memory management.
class Bamfile {
	public:
		samfile_t* in;
		bam_index_t* idx;
		//allocate memory
		Bamfile(const std::string& bampath){
			const char* cbampath = bampath.c_str();
			in = samopen(cbampath, "rb", 0);  
			if (in == 0) {  
				stop("Fail to open BAM file " + bampath);  
			}
			
			idx = bam_index_load(cbampath); // load BAM index  
			if (idx == 0) {  
				stop("BAM indexing file is not available for file " + bampath);
			}  
			bam_init_header_hash(in->header);
		}
		//deallocate
		void close(){
			bam_index_destroy(idx);  
			samclose(in);
		}
};

template <class TRegion>
static bool sortByStart(const TRegion& a, const TRegion& b){
	int ret = b.rid - a.rid;
	if (ret==0){ return b.loc > a.loc; }
	return ret > 0;
}

//interface for TRegion:
//int rid;
//int loc; start position
//int end(); end position
//interface for TPileup:
//void pileup(TRegion&, const bam1_t*, int, int)
template <class TRegion, class TPileup>
static void overlapAndPileup(Bamfile& bfile, std::vector<TRegion>& ranges, int mapqual, int shift, TPileup& pileupper, int maxgap){
	
	//sorting intervals according to start coordinate (and ref id of course)
	std::sort(ranges.begin(), ranges.end(), sortByStart<TRegion>);
	
	//trade-off between querying a new region and processing unnecessary reads
	
	const int MAX_GAP = maxgap;
	//variables processed, chunk_start, chunk_end, curr_range and range are indices for the vector ranges
	unsigned int processed = 0;
	bam1_t* read; read = bam_init1();
	//process one chunk of nearby ranges at a time
	while (processed < ranges.size()){
		unsigned int chunk_start = processed;
		int rid = ranges[chunk_start].rid;
		int start = ranges[chunk_start].loc - shift;
		//find out how many regions to process together
		unsigned int chunk_end = chunk_start+1;
		for (; chunk_end < ranges.size(); ++chunk_end){
			if (ranges[chunk_end].rid != rid || ranges[chunk_end].loc - ranges[chunk_end-1].end() - 2*shift > MAX_GAP){
				break;
			}
		}
		int end = ranges[chunk_end-1].end() + shift;
		//perform query
		//Rcout << "querying: " << rid << ":" << start << "-" << end << std::endl;
		//Rcout << "for the chunk: " << chunk_start << ":" << chunk_end << std::endl;
		
		bam_iter_t iter = bam_iter_query(bfile.idx, rid, start, end);
		//all ranges behind curr_range should not overlap with the next reads anymore
		unsigned int curr_range = chunk_start;
		//loop through the reads
		while (bam_iter_read((bfile.in)->x.bam, iter, read) >= 0){
			if ((read->core).qual >= mapqual){
				//++deleteme2;
				int r_start = (read->core).pos;
				//skip non-overlapping regions at the beginning
				while (curr_range < chunk_end && r_start >= ranges[curr_range].end() + shift) ++curr_range;
				//should never happen, unless the last reads returned by the iterator no not overlap with the queried interval
				if (curr_range == chunk_end) break; 
				
				int r_end = bam_calend(&(read->core), bam1_cigar(read)) -1;
				
				for (unsigned int range = curr_range; range < chunk_end && ranges[range].loc - shift <= r_end; ++range){
					pileupper.pileup(ranges[range], read, r_start, r_end);
				}
			}
			
		}
		
		bam_iter_destroy(iter);
		processed = chunk_end;
	}
	//Rcout << "found: " << deleteme << " overlaps out of " << deleteme2 << " scanned reads!!!" << std::endl;
	bam_destroy1(read);
}

class Coverager{
	public:
	
	void pileup(GArray& range, const bam1_t* read, int start, int end){
		//check if the read really overlaps
		if (start < range.end() && end >= range.loc){
			if (range.strand >= 0){
				//range on the reference strand
				int pos = start-range.loc;
				++range.array[pos>0?pos:0];
				pos = end + 1 - range.loc;
				if (pos < range.len){
					--range.array[pos];
				}
			} else {
				//range on the negative strand
				int pos = range.end() - 1 - end;
				++range.array[pos>0?pos:0];
				pos = range.end() - start;
				if (pos < range.len){
					--range.array[pos];
				}
			}
		}
	}
};

class Pileupper{
	public:
	
	int binsize;
	int shift;
	bool ss;
	
	Pileupper(int abinsize, int ashift, bool ass){
		binsize = abinsize;
		shift = ashift;
		ss = ass;
	}
	
	//can be accelerated by rewriting a similar one for special cases:
	//ss=false, ignoreRegionStrand, binsize=1
	//and storing in range extra information: range.loc - shift and range.loc + shift
	//this could be probably done by templating this class (so everything would be written just once).
	void pileup(GArray& range, const bam1_t* read, int start, int end){
		bool negstrand = isNegStrand(read);
		//relative position from the start of the range
		int pos = negstrand?end - range.loc - shift:start + shift - range.loc;
		//check overlap with the region
		if (pos < 0 || pos >= range.len) return;
		//strand of the region defines direction and sense and antisense
		if (range.strand < 0){
			pos = range.len - pos - 1;
			negstrand = !negstrand;
		}
		//even positions of the array are sense-reads, odd are antisense
		if (ss){ ++range.array[2*(pos/binsize) + (negstrand?1:0)]; }
		else {++range.array[pos/binsize]; } 
	}
};

// [[Rcpp::export]]
List pileup_core(RObject gr, std::string bampath, int mapqual=0, int binsize=1, int shift=0, bool ss=false, int maxgap=16385){
	std::vector<GArray> ranges;
	//opening bamfile and index
	Bamfile bfile(bampath);
	//adding regions to the vector
	parseRegions(ranges, gr, bfile.in);
	//allocate memory
	List ret = allocateAndDistributeMemory(ranges, binsize, ss);
	//pileup
	Pileupper p(binsize, shift, ss);
	overlapAndPileup(bfile, ranges, mapqual, shift, p, maxgap);
	//close bamfile and index
	bfile.close();
	
	return ret;
}

// [[Rcpp::export]]
List coverage_core(RObject gr, std::string bampath, int mapqual=0, int maxgap=16385){
	std::vector<GArray> ranges;
	//opening bamfile and index
	Bamfile bfile(bampath);
	//adding regions to the vector
	parseRegions(ranges, gr, bfile.in);
	//allocate memory
	List ret = allocateAndDistributeMemory(ranges, 1, false);
	//pileup
	Coverager c;
	overlapAndPileup(bfile, ranges, mapqual, 0, c, maxgap);
	//close bamfile and index
	bfile.close();
	//do cumsum on all the ranges
	typedef std::vector<GArray>::iterator IGArray;
	IGArray i_ranges = ranges.begin(); IGArray e_ranges = ranges.end(); 
	for (; i_ranges < e_ranges; ++i_ranges){
		int acc = i_ranges->array[0];
		int* e_array = i_ranges->array + i_ranges->len;
		for (int* i_array = i_ranges->array + 1; i_array < e_array; ++i_array){
			*i_array = acc += *i_array;
		}
	}
	
	return ret;
}

void loop(std::string bampath){
	//opening bamfile and index
	Bamfile bfile(bampath);
	
	int sum = 0;
	bam1_t* read = bam_init1();
	
	while (bam_read1((bfile.in)->x.bam, read) >= 0){ ++sum; }
	Rcout << "There are " << sum << " total reads in the bam file" << std::endl;
	
	bam_destroy1(read);
	bfile.close();
}