
### ------------------------------------------------------------

# This function accepts arguments from update_breakpoints() in order
# to launch the SV realignment pipeline which involves the following steps:

# - Assembly of the Nanopore reads from the region of interest with wtdbg2
# - Polishing of the resulting assembly using the Nanopore reads themselvs
#   aligned to the assembled contig with minimap2
# - Alignement of the resulting contig to the reference genome using AGE

# This function is merely a wrapper around the bash script age_realign.sh

# The function creates some temporary files that are used by the bash script.
# It also manages those files by removing those that are no longer used at the
# end and returning the names of the alignment file (.age.txt) and contig fasta
# file (.ontpolish.fa) to the function that calls it so it knows where to find them

# This function also creates the command line that will be parsed using getopts
# by the bash script.

call_age <- function(chr, svtype, start, svlen, reference_window, reads_window,
		     age_script, samtools, minimap2, age, wtdbg2, wtpoa_cns,
		     refgenome, nanopore_bam) {

	# Creating temporary files that will be used by the bash script

	# The file to which the Nanopore reads will be written
	reads_fasta <- tempfile("reads_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(reads_fasta)
	on.exit(file.remove(reads_fasta), add = TRUE)

	# The file to which the assembled contig will be written
	contig_fasta <- tempfile("contig_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(contig_fasta)
	on.exit(file.remove(contig_fasta), add = TRUE)

	# The file to which the polished contig will be written
	polished_fasta <- tempfile("polished_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(polished_fasta)

	# The file to which the reference sequence in the region will be written
	reference_fasta <- tempfile("reference_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(reference_fasta)
	on.exit(file.remove(reference_fasta), add = TRUE)

	# The file to which the age alignments will be written
	age_file <- tempfile("age", tmpdir = getwd(), fileext = ".txt")
	file.create(age_file)

	# The prefix that fill be used
	prefix <- tempfile("prefix", tmpdir = "")
	prefix <- sub("/", "", prefix)

	command <- paste0(age_script,
			  " -s ", samtools,
			  " -m ", minimap2,
			  " -a ", age,
			  " -w ", wtdbg2,
			  " -c ", wtpoa_cns,
			  " -t ", svtype,
			  " -h ", chr,
			  " -e ", start,
			  " -l ", svlen,
			  " -r ", refgenome,
			  " -b ", nanopore_bam,
			  " -n ", reads_fasta,
			  " -f ", contig_fasta,
			  " -p ", polished_fasta,
			  " -z ", reference_fasta,
			  " -x ", age_file,
			  " -y ", prefix,
			  " -u ", reads_window,
			  " -d ", reference_window)

	# Launching the command itself
	system(command)

	# Returning the names of the files to be used by update_breakpoints()
	return(c(age_file, polished_fasta))
}

### ------------------------------------------------------------

# This function takes the output of a call to parse_svinfo and
# that of a call to parse_age for the same SV and formats
# the data to be used for comparing the two ranges of coordinates
# into a data.frame

gather_align_data <- function(svinfo, age, window_range, contig_fasta) {

  # We need to extract the data from age in case more than one alignment was found
  excise_ref_len <- sapply(age$excise, function(x) x$ref$len)
  excise_ref_start <- sapply(age$excise, function(x) x$ref$start)
  excise_ref_end <- sapply(age$excise, function(x) x$ref$end)

  excise_contig_len <- sapply(age$excise, function(x) x$contig$len)
  excise_contig_start <- sapply(age$excise, function(x) x$contig$start)
  excise_contig_end <- sapply(age$excise, function(x) x$contig$end)
  
  # Checking that the input is suitable
  if(length(excise_ref_len)      < 1) return(NULL)
  if(length(excise_ref_start)    < 1) return(NULL)
  if(length(excise_ref_end)      < 1) return(NULL)
  if(length(excise_contig_len)   < 1) return(NULL)
  if(length(excise_contig_start) < 1) return(NULL)
  if(length(excise_contig_end)   < 1) return(NULL)

  output <- data.frame(chr = svinfo$chr,
		       svtype = svinfo$svtype,
  		       start = svinfo$start,
		       svlen = svinfo$svlen,
		       alt = svinfo$alt,
		       imprecise = svinfo$imprecise,
		       age_refname = age$reference$name,
		       age_refstart = age$reference$range[1],
		       age_refend = age$reference$range[2],
		       age_reflen = age$reference$length,
		       age_contigname = age$contig$name,
		       age_contigstart = age$contig$range[1],
		       age_contigend = age$contig$range[2],
		       age_contiglen = age$contig$length,
		       score = age$score,
		       identity = age$identity,
		       gaps = age$gaps,
		       time = age$time,
		       excise_ref_len = excise_ref_len,
		       excise_ref_start = excise_ref_start,
		       excise_ref_end = excise_ref_end,
		       excise_contig_len = excise_contig_len,
		       excise_contig_start = excise_contig_start,
		       excise_contig_end = excise_contig_end,
		       stringsAsFactors = FALSE)

  # Adding a series of variables to make the SVs truly comparable
  output$age_start  <- output$excise_ref_start - window_range + output$start
  output$age_length <- abs(excise_contig_len - excise_ref_len)

  # Getting the alternate sequence of insertions from the contig fasta
  if(all(output$svtype == "INS")) {
    ins_range <- GenomicRanges::GRanges(seqnames = output$age_contigname,
					IRanges::IRanges(start = pmin(output$excise_contig_start, 
							     output$excise_contig_end),
						end = pmax(output$excise_contig_start,
							   output$excise_contig_end)))
    Rsamtools::indexFa(contig_fasta)
    output$age_alt <- as.character(Rsamtools::scanFa(contig_fasta, ins_range))

    # We have to set it to reverse complement if it was aligned in reverse
    for(i in 1:nrow(output)) {
      if(output[i, "age_contigstart"] > output[i, "age_contigend"]) {
	output[i, "age_alt"] <- revcomp(output[i, "age_alt"])
      }

    # We also compute the Levenshtein distance
    output[i, "lev_dist"] <- adist(output[i, "alt"], output[i, "age_alt"])
    }

  # Then we compute the distance ratio as the Levensthein distance divided by the largest of the two sequences
    output$dist_ratio <- output$lev_dist / pmax(output$svlen, output$excise_contig_len)

    # We also determine the distance between the two insertion points
    output$offset <- abs(output$start - output$age_start)

  } else {
    output$age_alt <- NA
    output$lev_dist <- NA
    output$dist_ratio <- NA
    output$offset <- NA
  }

  # For deletions we want to compute the reciprocal overlap between the two deletions
  if(all(output$svtype == "DEL")) {
	  # Creating IRanges objects for each of the two sets
	  sniffles_ranges <- IRanges::IRanges(start = output$start, width = output$svlen)
	  age_ranges      <- IRanges::IRanges(start = output$age_start, width = output$excise_ref_len)

	  # Computing the overlapping ranges between the two
	  ol_ranges <- IRanges::pintersect(sniffles_ranges, age_ranges, resolve.empty = "start.x")
	  # Computing a column for the reciprocal overlap (minimum of the two overlaps)
	  output$r_overlap <- pmin(IRanges::width(ol_ranges) / IRanges::width(sniffles_ranges), 
				   IRanges::width(ol_ranges) / IRanges::width(age_ranges))
	  
	  # If there is no deletion, we set the overlap to 0 (otherwise it is NaN)
	  output$r_overlap[IRanges::width(age_ranges) == 0] <- 0

  } else {
	  output$r_overlap <- NA
  }
  
  return(output)
}

### ------------------------------------------------------------

# A function used to parse the output of the AGE aligner

# It takes the name of the file to be parsed (a string) as input and returns
# a list containing useful information extracted from the file.

parse_age <- function(filename) {
  # Reading the whole file in memory using scan (file is overall small)
  age <- scan(filename, what = character(), sep = "\n")

  # Counting the number of alignments in the file
  # This is done by counting the number of matches to the regexp "^MATCH"
  block_start <- grep("^MATCH", age)
  n_align <- length(block_start)

  # Checking that any line matched
  if(n_align == 0) {
	  warning("No alignment returned for ", filename, ", returning NULL")
	  return(NULL)
  }


  # Otherwise checking the number of alignments
  # And selecting the first one if more than one is found
  if(n_align > 1) {

	  warning(n_align, " contigs aligned for ", filename, ", using only first one")

	  # Extracting the first part of the alignment (region between the two first ^MATCH)
	  age <- age[block_start[1]:(block_start[2] - 1)]
  }

  # Checking that there has been an excised fragment
  if(!any(grepl("^EXCISED", age))) {
	  warning("No excised fragment for ", filename, ", returning NULL")
	  return(NULL)
  }

  # Extracting the lines pertaining to the reference and assembled contig
  ref    <- grep("^First  seq", age, value = TRUE)
  contig <- grep("^Second seq", age, value = TRUE)
  stopifnot(length(ref) == length(contig))


  # Extracting some info from those lines
  # Sequence name (from fasta) is truncated at blank so it can be parsed by fasta
  ref_start <- as.numeric(sub("^.*\\[\\s*(\\d+),.*", "\\1", ref))
  ref_end   <- as.numeric(sub("^.*\\[\\s*\\d+,\\s*(\\d+)\\].*", "\\1", ref))
  ref_len   <- as.numeric(sub("^.*=>\\s*(\\d+)\\s+nucs.*", "\\1", ref))
  ref_name  <- sub(".*'([[:graph:]]+).*'", "\\1", ref)

  contig_start <- as.numeric(sub("^.*\\[\\s*(\\d+),.*", "\\1", contig))
  contig_end   <- as.numeric(sub("^.*\\[\\s*\\d+,\\s*(\\d+)\\].*", "\\1", contig))
  contig_len   <- as.numeric(sub("^.*=>\\s*(\\d+)\\s+nucs.*", "\\1", contig))
  contig_name  <- sub(".*'([[:graph:]]+).*'", "\\1", contig)

  # Extracting some metadata
  score <- grep("^Score", age, value = TRUE)
  stopifnot(length(score) == 1)
  score <- as.numeric(regmatches(score, regexpr("\\d+", score)))

  identity <- grep("^Identic", age, value = TRUE)
  stopifnot(length(identity) == 1)
  identity <- as.numeric(sub("^.*\\(\\s*(\\d+)%\\)\\s*nucs.*", "\\1", identity))

  gaps <- grep("^Gaps", age, value = TRUE)
  stopifnot(length(gaps) == 1)
  gaps <- as.numeric(sub("^.*\\(\\s*(\\d+)%\\)\\s*nucs.*", "\\1", gaps))

  time <- grep("^Alignment time is", age, value = TRUE)
  stopifnot(length(time) == 1)
  time <- as.numeric(regmatches(time, regexpr("\\d+\\.?\\d*", time)))

  # We have to delimit the lines in which to search for excised regions
  ex_range <- c(grep("^EXCISED REGION\\(S\\):$", age), grep("^Identity at breakpoints", age))
  stopifnot(length(ex_range) == 2 && diff(ex_range) >= 1)
  ex_lines <- age[ex_range[1]:ex_range[2]]

  # Now we initialize a list() that will contain the info on excised regions
  excised <- list()

  # We extract the lines containing the excision info
  ref_excised    <- grep("^ first  seq", ex_lines, value = TRUE)
  contig_excised <- grep("^ second seq", ex_lines, value = TRUE)
  stopifnot(length(ref_excised) >= 1 && length(ref_excised) == length(contig_excised))

  for(i in 1:length(ref_excised)) {
    # Initializing the sub-element of the list
    excised[[i]] <- list()
    excised[[i]]$ref$len <- as.numeric(sub("^.*\\s+(\\d+)\\s+nucs.*", "\\1", ref_excised[i]))
    excised[[i]]$ref$start <- as.numeric(sub("^.*nucs\\s+\\[\\s*(\\d+),.*", "\\1", ref_excised[i]))
    excised[[i]]$ref$end <- as.numeric(sub("^.*\\[.*,(\\d+)\\s*\\].*", "\\1", ref_excised[i]))

    excised[[i]]$contig$len <- as.numeric(sub("^.*\\s+(\\d+)\\s+nucs.*", "\\1", contig_excised[i]))
    excised[[i]]$contig$start <- as.numeric(sub("^.*nucs\\s+\\[\\s*(\\d+),.*", "\\1", contig_excised[i]))
    excised[[i]]$contig$end <- as.numeric(sub("^.*\\[.*,(\\d+)\\s*\\].*", "\\1", contig_excised[i]))
  }

  # Preparing the list object to be returned
  output <- list(reference = list(name = ref_name, range = c(ref_start, ref_end), length = ref_len),
		 contig = list(name = contig_name, range = c(contig_start, contig_end), length = contig_len),
		 score = score,
		 identity = identity,
		 gaps = gaps,
		 time = time,
		 excised = excised)

  return(output)
}

### ------------------------------------------------------------

# A function that parses a line in a VCF file in order to extract metadata
# about the structural variants it represents

parse_svinfo <- function(vcf_line) {

  stopifnot(length(vcf_line) == 1)

  # Splitting the vcf line into its fields
  split_line <- strsplit(vcf_line, "\\s+")[[1]]

  # Assigning some elements directly
  chr <- split_line[1]
  start <- as.numeric(split_line[2])
  alt <- split_line[5]

  # Some require a bit of processing
  svtype <- sub(".*SVTYPE=([[:upper:]]+).*", "\\1", split_line[8])
  svlen  <- as.numeric(sub(".*SVLEN=-?(\\d+).*", "\\1", split_line[8]))
  imprecise <- grepl("IMPRECISE", split_line[8])

  list(chr = chr, svtype = svtype, start = start, svlen = svlen, alt = alt, imprecise = imprecise)
}

### ------------------------------------------------------------

# This function takes the name of a vcf file as input and processes 
# the input line by line to write the output back to another vcf
# file.

# Many other functions work under the hood to allow this function
# to process the vcf lines in parallel and call the external commands
# that produce the alignment:

# - update_breakpoints()

refine_breakpoints <- function(input_vcf, output_vcf, nanopore_bam, refgenome,
			       ncores = 1, reference_window = 500, reads_window = 200,
			       min_overlap = 0.5, min_identity = 85, max_gaps = 15,
			       max_distance = 0.5, max_offset = 50, max_svlen = 50000,
			       age_script = "./age_realign.sh", samtools = "samtools", 
			       minimap2 = "minimap2", age = "age_align", wtdbg2 = "wtdbg2",
			       wtpoa_cns = "wtpoa-cns") {

	# Reading the input vcf_file using the scan function
	vcf_lines <- scan(input_vcf, what = character(), sep = "\n", quiet = TRUE)

	# Extracting the header lines and the body lines of the vcf
	header_lines <- grep("^#", vcf_lines, value = TRUE)
	body_lines <- grep("^#", vcf_lines, value = TRUE, invert = TRUE)

	# A sanity check to make sure that the splitting worked properly
	stopifnot(length(vcf_lines) == length(header_lines) + length(body_lines))

	# Finding the last ##INFO line
	last_info <- max(grep("^##INFO", header_lines))

	# Creating a character vector with the lines to add
	info_lines <- c(
	'##INFO=<ID=MAX_SVLEN,Number=0,Type=Flag,Description="Breakpoints not refined because SVLEN exceeds the threshold set">',
	'##INFO=<ID=NO_ALIGNMENT,Number=0,Type=Flag,Description="Breakpoints not refined because no alignement was produced">',
	'##INFO=<ID=GATHER_DATA_FAILED,Number=0,Type=Flag,Description="Breakpoints not refined because gathering alignment failed">',
	'##INFO=<ID=BELOW_THRESHOLDS,Number=0,Type=Flag,Description="Breakpoints not refined because quality thresholds not met">',
	'##INFO=<ID=REALIGNED,Number=0,Type=Flag,Description="Breakpoints refined by the refine_breakpoints() R function">',
	'##INFO=<ID=ORIGINAL_START,Number=1,Type=Integer,Description="Start position of variant before breakpoints were refined">',
	'##INFO=<ID=ORIGINAL_LEN,Number=1,Type=Integer,Description="Length of variant before breakpoints were refined">',
	'##INFO=<ID=ORIGINAL_ALT,Number=1,Type=String,Description="ALT allele of the variant before breakpoints were refined">',
	'##INFO=<ID=CONTIG_LEN,Number=1,Type=Integer,Description="Length of the contig used for breakpoint refinement">',
	'##INFO=<ID=IDENTITY,Number=1,Type=Integer,Description="Percent identity of the alignment used to refine the breakpoints">',
	'##INFO=<ID=GAPS,Number=1,Type=Integer,Description="Percentage of gaps in the alignment used to refine the breakpoints">',
	'##INFO=<ID=OVERLAP,Number=1,Type=Float,Description="Reciprocal overlap between original and breakpoint-refined variant">',
	'##INFO=<ID=DISTANCE_RATIO,Number=1,Type=Float,Description="Edit distance ratio between original and refined allele">',
	'##INFO=<ID=OFFSET,Number=1,Type=Integer,Description="Absolute difference (bp) between original allele and refined allele">'
	)

	# Updating the vcf_lines with those metainfo lines
	header_lines <- c(header_lines[1:last_info], info_lines, header_lines[(last_info + 1):length(header_lines)])

	# We also write the header to the output file
	output_file <- file(output_vcf, open = "w")
	on.exit(close(output_file), add = TRUE)
	cat(header_lines, file = output_file, sep = "\n")

	# We coerce the lines of the vcf body to a list for processing with parallel::mclapply()
	body_lines <- as.list(body_lines)

	# Processing with mclapply()
	realigned_lines <- parallel::mclapply(body_lines, FUN = update_breakpoints, 
					      mc.cores = ncores,
					      reference_window = reference_window,
					      reads_window = reads_window,
					      min_overlap = min_overlap,
					      min_identity = min_identity,
					      max_gaps = max_gaps,
					      max_distance = max_distance,
					      max_offset = max_offset,
					      max_svlen = max_svlen,
					      age_script = age_script,
					      samtools = samtools,
					      minimap2 = minimap2,
					      age = age,
					      wtdbg2 = wtdbg2,
					      wtpoa_cns = wtpoa_cns,
					      refgenome = refgenome,
					      nanopore_bam = nanopore_bam)

	# Turning the list back into a character vector
	realigned_lines <- as.character(realigned_lines)

	# Writing the results to the output file
	cat(realigned_lines, file = output_file, sep = "\n")

	# No useful value is returned because the output is written to file
	return(invisible(NULL))
}

### ------------------------------------------------------------

# A function that computes the reverse-complement of the input string

revcomp <- function(sequence) {

	# Checking that only one sequence is provided
	stopifnot(length(sequence) == 1)

	# The lookup table that will be used for replacement
	rep_table <- c("A" = "T",
		       "T" = "A",
		       "G" = "C",
		       "C" = "G",
		       "N" = "N")

	# Splitting the sequence into its constituent nucleotides
	sequence <- strsplit(sequence, "")[[1]]

	# Replacing the nucleotides by their complement
	sequence <- rep_table[sequence]

	# Returning the inverted sequence
	paste0(rev(sequence), collapse = "")
}

### ------------------------------------------------------------

# This function accepts a vcf line of a structural variant as input
# and returns a line with the vcf coordinates optionally updated following
# realignment with AGE.

# The function is meant to be called via the refine_breakpoints() function
# which will apply it in parallel to every line of a vcf file using mclapply()

# This function requires the following functions to be available in the globalenv()
# - parse_svinfo() to parse the informations from the vcf line
# - call_age() which is a wrapper around the bash script age_realign.sh
# - parse_age() which gets in the information from the age file into R
# - gather_align_data() which makes a data.frame from the age and svinfo data

update_breakpoints <- function(vcf_line, reference_window, reads_window,
			       min_overlap, min_identity, max_gaps,
			       max_distance, max_offset, max_svlen,
			       age_script, samtools, minimap2, age,
			       wtdbg2, wtpoa_cns, refgenome, nanopore_bam) {

	# Parsing the information from the vcf line
	svinfo <- parse_svinfo(vcf_line)

	# We do not realign duplications and inversions; these are returned unmodified
	if(svinfo$svtype %in%  c("DUP", "INV")) {
		return(vcf_line)
	}

	# At this point there should only be insertions and deletions left
	stopifnot(svinfo$svtype %in% c("DEL", "INS"))

	# We do not process this SV if if is larger than max_svlen
	if(svinfo$svlen > max_svlen) {
		warning("Following line not processed, SV length > ", max_svlen, "\n", vcf_line)
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";MAX_SVLEN")
		return(paste0(sline, collapse = "\t"))
	}

	# Creating two boolean values to make the code clearer
	del <- svinfo$svtype == "DEL"
	ins <- svinfo$svtype == "INS"

	# Passing the information to the call_age() function
	# This function returns a character vector of length 2 with
	# the information needed for the rest of the function :
	# - [1] the name of the .age.txt file created
	# - [2] the name of the .ontpolish.fa file containing the assembled and polished contig
	age_files <- call_age(chr = svinfo$chr, svtype = svinfo$svtype, start = svinfo$start,
			      svlen = svinfo$svlen, reference_window = reference_window, 
			      reads_window = reads_window, age_script = age_script,
			      samtools = samtools, minimap2 = minimap2, age = age,
			      wtdbg2 = wtdbg2, wtpoa_cns = wtpoa_cns, refgenome = refgenome,
			      nanopore_bam = nanopore_bam)

	# Parsing the age file
	age_data <- parse_age(age_files[1])

	# We check if age returned anything (will be NULL in that case)
	if(is.null(age_data)) {
		# In that case we return the vcf line (almost) unmodified and remove files
		file.remove(age_files[1], age_files[2])
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";NO_ALIGNMENT")
		return(paste0(sline, collapse = "\t"))
	}

	# Otherwise we generate a data.frame with the data pertaining to that line
	align_data <- gather_align_data(svinfo, age_data, reference_window, age_files[2])

	# We no longer need the age and contig files
	file.remove(age_files[1], age_files[2])

	# We also remove the .fai index if it exists
	if(file.exists(fai <- paste0(age_files[2], ".fai"))) file.remove(fai)

	# We check if the call to gather_align_data returned anything useful
	if(is.null(align_data)) {
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";GATHER_DATA_FAILED")
		return(paste0(sline, collapse = "\t"))
	}

	# Based on these data we decide how to realign the variant
	# There may be more than one row in the data.frame so we have to loop over
	# all of them ; we stop looping once we found a suitable one

	# This variable determines whether there has been realignment
	realigned <- FALSE

	for(i in 1:nrow(align_data)) {

		# Extract some data which is the same no matter if INS or DEL
		identity <- align_data[i, "identity"]
		gaps <- align_data[i, "gaps"]

		# We apply different filters if they are deletions or insertions
		if(del) {
			# We extract the reciprocal overlap
			rol <- align_data[i, "r_overlap"]

			if(rol >= min_overlap && identity >= min_identity && gaps <= max_gaps) {
				realigned <- TRUE
				irow <- i
				break
			}

		} else if (ins) {
			# Extracting the Levenshtein distance ratio and absolute offset
			dist_ratio <- align_data[i, "dist_ratio"]
			offset <- abs(align_data[i, "offset"])

			if(dist_ratio <= max_distance && offset <= max_offset && identity >= min_identity && gaps <= max_gaps) {
				realigned <- TRUE
				irow <- i
				break
			}
		}
	}

	# We return the original line if it is not to be updated
	if(!realigned) {
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";BELOW_THRESHOLDS")
		return(paste0(sline, collapse = "\t"))
	}

	# Otherwise we need to update the vcf line
	# We first split the vcf line into its fields
	sline <- strsplit(vcf_line, "\\s+")[[1]]

	# We first update the start position, SVLEN, and ALT allele (if an insertion)
	sline[2] <- align_data[irow, "age_start"]
	if(ins) sline[5] <- align_data[irow, "age_alt"]
	# REF and ALT allele for realigned deletions are changed to N and <DEL> by default
	if(del) {sline[4] <- "N"; sline[5] <- "<DEL>"}

	# The SVLEN modification will depend on whether it is a deletion or an insertion
	if(del) {
		sline[8] <- sub("SVLEN=-\\d+",
				paste0("SVLEN=", as.character(-align_data[irow, "excise_ref_len"])),
				sline[8])
	} else if(ins) {
		sline[8] <- sub("SVLEN=\\d+",
				paste0("SVLEN=", as.character(align_data[irow, "excise_contig_len"])),
				sline[8])
	}

	# We also need to update INFO/END, which depends on whether it is INS or DEL
	# The way END is computed for deletions is not how I think about it
	# (I would have done start + length - 2)
	# But it is how Sniffles does, so I do the same out of consistency
	if(del) {
		sline[8] <- sub("([^[:alpha:]])END=\\d+",
				paste0("\\1END=", as.character(align_data[irow, "age_start"] + 
							    align_data[irow, "excise_ref_len"])),
				sline[8])
	} else if(ins) {
		sline[8] <- sub("([^[:alpha:]])END=\\d+",
				paste0("\\1END=", as.character(align_data[irow, "age_start"])),
				sline[8])
	}

	# Now we can add tags to the INFO field
	# The INFO field is the 8th one so we must update index 8
	sline[8] <- paste0(sline[8], ";REALIGNED")
	sline[8] <- paste0(sline[8], ";ORIGINAL_START=", svinfo$start)
	sline[8] <- paste0(sline[8], ";ORIGINAL_LEN=", ifelse(del, -svinfo$svlen, svinfo$svlen))
	if(ins) sline[8] <- paste0(sline[8], ";ORIGINAL_ALT=", svinfo$alt)
	sline[8] <- paste0(sline[8], ";CONTIG_LEN=", align_data[irow, "age_contiglen"])
	sline[8] <- paste0(sline[8], ";IDENTITY=", as.character(identity))
	sline[8] <- paste0(sline[8], ";GAPS=", as.character(gaps))
	if(del) sline[8] <- paste0(sline[8], ";OVERLAP=", as.character(rol))
	if(ins) sline[8] <- paste0(sline[8], ";DISTANCE_RATIO=", as.character(dist_ratio))
	if(ins) sline[8] <- paste0(sline[8], ";OFFSET=", as.character(offset))

	# By now we are able to return the re-assembled line
	return(paste0(sline, collapse = "\t"))
}

