# Breakpoint refinement pipeline for Oxford Nanopore-discovered structural variants

This repository contains a pipeline used to update the breakpoint location and sequence content (for insertions) for structural variants (SVs) called from Oxford Nanopore data.
This pipeline has only been tested using output from the [Sniffles](https://github.com/fritzsedlazeck/Sniffles) SV caller run on Oxford Nanopore data.
We may add support for other SV callers or sequencing platforms if there is interest to do so from the community.

This pipeline is mostly written R, but uses the `system` function to pass parameters to a bash script that calls the external software that does the bulk of the data processing.
The pipeline is not guaranteed to run on a non-Linux system or on a shell other than bash.

# Installation

The scripts can be downloaded from the repository by running the following command:

	git clone https://github.com/malemay/breakpoint_refinement.git

The pipeline has not been formalized into a proper `R` package, but might be in the future if there is community interest.

# Software dependencies

The following programs should be installed for the pipeline to run.
The versions we used for development are indicated in parentheses.

* [AGE](https://github.com/abyzovlab/AGE) (commit 6fa60999f573998a95b6ef751b454e6719b1849d)
* [minimap2](https://github.com/lh3/minimap2) (commit c9874e2dc50e32bbff4ded01cf5ec0e9be0a53dd)
* [`R` programming language](https://cran.r-project.org/) (version 4.0.3)
* [samtools](https://github.com/samtools/samtools) (commit 26d7c73c690d298c3d4f6979224933d2e2d102cf)
* [wtdbg2](https://github.com/ruanjue/wtdbg2) (commit 79334f4d92084f5f5ff81b48f6b993ae14b1d88d)

The `parallel` package also needs to be installed for the pipeline to be able to run in parallel on multiple cores.
It can be installed by running this command in `R`:

	install.packages("parallel")

# Arguments to `refine_breakpoints`

`refine_breakpoints` is the only function that most users should really need to run the pipeline.
It takes as input a `.vcf` file containing SV calls originating from a single sample as well as a `.bam` file containing the Oxford Nanopore reads for that sample.
Its arguments are documented below.

* `input_vcf`: name of the input VCF file (character vector of length 1) 
* `output_vcf`: name of the output VCF file (character vector of length 1)
* `ncores`: number of cores for the for parallel processing (default = 1)
* `reference_window`: number of bases to span on either side of SV position to extract from reference genome for alignment (default = 500)
* `reads_window`: number of bases to span on either side of SV position to extract reads from the `.bam` file for alignment (default = 200)
* `min_overlap`: minimum relative overlap between original and proposed refined deletion for the SV to be updated in the output VCF (default = 0.5)
* `min_identity`: minimum percent identity between aligned sequences flanking the SV for the SV to be updated in the output VCF (default = 85) 
* `max_gaps`: maximum percent gaps between aligned sequences flanking the SV for the SV to be updated in the output VCF (default = 15)
* `max_distance`: maximum relative (Levensthein) edit distance between the original insertion sequence and proposed refined sequence for the SV to be updated in the output VCF (default = 0.5)
* `max_offset`: maximum distance (in bases) between the original insertion position and proposed refined position for the SV to be updated in the output VCF (default = 50)
* `max_svlen`: maximum length (in bases) of the SV for consideration by the pipeline; larger values considerable increase memory and computing time requirements (default = 5000)
* `age_script`: path to the shell script that will be used to launch the external programs (default = "age_realign.sh")
* `samtools`: path to the samtools executable (default = "samtools")
* `minimap2`: path to the minimap2 executable (default = "minimap2")
* `age`: path to the age_align executable (default = "age_align")
* `wtdbg2`: path to the wtdbg2 executable (default = "wtdbg2")
* `wtpoa_cns`: path to the wtpoa-cns executable (default = "wtpoa-cns")
* `refgenome`: path to the reference genome
* `nanopore_bam`: path to the `.bam` file containing the reads used for alignment and polishing (should match the sample from which the VCF file originates)


# Testing the installation

Sourcing the file `breakpoint_refinement.R` in R gives access to the `refine_breakpoints` function and other functions it needs.

	source("breakpoint_refinement.R")

The following command can be run in R to test the installation on test data, provided that the programs listed above are in your `$PATH`:

	refine_breakpoints("input_test.vcf", "output_test.vcf", "test.bam", "refgenome.fa", ncores = 4)


