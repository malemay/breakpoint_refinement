# Breakpoint refinement pipeline for Oxford Nanopore-discovered structural variants

This repository contains a pipeline used to update the breakpoint location and sequence content (for insertions) for structural variants (SVs) called from Oxford Nanopore data.
This pipeline has only been tested using output from the [Sniffles](https://github.com/fritzsedlazeck/Sniffles) SV caller.
We may add support for other SV callers if there is interest to do so from the community.

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

`refine_breakpoints` is the only function that most user should really need to use to run the pipeline.
Its arguments are documented below.

input_vcf
output_vcf
ncores = 1
reference_window
reads_window,
min_overlap, min_identity, max_gaps,
max_distance, max_offset, max_svlen,
age_script, samtools, minimap2,
age, wtdbg2, wtpoa_cns, refgenome,
nanopore_bam) {


# Testing the installation

