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

* [`R` programming language](https://cran.r-project.org/) (we recommend version 4 minimally)
* [AGE](https://github.com/abyzovlab/AGE)
* [samtools](https://github.com/samtools/samtools)

# Arguments to `refine_breakpoints`

# Testing the installation

