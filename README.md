# Breakpoint refinement pipeline for Oxford Nanopore-discovered structural variants

This repository contains a pipeline used to update the breakpoint location and sequence content (for insertions) for structural variants (SVs) called from Oxford Nanopore data.
This pipeline has only been tested using output from the [Sniffles](https://github.com/fritzsedlazeck/Sniffles) SV caller.
We may add support for other SV callers if there is interest to do so from the community.

This pipeline is mostly written R, but uses the `system` function to pass parameters to a bash script that calls the external software that does the bulk of the data processing.
The pipeline is not guaranteed to run on a non-Linux system or on a shell other than bash.

# Installation

The scripts can be downloaded from the repository by running the following command:

	git clone https://github.com/malemay/breakpoint_refinement.git

# Software dependencies

# Arguments to `refine_breakpoints`

# Testing the installation

