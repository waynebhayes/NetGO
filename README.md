# NetGO
The scripts in this repo---written mostly in AWK and BASH (the Bourne Again Shell)---are designed to evaluate network alignments (for example those produced by https://github.com/waynebhayes/SANA, of which this repo is a submodule) from the standpoint of shared Gene Ontology (GO) terms.  Here's a brief description of some of the files.

- misc.awk is a very large "library" of functions and code, written in AWK, that is designed to be inculded in all AWK programs here.
- hawk: "Hayes AWK", which is just AWK with the above misc.awk included so all its functions can be considered as standard.
- awkcel: a blantant rip-off on the name "Excel" (Copyright Microsoft), it takes "tab-separated files" (by convention with the extension .tsv, like Excel's .csv), in which the top line is a _header_ containing the tab-separated names of each column, where each name _must_ be a valid AWK variable name. Subsequently, each line assigns each variable name from the column values, so that column values can be referred to by their appropriate variable name.
- ebm.sh: implementation of the _Empirical Brown's Method_ for combining (possibly correlated) p-values into one holistic p-value; see Poole (2016, https://academic.oup.com/bioinformatics/article/32/17/i430/2450768).
- GOeval.sh - method for evaluating network alignments using GO terms, based on the paper _Exact p-values for global network alignments via combinatorial analysis of shared GO terms_, (submitted; preprint available at https://arxiv.org/abs/2010.06415)
- NetGO.awk - actually a bash script, implementing an older method described in Hayes & Mamano (2018; https://academic.oup.com/bioinformatics/article/34/8/1345/4708230).

## old comments
NetGO scoring for sets of (sets of) proteins, eg as exist in many-to-many or multiple network alignments

The GNU awk version will usually be the most up-to-date, though the other language versions will be updated periodically to be inline with the awk version.
