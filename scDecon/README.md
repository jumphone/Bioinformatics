https://cibersort.stanford.edu/runcibersort.php

Deconvolution analysis
CIBERSORT was used to perform the deconvolution analysis of the bulk and scRNA-seq tumour data against the mouse clusters52. The full transcriptomes of the tumour data were used as the input mixture, and the signature gene input was the mouse cluster expression matrix after removal of all genes that could introduce bias in the deconvolution process, including 1,400 cell cycle genes, 300 genes associated with ribosome biogenesis, and around 100 mitochondrial and apoptosis-related genes. Quantile normalization was disabled and 100–500 permutations were run. To test CIBERSORT on our datasets, we created synthetic bulk mixtures from the mouse clusters, and selected known amounts of reads from various clusters. CIBERSORT roughly yielded the expected relative abundances. In order to generate reliable input expression profiles, tumour clusters with a very low number of cells were discarded from the analysis. To validate our mouse cluster signatures, we obtained human and mouse data of brain cell types from published datasets33,34 (Supplementary Table 1, sheet 2), and deconvoluted them against our mouse signatures to ensure that our expected abundances had similar values to our cell of origin matches.

