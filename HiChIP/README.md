# Enhancer connectome in primary human cells identifies target genes of disease-associated DNA elements

# HiChIP data processing. 

HiChIP paired-end reads were aligned to the hg19
or mm9 genome using the HiC-Pro pipeline51. Default settings were used to
remove duplicate reads, assign reads to MboI restriction fragments, filter for
valid interactions, and generate binned interaction matrices. HiC-Pro filtered
reads were then processed into a .hic file using the hicpro2juicebox function.
The Juicer pipeline HiCCUPS tool was used to identify high-confidence
loops4 using the same parameters as for the GM12878 in situ Hi-C map: hiccups
-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 .hic_input
HiCCUPS_output. For T cell Juicer loops, performing the default Juicer calls
resulted in a high rate of false positives upon visual inspection of the interaction
matrix. We therefore called loops with the same HiCCUPS parameters
in two biological replicates for each T cell subtype and then filtered loops for
those that were reproducibly called in both replicates. In addition, we removed
all loops greater than 1 Mb in length.
1D signal enrichment and peak calling were generated from the HiC-Pro
filtered contacts file. Intrachromosomal contacts were filtered, and both
anchors were extended by 75 bp. The combined bed file containing both
anchors was then used to generate bigwigs for visualization in the WashU
Epigenome Browser or call peaks using MACS2.
Allele-specific HiChIP data processing was achieved using HiC-Pro’s allelespecific
analysis features51. First, HCASMC phasing data41 were used to mask
the hg19 genome and make indexes. HiC-Pro settings were similar to those
described above, with the exception that reads were aligned to the masked
genome and then assigned to a specific allele on the basis of phasing data.


51. Servant, N. et al. HiC-Pro: an optimized and flexible pipeline for Hi-C data processing. Genome Biol. 16, 259 (2015).

4. Rao, S.S.P. et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159, 1665–1680 (2014)

41. Miller, C.L. et al. Integrative functional genomics identifies regulatory mechanisms at coronary artery disease loci. Nat. Commun. 7, 12092 (2016).
