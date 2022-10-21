# TopDeck
TopDeck is a novel method for identifying fusion gene samples from The Cancer Genome Atlas (TCGA) using recount3. Developed to create a subset of potential fusion gene samples that can then be run on more accurate fusion finders, TopDeck combines a gene overexpression analysis with a novel exon comparison analysis to identify a list of potential fusion genes. 

TopDeck allows for the analysis of a single 3' fusion gene partner, by identifying if it is overexpressed, and if there is a significant increase in exon expression after the expected fusion gene breakpoint. Developed on 10 genes which act as common 3' fusion genes in TCGA (ALK, ARHGAP26, ERG, ETV1, MAML3, NTRK3, RARA, RET, TACC3, TFE3), TopDeck had an overall sensitivity of 36.83%.

TopDeck is the culmination of an Honours thesis for Monash University in collaboration with the Peter MacCallum Cancer Centre.
