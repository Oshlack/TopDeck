# TopDeck
TopDeck is a novel method for identifying fusion gene samples from The Cancer Genome Atlas (TCGA) using recount3. Developed to create a subset of potential fusion gene samples that can then be run on more accurate fusion finders, TopDeck combines a gene overexpression analysis with a novel exon comparison analysis to identify a list of potential fusion genes. 

TopDeck allows for the analysis of a single 3' fusion gene partner, by identifying if it is overexpressed, and if there is a significant increase in exon expression after the expected fusion gene breakpoint. Developed on 10 genes which act as common 3' fusion genes in TCGA (ALK, ARHGAP26, ERG, ETV1, MAML3, NTRK3, RARA, RET, TACC3, TFE3), TopDeck had an overall sensitivity of 36.83%.

TopDeck is the culmination of an Honours thesis for Monash University in collaboration with the Peter MacCallum Cancer Centre.

## Change in Exon Expression
A demonstration of the change in exon expression method using 3' gene TACC3 on bladder urothelial carcinoma samples is available in file tacc3_demonstration.Rmd, with functions used explained.

Two different change in exon expression methods were used, from change in individual exon expression and change in average exon expression.
- The Z-scores for all samples investigated using the change in individual exon expression are available in file: all_diff_objects.Rdata
- The Z-scores for all samples investigated using the change in average exon expression are available in file: prop_cum_all.Rdata

The true positive rail_ids are available in files:
- truepos_diff.csv (individual exon expression)
- truepos_prop_cum.csv (average exon expression)

The false positive rail_ids are available in files:
- falsepos_diff.csv (individual exon expression)
- falsepos_prop_cum.csv (average exon expression)


The change in exon expression method uses functions from exon_expression_functions.R.

## Bibliography
A bibliography for R packages used is available in bibliography.R
