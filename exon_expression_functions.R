library(recount3)
library(snapcount)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(workflowr)
library(patchwork)

#pulling out known fusions
load("~/Documents/recount_challenge_project/data/all_tcga_fusion_list.RData")
#gives tcga_full_fusion_merge


z_calc <- function (df) {
  mean_df <- apply(df, 1, mean)
  sd_df <- apply(df, 1, sd)
  new_df <- (df-mean_df)/sd_df
  rownames(new_df) <- rownames(df)
  new_df
}

#For 3' genes on the negative strand
get_all_sig_z_diff <- function(gene, project_id) {
  df <- raw_exon_counts_neg(gene)
  proj_df <- subset_project_neg(gene, project_id, df)
  kbp_df <- norm_length_neg(gene, proj_df)
  diff_df <- diff_calc(kbp_df)
  diff_auc_df <- norm_auc_neg(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_z
}

#For 3' genes on the positive strand
get_all_sig_z_diff_pos <- function(gene, project_id) {
  df <- raw_exon_counts(gene)
  proj_df <- subset_project(gene, project_id, df)
  kbp_df <- norm_length(gene, proj_df)
  diff_df <- diff_calc(kbp_df)
  diff_auc_df <- norm_auc(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_z
}

#doing it all with one function
#title, xlabel, ylabel, legend
condensed_cum_function_pos <- function(gene, project_id) {
  df <- raw_exon_counts(gene)
  proj_df <- subset_project(gene, project_id, df)
  cumsum_df <- cumsum_calc(gene, proj_df)
  cumsum_auc_df <- norm_auc(gene, project_id, cumsum_df)
  #z_calc(cumsum_auc_df)
  #graph_log(gene, project_id, cumsum_auc_df, title, xlabel, ylabel, legend)
}

#creating a difference graph with one function
difference_cum_pos <- function(gene, project_id, title, xlabel, ylabel, legend) {
  df <- raw_exon_counts(gene)
  proj_df <- subset_project(gene, project_id, df)
  cumsum_df <- cumsum_calc(gene, proj_df)
  diff_df <- diff_calc(cumsum_df)
  diff_auc_df <- norm_auc(gene, project_id, diff_df)
  #graph_notlog(gene, project_id, diff_auc_df, title, xlabel, ylabel, legend)
  z_calc(diff_auc_df)
}

calc_significant <- function(df) {
  df <- df
  calc <- function(n) {
    sum(df[n,] > 3)
  }
  a <- nrow(df)
  final_df <- as.data.frame(sapply(1:a, calc))
  rownames(final_df) <- rownames(df)
  final_df
}


summary_significant <- function(gene, project_id, df) {
  all_sig <- calc_significant(df) #Count total number of significant samples at each exon
  known_z <- subset_known(gene, project_id, df) #Z scores of samples known to contain fusion genes with that 3' gene
  known_sig <- calc_significant(known_z) #Count number of known samples that are significant at each exon
  final_df <- cbind(all_sig, known_sig) #Merge All and Known dataframes
  colnames(final_df) <- c("All", "Known") #Rename columns to communicate
  final_df #df output
}

#i agress that actually should have the sb separate 
set_sb <- function(gene) {
  sb <- QueryBuilder(compilation = "tcga", regions = gene)
  sb <- set_coordinate_modifier(sb, Coordinates$Within)
  sb <- set_row_filters(sb, strand == "+")
  query_exon(sb)
}

#extract raw exon counts for all tcga samples
raw_exon_counts <- function(gene) {
  exon <- set_sb(gene)
  as.data.frame(as.matrix(assays(exon)$counts))
}

#and then we want to subset to project using these functions:
tcga_info <- function(gene) {
  info <- set_sb(gene)
  tcga_info <- cbind(info@colData@rownames, info@colData@listData[["gdc_cases.project.project_id"]], info@colData@listData[["gdc_cases.samples.submitter_id"]], info@colData@listData[["gdc_file_id"]], info@colData@listData[["gdc_cases.samples.sample_type"]])
  tcga_info <- as.data.frame(tcga_info)
  colnames(tcga_info) <- c("rail_id", "project", "tcga_id", "UCSC_id", "sample_type")
  tcga_info$auc <- info@colData@listData[["auc"]]
  tcga_info$tmc <- info@colData@listData[["mapped_read_count"]]
  tcga_info
}


subset_project <- function(gene, project_id, df) {
  info <- tcga_info(gene)
  subset <- filter(info, project == project_id)
  ids <- subset$rail_id
  df[,ids]
}

#to get the metadata
gene_metadata <- function(gene) {
  info <- set_sb(gene)
  as.data.frame(rowRanges(info))
}

#function necessary to set up the df
norm_df <- function(gene, df) {
  new_df <- df
  metadata <- gene_metadata(gene)
  new_df$width <- metadata$width
  new_df
}



noncum_comp_z_score <- function(gene, project_id) {
  df <- raw_exon_counts(gene)
  proj_df <- subset_project(gene, project_id, df)
  kbp_df <- norm_length(gene, proj_df)
  diff_df <- diff_calc(kbp_df)
  diff_auc_df <- norm_auc(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_sig <- calc_significant(all_z)
  known_z <- subset_known(gene, project_id, all_z)
  known_sig <- calc_significant(known_z)
  final_df <- cbind(all_sig, known_sig)
  colnames(final_df) <- c("All", "Known")
  final_df
}


noncum_comp_z_score_neg <- function(gene, project_id) {
  df <- raw_exon_counts_neg(gene)
  proj_df <- subset_project_neg(gene, project_id, df)
  kbp_df <- norm_length_neg(gene, proj_df)
  diff_df <- diff_calc(kbp_df)
  diff_auc_df <- norm_auc_neg(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_sig <- calc_significant(all_z)
  known_z <- subset_known(gene, project_id, all_z)
  known_sig <- calc_significant(known_z)
  final_df <- cbind(all_sig, known_sig)
  colnames(final_df) <- c("All", "Known")
  final_df
}
#calculative cumulative counts, normalised by exon length if that exon is expressed
cumsum_calc <- function(gene, df) {
  a <- ncol(df)
  df <- norm_df(gene, df)
  filt_norm <- function(n) {
    subset_df <- filter(df, df[,n] > 0)
    cumsum_df <- cumsum(subset_df)
    norm_df <- as.data.frame(cumsum_df[,n]/cumsum_df$width*1000)
    rownames(norm_df) <- rownames(cumsum_df)
    colnames(norm_df) <- "rail_id"
    norm_df$rn1 <- as.numeric(rownames(norm_df))
    temp_df <- data.frame(rn1 = as.numeric(rownames(df)))
    new_df <- left_join(temp_df, norm_df)
    final_df <- as.data.frame(as.numeric(new_df[,2]))
    rownames(final_df) <- new_df$rn1
    colnames(final_df) <- colnames(df)[n]
    final_df
  }
  na_df <- as.data.frame(sapply(1:a, filt_norm))
  rownames(na_df) <- rownames(df)
  filled_df <- na_df %>%
    fill(names(na_df))
  filled_df[is.na(filled_df)] <- 0
  filled_df
}

#normalise by auc, subsetted by project
norm_auc <- function(gene, project_id, df) {
  info <- tcga_info(gene)
  info_subset <- filter(info, project == project_id)
  auc <-  info_subset$auc/1000000
  new_df <- as.data.frame(mapply('/', df, auc))
  rownames(new_df) <- rownames(df)
  new_df
}

#get known ids from tcga_full_fusion_merge
get_known_ids <- function(gene, project_id) {
  subset <- filter(tcga_full_fusion_merge, project == project_id)
  subset_further <- filter(subset, `Second Gene` == gene)
  print(subset_further$rail_id)
}

#and then subset known
subset_known <- function(gene, project_id, df) {
  ids <- get_known_ids(gene, project_id)
  final_df <- as.data.frame(df[,ids])
  colnames(final_df) <- ids
  final_df
}

#true positives
get_tp <- function(gene, df) {
  ids <- get_known_ids_new(gene)
  subset(df, rownames(df) %in% ids)
}

#For 3' genes on the negative strand
function(gene, project_id) {
  df <- raw_exon_counts_neg(gene)
  proj_df <- subset_project_neg(gene, project_id, df)
  kbp_df <- norm_length_neg(gene, proj_df)
  diff_df <- diff_calc(kbp_df)
  diff_auc_df <- norm_auc_neg(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_z
}

#For 3' genes on the positive strand
get_all_sig_z_diff_pos <- function(gene, project_id) {
  df <- raw_exon_counts(gene)
  proj_df <- subset_project(gene, project_id, df)
  kbp_df <- norm_length(gene, proj_df)
  diff_df <- diff_calc(kbp_df)
  diff_auc_df <- norm_auc(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_z
}

#For 3' genes on the negative strand
get_all_sig_z_prop_cum <- function(gene, project_id) {
  df <- raw_exon_counts_neg(gene)
  proj_df <- subset_project_neg(gene, project_id, df)
  cumsum_df <- cumsum_calc_neg(gene, proj_df)
  diff_df <- diff_calc(cumsum_df)
  diff_auc_df <- norm_auc_neg(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_z
}

#For 3' genes on the positive strand
get_all_sig_z_prop_cum_pos <- function(gene, project_id) {
  df <- raw_exon_counts(gene)
  proj_df <- subset_project(gene, project_id, df)
  cumsum_df <- cumsum_calc(gene, proj_df)
  diff_df <- diff_calc(cumsum_df)
  diff_auc_df <- norm_auc(gene, project_id, diff_df)
  all_z <- z_calc(diff_auc_df)
  all_z
}

#and then the cumulative graph - log2
graph_log <- function(gene, project_id, df, title, xlabel, ylabel, legend) {
  known_df <- as.data.frame(subset_known(gene, project_id, df))
  all_df <- df
  n <- nrow(df)
  known_df$Exon <- 1:n
  known_long <- melt(known_df, id.vars = "Exon")
  all_df$Exon <- 1:n
  all_long <- melt(all_df, id.vars = "Exon")
  
  ggplot(data = all_long, aes(x = Exon, y = log2(value), group = variable)) +
    geom_line(color = "grey") +
    geom_line(data = known_long, aes(x = Exon, y = log2(value), group = variable, color = variable)) +
    geom_point(data = known_long, aes(x = Exon, y = log2(value), group = variable, color = variable)) +
    ylab(ylabel) +
    xlab(xlabel) +
    theme(legend.position = legend) +
    ggtitle(title)
}

#and then the cumulative graph - not logged
graph_notlog <- function(gene, project_id, df, title, xlabel, ylabel, legend) {
  known_df <- subset_known(gene, project_id, df)
  all_df <- df
  n <- nrow(df)
  known_df$Exon <- 1:n
  known_long <- melt(known_df, id.vars = "Exon")
  all_df$Exon <- 1:n
  all_long <- melt(all_df, id.vars = "Exon")
  
  ggplot(data = all_long, aes(x = Exon, y = (value), group = variable)) +
    geom_line(color = "grey") +
    geom_line(data = known_long, aes(x = Exon, y = (value), group = variable, color = variable)) +
    geom_point(data = known_long, aes(x = Exon, y = (value), group = variable, color = variable)) +
    ylab(ylabel) +
    xlab(xlabel) +
    theme(legend.position = legend) +
    ggtitle(title)
}

norm_length <- function(gene, df) {
  metadata <- gene_metadata(gene)
  width <- metadata$width
  as.data.frame(df/(width/1000))
}

#calculate the difference between consecutive rows
diff_calc <- function(df) {
  new_df <- as.data.frame(mapply(diff, df))
  a <- nrow(df)
  rownames(new_df) <- rownames(df)[2:a]
  new_df
}

#now I need negative versions of everything
set_sb_neg <- function(gene) {
  sb <- QueryBuilder(compilation = "tcga", regions = gene)
  sb <- set_coordinate_modifier(sb, Coordinates$Within)
  sb <- set_row_filters(sb, strand == "-")
  query_exon(sb)
}

raw_exon_counts_neg <- function(gene) {
  exon <- set_sb_neg(gene)
  df <- as.data.frame(as.matrix(assays(exon)$counts))
  df[dim(df)[1]:1,]
}

#and then we want to subset to project using these functions:
tcga_info_neg <- function(gene) {
  info <- set_sb_neg(gene)
  tcga_info <- cbind(info@colData@rownames, info@colData@listData[["gdc_cases.project.project_id"]], info@colData@listData[["gdc_cases.samples.submitter_id"]], info@colData@listData[["gdc_file_id"]], info@colData@listData[["gdc_cases.samples.sample_type"]])
  tcga_info <- as.data.frame(tcga_info)
  colnames(tcga_info) <- c("rail_id", "project", "tcga_id", "UCSC_id", "sample_type")
  tcga_info$auc <- info@colData@listData[["auc"]]
  tcga_info$tmc <- info@colData@listData[["mapped_read_count"]]
  tcga_info
}


subset_project_neg <- function(gene, project_id, df) {
  info <- tcga_info_neg(gene)
  subset <- filter(info, project == project_id)
  ids <- subset$rail_id
  df[,ids]
}

#to get the metadata
gene_metadata_neg <- function(gene) {
  info <- set_sb_neg(gene)
  df <- as.data.frame(rowRanges(info))
  df[dim(df)[1]:1,]
}

norm_length_neg <- function(gene, df) {
  metadata <- gene_metadata_neg(gene)
  width <- metadata$width
  as.data.frame(df/(width/1000))
}


#function necessary to set up the df
norm_df_neg <- function(gene, df) {
  new_df <- df
  metadata <- gene_metadata_neg(gene)
  new_df$width <- metadata$width
  new_df
}

#calculative cumulative counts, normalised by exon length if that exon is expressed
cumsum_calc_neg <- function(gene, df) {
  a <- ncol(df)
  df <- norm_df_neg(gene, df)
  filt_norm <- function(n) {
    subset_df <- filter(df, df[,n] > 0)
    cumsum_df <- cumsum(subset_df)
    norm_df <- as.data.frame(cumsum_df[,n]/cumsum_df$width*1000)
    rownames(norm_df) <- rownames(cumsum_df)
    colnames(norm_df) <- "rail_id"
    norm_df$rn1 <- as.numeric(rownames(norm_df))
    temp_df <- data.frame(rn1 = as.numeric(rownames(df)))
    new_df <- left_join(temp_df, norm_df)
    final_df <- as.data.frame(as.numeric(new_df[,2]))
    rownames(final_df) <- new_df$rn1
    colnames(final_df) <- colnames(df)[n]
    final_df
  }
  na_df <- as.data.frame(sapply(1:a, filt_norm))
  rownames(na_df) <- rownames(df)
  filled_df <- na_df %>%
    fill(names(na_df))
  filled_df[is.na(filled_df)] <- 0
  filled_df
}

#normalise by auc, subsetted by project
norm_auc_neg <- function(gene, project_id, df) {
  info <- tcga_info_neg(gene)
  info_subset <- filter(info, project == project_id)
  auc <-  info_subset$auc/1000000
  new_df <- as.data.frame(mapply('/', df, auc))
  rownames(new_df) <- rownames(df)
  new_df
}

#doing it all with one function
#title, xlabel, ylabel, legend
condensed_cum_function_neg <- function(gene, project_id, title, xlabel, ylabel, legend) {
  df <- raw_exon_counts_neg(gene)
  proj_df <- subset_project_neg(gene, project_id, df)
  cumsum_df <- cumsum_calc_neg(gene, proj_df)
  cumsum_auc_df <- norm_auc_neg(gene, project_id, cumsum_df)
  #z_calc(cumsum_auc_df)
  #graph_log(gene, project_id, cumsum_auc_df, title, xlabel, ylabel, legend)
}

#creating a difference graph with one function
difference_cum_neg <- function(gene, project_id, title, xlabel, ylabel, legend) {
  df <- raw_exon_counts_neg(gene)
  proj_df <- subset_project_neg(gene, project_id, df)
  cumsum_df <- cumsum_calc_neg(gene, proj_df)
  diff_df <- diff_calc(cumsum_df)
  diff_auc_df <- norm_auc_neg(gene, project_id, diff_df)
  #graph_notlog(gene, project_id, diff_auc_df, title, xlabel, ylabel, legend)
  z_calc(diff_auc_df)
}



