#Code for DTU analysis using a combination of isoformSwitchAnalyzeR, DEXseq and stageR

#Procedure
#1) take in sample information with file location and read in files to create switchlist
#2) prefilter this list
#3) read in object for DEXseq and subset and run analysis
#4) run stageR object
#5) generate boxplot
#6) generate switchplot after taking in consequence into consideration

#check necessary libraries have been loaded
require(dplyr)
require(tibble)
require(tidyr)
require(IsoformSwitchAnalyzeR)
require(DEXSeq)
require(stageR)
require(tximport)
require(DRIMSeq)
require(stringr)
require(ggplot2)

#Functions
#Take in input and filter out low expression isoforms based on isoformSwitchAnalyzeR

prefilter <- function(sample_vector = NULL, metadata_info = NULL, sample_id_column_name = NULL, output_suffix = NULL){
  print("Starting analysis...")
  
  setwd("~")
  salmon_quant <- IsoformSwitchAnalyzeR::importIsoformExpression(sampleVector = sample_vector)
  
  mydesign <- data.frame(sampleID = colnames(salmon_quant$abundance)[-1])
  mydesign <- merge(mydesign, metadata_info, by.x = "sampleID", by.y = as.character(sample_id_column_name)) %>% 
    dplyr::select(sampleID, condition:colnames(metadata_info)[ncol(metadata_info)])
  print("Summary of the design that will be used is:")
  str(mydesign)
  
  print("Creating switch object...")
  switchlist <- IsoformSwitchAnalyzeR::importRdata(isoformCountMatrix = salmon_quant$counts,
                                                   isoformRepExpression = salmon_quant$abundance,
                                                   designMatrix = mydesign,
                                                   isoformExonAnnoation = "~/../Desktop/salmon/Homo_sapiens.GRCh38.84.chr_patch_hapl_scaff.gtf",
                                                   isoformNtFasta = "~/../Desktop/salmon/Homo_sapiens.GRCh38.all_tsc.fa",
                                                   removeNonConvensionalChr = T,
                                                   ignoreAfterPeriod = T,
                                                   quiet = F)
  print("Summary of the created switchlist:")
  summary(switchlist)
  
  print("Summary of the filtered switchlist:")
  switchlist <- IsoformSwitchAnalyzeR::preFilter(switchlist,
                          geneExpressionCutoff = 3,
                          isoformExpressionCutoff = 0,
                          removeSingleIsoformGenes = TRUE,
                          quiet = F)
  summary(switchlist)
  
  print(paste("Your switchlist has been created with the name: switchlist", output_suffix, sep = "_"))
  assign(switchlist, x = paste("switchlist", output_suffix, sep = "_"), envir = .GlobalEnv)
  print("~~~DONE~~~")
}




# Run stageR

run_stager <- function(transcript_pval = NULL, gene_qval = NULL, dataset_name = NULL, output_suffix = NULL, cutoff = 0.05, allow_NAs  = FALSE){
  pConfirmation = matrix(transcript_pval$pvalue, ncol=1)
  dimnames(pConfirmation) = list(transcript_pval$featureID,"transcript")
  
  pScreen = gene_qval$qval %>% as.vector()
  names(pScreen) = gene_qval$gene
  
  if(any(is.na(gene_qval$qval) == F)){
    stageRObj = stageRTx(pScreen = pScreen, 
                         pConfirmation = pConfirmation, 
                         pScreenAdjusted = TRUE, 
                         tx2gene = transcript_pval[1:2])

    stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = as.numeric(cutoff), allowNA = allow_NAs)
    dex.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = F)
    dex.padj = merge(transcript_pval, dex.padj, by.x = c("groupID","featureID"), by.y = c("geneID","txID"))
  } else {
    print("qvals were found to be zero, please rectify before running this function")
  }
  assign(dex.padj, x = paste(dataset_name, "dtu_result", output_suffix, sep = "_"), envir = .GlobalEnv)
  
  print(paste("p adjusted values were generated and can be found within: ", 
              "dtu_result_", 
              output_suffix, 
              sep = ""))
        
  print("Result summary:")
  print(paste("Significant genes after stageR screening stage: ", length(unique(dex.padj[dex.padj$gene < as.numeric(cutoff),]$groupID))))
  print("Significant transcripts after stageR confirmation stage: ")
  return(table(dex.padj$transcript < as.numeric(cutoff)))
  
  print("~~~DONE~~~")
}

mod_for_stager <- function(dtu_result = NULL){
  tsc_p_val <- dtu_result %>% dplyr::select(featureID, groupID, pvalue)
  gene_q_val <- dtu_result %>% dplyr::select(groupID, gene) %>% setNames(c("gene", "qval")) %>% unique()
  assign(tsc_p_val, x = "tsc_p_val", envir = .GlobalEnv)
  assign(gene_q_val, x = "gene_q_val", envir = .GlobalEnv)
}

# Visualization: boxplot for significant switches

plot_expression <- function(countData = NULL, dtu_result = NULL, geneID = NULL, sample_info = NULL, logTransform = FALSE) {
  
  colnames(countData)[1:2] <- c("gid","tid")
  sub <- subset(countData, gid == geneID)
  sub <- reshape2::melt(sub, id = c("gid", "tid"))
  sub <- merge(sample_info, sub, by.x = "sample_id", by.y = "variable")
  sub <- merge(sub, dtu_result, by.x = "tid", by.y = "featureID")
  y.pos <- sub %>% group_by(tid) %>% summarise(max_exp = max(value)) %>% mutate(y.position = max_exp + (max_exp/4))
  sub <- sub %>% mutate(p.adj.signif = case_when(transcript > 0.05 ~ "NS", 
                                                 transcript <0.05 & transcript >= 0.01 ~ "*", 
                                                 transcript < 0.01 & transcript >= 0.001 ~ "**", 
                                                 transcript < 0.001 & transcript >= 0.0001 ~ "***", 
                                                 transcript < 0.0001 ~ "****", TRUE ~ "NS")) %>% 
    merge(., y.pos, by = "tid") 
  
  if(logTransform) {
    sub$value = log10(sub$value + 0.01)
    sub$y.position = log10(sub$y.position + 0.01)
  }
  
  clrs = c("dodgerblue3", "cadetblue")
  
  p <- ggplot(sub, aes(tid, value, color = condition, fill = condition)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.8, lwd = 0.5) + 
    geom_jitter(position = position_jitterdodge()) +
    stat_summary(fun = mean, geom = "point", color = "black", shape = 2, size = 3, position=position_dodge(width = 0.8)) + 
    scale_color_manual(values = clrs) + 
    scale_fill_manual(values = clrs) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(geneID) + 
    xlab("Transcripts") + 
    geom_text(data = (sub %>% filter(p.adj.signif != "NS")), aes(x = tid, y = y.position, label = p.adj.signif), 
              size = 4, colour = "red4")
  
  if(logTransform) {
    p = p + ylab("log10(Expression + 0.01)")
  } else {
    p = p + ylab("Normalised Expression")
  }
  p
}


plot_dtu_expression <- function(dxd_object = NULL, dtu_result = NULL, sample_info = NULL, cutoff = 0.05, n = NULL, gene_id = NULL, outname = "dtu_result", logtransform = FALSE){
  samp_num <- dplyr::n_distinct(dxd_object@modelFrameBM$sample_id) %>% as.numeric()
  dex.norm <- cbind(as.data.frame(stringr::str_split_fixed(rownames(counts(dxd_object)), ":", 2)), 
                    as.data.frame(counts(dxd_object, normalized = TRUE))[,1:samp_num])
  colnames(dex.norm) <- c("groupID", "featureID", as.character(colData(dxd_object)$sample_id)[1:samp_num])
  row.names(dex.norm) <- NULL
  
  if(is.null(gene_id) == TRUE){
    
    dtu_result_sub <- dtu_result %>% filter(gene < as.numeric(cutoff)) %>% arrange(transcript, gene)
    gene_ids <- unique(dtu_result_sub$groupID) %>% as.vector()
    
    if(is.null(n) == TRUE){
      plots <- c()
      for(i in 1:n_distinct(dtu_result_sub$groupID)){
        plots[[i]] = plot_expression(countData = dex.norm, dtu_result = dtu_result, geneID = gene_ids[i], sample_info = sample_info, logTransform = logtransform)
      }
    } else {
      if(n <= dplyr::n_distinct(dtu_result_sub$groupID)){
        plots <- c()
        for(i in 1:n){
          plots[[i]] = plot_expression(countData = dex.norm, dtu_result = dtu_result, geneID = gene_ids[i], sample_info = sample_info, logTransform = logtransform)
        }
      } else {
      return(print(paste("Requested n is too large, choose a number less than or equal to", n_distinct(dtu_result_sub$groupID), sep = " ")))
      }
    }
  } else {
    gene_id <- as.vector(gene_id) 
    plots <- c()
    for(i in seq(1, length(gene_id))){
      plots[[i]] = plot_expression(countData = dex.norm, dtu_result = dtu_result, geneID = gene_id[i], sample_info = sample_info, logTransform = logtransform)
    }
  }
  
  pdf(paste(outname, "pdf", sep = "."), width = 8, height = 5, onefile = T)
  for(i in plots){
    print(i)
  }
  dev.off()
  
  return(print(paste("Plots were generated within ", outname, ".pdf", sep = "")))
  
  print("~~~DONE~~~")
}

# Visualization: switchplot for genes of interest

switchlist_sub <- function(switchlist = NULL, dtu_result = NULL, outname ="sub"){
  originalName <- deparse(substitute(switchlist))
  dtu_to_keep <- dtu_result %>% filter(gene < 0.05)
  gene_ids_ver <- switchlist[["isoformFeatures"]][["gene_id"]] %>% as.data.frame() %>% setNames("gene_id_ver")
  gene_ids_ver <- gene_ids_ver %>% separate(gene_id_ver, into = c("gene_id", NA), sep = "\\.", remove = F)
  gene_ids_ver <- gene_ids_ver %>% filter(gene_id %in% dtu_to_keep$groupID) %>% unique()
  switchlist_sub <- subsetSwitchAnalyzeRlist(switchlist, 
                                             subset = switchlist$isoformFeatures$gene_id %in% gene_ids_ver$gene_id_ver)
  assign(switchlist_sub, x = paste(originalName, outname, sep = "_"), envir = .GlobalEnv)
}

#Part I: generating fasta files for external run 

switchplot_part1 <- function(switchlist = NULL, dtu_result = NULL, fasta_suffix = "extracted_sequence"){
  
  originalName <- deparse(substitute(switchlist))
  switchlist_gene_ids <- switchlist$isoformFeatures$gene_id %>% 
    as.data.frame() %>% 
    setNames("gene_id_ver") %>%  
    tidyr::separate(., col = "gene_id_ver", into = c("gene_id", NA), sep = "\\.", remove = F) %>% 
    unique()
  
  gene_ids_to_keep <- dtu_result %>% filter(gene < 0.05) %>% distinct(groupID)
  switchlist_gene_ids_to_keep <- switchlist_gene_ids %>% 
    filter(gene_id %in% gene_ids_to_keep$groupID) %>% 
    dplyr::select(gene_id_ver)
  
  switchlist_tsc_ids <- switchlist[["isoformFeatures"]][["isoform_id"]] %>% 
    as.data.frame() %>% 
    setNames("tsc_id")
  
  if(identical(dtu_result$featureID, switchlist_tsc_ids$tsc_id) == TRUE){

    switchlist$isoformFeatures$gene_switch_q_value <- dtu_result$gene
    switchlist$isoformFeatures$isoform_switch_q_value <- dtu_result$transcript
    switchlist <- subsetSwitchAnalyzeRlist(switchlist, 
                                           subset = switchlist$isoformFeatures$gene_id %in% switchlist_gene_ids_to_keep$gene_id_ver)
    
    
    switchlist <- extractSequence(switchlist, 
                                  onlySwitchingGenes = F, 
                                  removeLongAAseq = T, 
                                  alsoSplitFastaFile = T, 
                                  outputPrefix = as.character(fasta_suffix))
    
  } else {
    return(print("Switchlist and dtu result entries are not in the same order, Please check!"))
  }
  
  print('Fasta files were saved in the current directory')
  assign(switchlist, x = originalName, envir = .GlobalEnv)
  print("~~~DONE~~~")
}

#Part II: generating switch plots:
switchplot_part2 <- function(switchlist = NULL, n = NULL, gene_ids = NULL, gene_names = NULL, outname = "switchplots",
                             analyze_cpc2 = TRUE, analyze_pfam = TRUE, analyze_iupred = TRUE, analyze_signalp = TRUE, 
                             pathtocpc2file = NULL, pathtoiupred = NULL, pathtopfam = NULL, pathtosignalp = NULL, 
                             consequences = c('intron_retention',
                             'coding_potential',
                             'ORF_seq_similarity',
                             'NMD_status',
                             'domains_identified',
                             'IDR_identified',
                             'IDR_type',
                             'signal_peptide_identified'), analyzeSplicing = TRUE){
  
  source("switchplottranscript_2_mod.R")
  source("analyzepfam.R")
  originalName <- deparse(substitute(switchlist))
  if(analyze_cpc2 == TRUE){
    switchlist <- IsoformSwitchAnalyzeR::analyzeCPC2(switchlist, pathToCPC2resultFile = pathtocpc2file, 
                                                     removeNoncodinORFs = F, quiet = T)
  } 

  if(analyze_pfam == TRUE){
    switchlist <- analyzePFAM(switchlist, pathToPFAMresultFile = pathtopfam, 
                                                     showProgress = F, quiet = F)
  }
  
  if(analyze_signalp == TRUE){
    switchlist <- IsoformSwitchAnalyzeR::analyzeSignalP(switchlist, pathToSignalPresultFile = pathtosignalp, 
                                                        quiet = T)
  }
  
  if(analyze_iupred == TRUE){
    switchlist <- IsoformSwitchAnalyzeR::analyzeIUPred2A(switchlist, pathToIUPred2AresultFile = pathtoiupred,  
                                                     showProgress = F)
  }
  
  if(analyzeSplicing == TRUE){
    switchlist <- analyzeAlternativeSplicing(switchlist, onlySwitchingGenes = F, quiet = T, showProgress = F)
  } else {
    switchlist <- switchlist
  }
  
  
  n_max <- switchlist$isoformFeatures$gene_id %>% as.data.frame() %>%
    setNames("gene_id_vers") %>% dplyr::n_distinct()
  
  gene_id_vers <- switchlist$isoformFeatures$gene_id %>% as.data.frame() %>%
    setNames("gene_id_vers") %>% dplyr::distinct()
  
  if(is.null(gene_ids) == FALSE){
    gene_ids_for_plot <- as.character(gene_ids)
    plots <- c()
    for (i in 1:length(gene_ids_for_plot)){
      plots[[i]] <- switchPlotTranscript_2_mod(switchlist, gene = gene_ids_for_plot[i])
    }
  } else if (is.null(gene_names) == FALSE){
    gene_names_for_plot <- as.character(gene_names)
    plots <- c()
    for (i in 1:length(gene_names_for_plot)){
      plots[[i]] <- switchPlotTranscript_2_mod(switchlist, gene = gene_names_for_plot[i])
    }
  } else if (is.null(n) == FALSE){
    if(n > n_max){
      return(print(paste("Requested n is too large, choose a number less than or equal to", 
                           n_max, sep = " ")))
      } else {
        plots <- c()
        for (i in 1:n){
          plots[[i]] <- switchPlotTranscript_2_mod(switchlist, gene = gene_id_vers$gene_id_vers[i])
        }
      }
  } else {
    plots <- c()
    for (i in 1:nrow(gene_id_vers)){
      plots[[i]] <- switchPlotTranscript_2_mod(switchlist, gene = gene_id_vers$gene_id_vers[i])
    }
  }
  
  pdf(paste(outname, "pdf", sep = "."), width = 8, height = 5, onefile = T)
  for(i in plots){
    print(i)
  }
  dev.off()
  
  return(print(paste("Plots were generated within ", outname, ".pdf", sep = "")))
  assign(switchlist, x = originalName, envir = .GlobalEnv)
  
  print("~~~DONE~~~")
  
}

# Isoform fraction plot of expression significance


IF_plot <- function(IF_data = NULL, geneID = NULL, IF_expr = NULL, samp_info = NULL, clear = FALSE){
  sub <- subset(IF_data, groupID == geneID)
  sub <- sub %>% mutate(p.adj.signif = case_when(transcript > 0.05 ~ "NS", 
                                                 transcript <0.05 & transcript >= 0.01 ~ "*", 
                                                 transcript < 0.01 & transcript >= 0.001 ~ "**", 
                                                 transcript < 0.001 & transcript >= 0.0001 ~ "***", 
                                                 transcript < 0.0001 ~ "****", TRUE ~ "NS"))
  sub <- sub %>% dplyr::select(featureID, groupID, gene_name, p.adj.signif, IF1_val, IF2_val, gene)
  sub_exp <- sub %>% merge(., IF_expr, by.x = "featureID", by.y = "isoform_id", all.x = T) %>% 
    dplyr::select(-c(IF1_val, IF2_val, gene, p.adj.signif, gene_name, groupID))
  samp_cond <- samp_info %>% dplyr::select(sample_id_match, condition)
  
  sub <- sub %>% reshape2::melt(., id.vars = c("featureID", "groupID", "gene_name", "p.adj.signif", "gene")) %>% 
    mutate(condition = if_else(variable == "IF1_val", "CNT", "MDD")) 
  sub_exp <- sub_exp %>% reshape2::melt(., id.vars = c("featureID")) %>%
    merge(., samp_cond, by.x = "variable", by.y = "sample_id_match", all.x = T) 
  sub_exp <- sub_exp %>% filter(value != "NaN")
  
  y.pos <- sub_exp %>% group_by(featureID) %>% summarise(max_exp = max(value)) %>% mutate(y.position = max_exp + (max_exp/4))
  sub <- merge(sub, y.pos, by = "featureID")
  
  if(clear == FALSE){
    p <- ggplot(sub, aes(featureID, value, fill = as.factor(condition))) +
      geom_bar(stat = "identity", position = "dodge", width = 0.8)  + 
      geom_point(data = sub_exp, aes(x = featureID, y = value, color = as.factor(condition)), 
                 position = position_jitterdodge(jitter.width = 0.3, dodge.width=0.9), 
                 shape = 20, size = 0.3) + 
      scale_fill_manual(values = alpha(c("dodgerblue3", "cadetblue"), 0.5),name = "") +
      scale_color_manual(values = c("royalblue4", "palegreen4"), name = "") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90)) +
      labs(x = "", y = "Isoform Fraction (IF)", 
           title = paste(unique(sub$gene_name), " (", unique(sub$groupID), ") : ", "q-val = ", round(unique(sub$gene), 4), sep = "")) + 
      geom_text(data = (sub %>% filter(p.adj.signif != "NS")), 
                aes(x = featureID, y = y.position, label = p.adj.signif), 
                size = 4, colour = "red4")
  } else {
    p <- ggplot(sub, aes(featureID, value, fill = as.factor(condition))) +
      geom_bar(stat = "identity", position = "dodge", width = 0.8)  + 
      scale_fill_manual(values = c("dodgerblue3", "cadetblue"),name = "") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90)) +
      labs(x = "", y = "Isoform Fraction (IF)", 
           title = paste(unique(sub$gene_name), " (", unique(sub$groupID), ") : ", "q-val = ", round(unique(sub$gene), 4), sep = "")) + 
      geom_text(data = (sub %>% filter(p.adj.signif != "NS")), 
                aes(x = featureID, y = y.position, label = p.adj.signif), 
                size = 4, colour = "red4")
  }
  
  
  return(p)
  
}   

plot_IF_dtu <- function(switchlist = NULL, dtu_result = NULL, samp_info = NULL, cutoff = 0.05, 
                        n = NULL, gene_id = NULL, clear = FALSE, outname = "IF_dtu"){
  samp_num <- ncol(switchlist[["isoformCountMatrix"]] %>% as.data.frame())-1
  IF_data_frame <- data.frame(groupID = switchlist[["isoformFeatures"]][["gene_id"]],
                              gene_name = switchlist[["isoformFeatures"]][["gene_name"]],
                              featureID = switchlist[["isoformFeatures"]][["isoform_id"]],
                              IF1_val = switchlist[["isoformFeatures"]][["IF1"]], 
                              IF2_val = switchlist[["isoformFeatures"]][["IF2"]])
  IF_details <- IF_data_frame %>% merge(., dtu_result, by = "featureID")
  IF_expression <- switchlist[["isoformRepIF"]]
  
  if(is.null(gene_id) == TRUE){
    
    IF_sub <- IF_details %>% filter(gene < as.numeric(cutoff)) %>% arrange(transcript, gene) %>% separate(groupID.x, into = c("groupID", NA), sep = "\\.")
    gene_ids <- unique(IF_sub$groupID) %>% as.vector()
    
    if(is.null(n) == TRUE){
      plots <- c()
      for(i in 1:n_distinct(IF_sub$groupID)){
        plots[[i]] = IF_plot(IF_data = IF_sub, geneID = gene_ids[i], 
                             IF_expr = IF_expression, samp_info = samp_info, clear = clear)
      }
    } else {
      if(n <= dplyr::n_distinct(IF_sub$groupID)){
        plots <- c()
        for(i in 1:n){
          plots[[i]] = IF_plot(IF_data = IF_sub, geneID = gene_ids[i], 
                               IF_expr = IF_expression, samp_info = samp_info, clear = clear)
        }
      } else {
        return(print(paste("Requested n is too large, choose a number less than or equal to", n_distinct(IF_sub$groupID), sep = " ")))
      }
    }
  } else {
    gene_id <- as.vector(gene_id) 
    plots <- c()
    for(i in seq(1, length(gene_id))){
      plots[[i]] = IF_plot(IF_data = IF_sub, geneID = gene_id[i], 
                           IF_expr = IF_expression, samp_info = samp_info, clear = clear)
    }
  }
  
  pdf(paste(outname, "pdf", sep = "."), width = 8, height = 5, onefile = T)
  for(i in plots){
    print(i)
  }
  dev.off()
  
  return(print(paste("Plots were generated within ", outname, ".pdf", sep = "")))
  
  print("~~~DONE~~~")
}






