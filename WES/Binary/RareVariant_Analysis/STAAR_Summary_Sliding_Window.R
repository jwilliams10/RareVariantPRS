###########################################################
# Summarization and visualization of sliding window
# analysis results using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################
rm(list=ls())
gc()


Sliding_Window_Results_Summary <- function(agds_dir,jobs_num,input_path,output_path,sliding_window_results_name,
                                           obj_nullmodel,known_loci=NULL,
                                           method_cond=c("optimal","naive"),
                                           QC_label="annotation/filter",geno_missing_imputation=c("mean","minor"),variant_type=c("SNV","Indel","variant"),
                                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                           Use_annotation_weights=FALSE,Annotation_name=NULL,
                                           alpha=0.05,manhattan_plot=FALSE,QQ_plot=FALSE){
  
  ## evaluate choices
  method_cond <- match.arg(method_cond)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  results_sliding_window_genome <- c()
  
  for(chr in 1:22)
  {
    results_sliding_window_genome_chr <- c()
    
    if(chr>1)
    {
      jobs_num_chr <- sum(jobs_num$sliding_window_num[1:(chr-1)])
    }else
    {
      jobs_num_chr <- 0
    }
    
    for(i in 1:jobs_num$sliding_window_num[chr])
    {
      print(i + jobs_num_chr)
      results_sliding_window <- get(load(paste0(input_path,sliding_window_results_name,"_",i+jobs_num_chr,".Rdata")))
      
      results_sliding_window_genome_chr <- rbind(results_sliding_window_genome_chr,results_sliding_window)
    }
    
    results_sliding_window_genome <- rbind(results_sliding_window_genome,results_sliding_window_genome_chr)
  }
  
  rm(results_sliding_window_genome_chr)
  gc()
  
  save(results_sliding_window_genome,file=paste0(output_path,"results_sliding_window_genome.Rdata"))
  
  dim(results_sliding_window_genome)
  
  ### Significant Results
  # alpha <- 0.05
  results_sig <- results_sliding_window_genome[results_sliding_window_genome[,colnames(results_sliding_window_genome)=="STAAR-O"]<alpha,,drop=FALSE]
  
  save(results_sig,file=paste0(output_path,"sliding_window_sig.Rdata"))
  write.csv(results_sig,paste0(output_path,"sliding_window_sig.csv"))
  
  dim(results_sig)[1]
  
  if(length(known_loci)!=0)
  {
    results_sig_cond <- c()
    if(length(results_sig)!=0)
    {
      for(kk in 1:dim(results_sig)[1])
      {
        chr <- as.numeric(results_sig[kk,1])
        start_loc <- as.numeric(results_sig[kk,2])
        end_loc <- as.numeric(results_sig[kk,3])
        
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)
        
        res_cond <- Sliding_Window_cond(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                        start_loc=start_loc,end_loc=end_loc,known_loci=known_loci,method_cond=method_cond,
                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                        Annotation_name_catalog=Annotation_name_catalog,Annotation_dir=Annotation_dir,
                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        results_sig_cond <- rbind(results_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    save(results_sig_cond,file=paste0(output_path,"sliding_window_sig_cond.Rdata"))
    write.csv(results_sig_cond,paste0(output_path,"sliding_window_sig_cond.csv"))
  }
  
  ## manhattan plot
  if(manhattan_plot)
  {
    print("Manhattan plot")
    
    png(paste0(output_path,"sliding_window_manhattan.png"), width=12, height=8, units = 'in', res = 600)
    
    print(manhattan_plot(as.numeric(results_sliding_window_genome[,1]), (as.numeric(results_sliding_window_genome[,2])+as.numeric(results_sliding_window_genome[,3]))/2, as.numeric(results_sliding_window_genome[,colnames(results_sliding_window_genome)=="STAAR-O"]), col = c("blue4", "orange3"),sig.level=0.05/dim(results_sliding_window_genome)[1]))
    
    dev.off()
  }
  
  ## Q-Q plot
  if(QQ_plot)
  {
    print("Q-Q plot")
    
    ## STAAR-O
    observed <- sort(as.numeric(results_sliding_window_genome[,colnames(results_sliding_window_genome)=="STAAR-O"]))
    lobs <- -(log10(observed))
    
    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))
    
    png(paste0(output_path,"sliding_window_qq_staar_o.png"), width=8, height=8, units = 'in', res = 600)
    
    par(mar=c(5,6,4,4))
    
    plot(lexp,lobs,pch=20, cex=1, xlim = c(0, max(lexp)), ylim = c(0, max(lobs)),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)
    
    abline(0, 1, col="red",lwd=2)
    
    dev.off()
    
    ## S(1,25)
    observed <- sort(as.numeric(results_sliding_window_genome[,colnames(results_sliding_window_genome)=="SKAT(1,25)"]))
    lobs <- -(log10(observed))
    
    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))
    
    png(paste0(output_path,"sliding_window_qq_skat_1_25.png"), width=8, height=8, units = 'in', res = 600)
    
    par(mar=c(5,6,4,4))
    
    plot(lexp,lobs,pch=20, cex=1, xlim = c(0, max(lexp)), ylim = c(0, max(lobs)),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)
    
    abline(0, 1, col="red",lwd=2)
    
    dev.off()
    
    
    ## B(1,1)
    observed <- sort(as.numeric(results_sliding_window_genome[,colnames(results_sliding_window_genome)=="Burden(1,1)"]))
    lobs <- -(log10(observed))
    
    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))
    
    png(paste0(output_path,"sliding_window_qq_burden_1_1.png"), width=8, height=8, units = 'in', res = 600)
    
    par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=20, cex=1, xlim = c(0, max(lexp)), ylim = c(0, max(lobs)),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)
    abline(0, 1, col="red",lwd=2)
    
    dev.off()
    
  }
}







## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  ###########################################################
  #           User Input
  ###########################################################
  ## Number of jobs for each chromosome
  jobs_num <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/jobs_num.Rdata"))
  ## aGDS directory
  agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))
  ## Null model
  obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/nullmodels_staar/",trait,"_Train_Null_Model.RData")))
  ### Known loci
  known_loci <- NULL
  
  ## results path
  input_path <- "/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/"
  output_path <- input_path
  ## results name
  sliding_window_results_name <- paste0(trait,"_UKBB_WES_Sliding_Train")
  
  ## QC_label
  QC_label <- "annotation/info/QC_label"
  ## geno_missing_imputation
  geno_missing_imputation <- "mean"
  ## method_cond
  method_cond <- "optimal"
  ## variant_type
  variant_type <- "SNV"
  ## alpha level
  alpha <- 1
  
  ## Annotation_dir
  Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
  ## Annotation channel
  Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/Annotation_name_catalog.Rdata"))
  ## Use_annotation_weights
  Use_annotation_weights <- TRUE
  ## Annotation name
  Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                       "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
  
  ###########################################################
  #           Main Function 
  ###########################################################
  Sliding_Window_Results_Summary(agds_dir=agds_dir,jobs_num=jobs_num,
                                 input_path=input_path,output_path=output_path,sliding_window_results_name=sliding_window_results_name,
                                 obj_nullmodel=obj_nullmodel,known_loci=known_loci,
                                 method_cond=method_cond,
                                 QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                 alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)

  file.rename("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/sliding_window_sig.csv",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_sliding_window_sig.csv"))
  
}