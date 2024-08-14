rm(list = ls())

system("dx download -r UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/")
system("dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt")

pheno_train <- read.delim("All_Train.txt")

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =18),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

library(plotrix)
library(data.table)
library(RColorBrewer)
library(optparse)
library(dplyr)
library(ggplot2)
library(cowplot)

lambda_dat <- NULL

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  dat <- read.delim(paste0("GWAS_SummaryStats/",trait,"_sumstats.",trait,".glm.linear"), header=FALSE, comment.char="#")
  colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
  dat <- dat[dat$TEST == "ADD",]
  
  dat <- dat[,c("CHROM","ID","REF","POS","A1","BETA","P","A1_FREQ")]
  colnames(dat) <- c("CHR","SNP","REF","BP","A1","BETA","P","A1_FREQ") 
  dat$MAF <- ifelse(dat$A1_FREQ <= 0.5, dat$A1_FREQ,1-dat$A1_FREQ)
  dat <- dat[dat$MAF > 0.01,]
  
  
  x <- dat$P
  z <- qnorm(x / 2)
  lambda <- round(median(z^2) / qchisq(0.5,1), 3)
  lambda_1000 <- round(1+1000*(lambda-1)/sum(!is.na(pheno_train[,trait])) ,3)
  
  lambda_dat <- rbind(lambda_dat,data.frame(Trait = trait,Ancestry_Group = "European",lambda = lambda, lambda_1000 = lambda_1000, datasource = "UKB WES"))
  
  p.pwas <- 5E-08
  
  nCHR <- length(unique(dat$CHR))
  dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(dat$CHR)){
    nbp[i] <- max(dat[dat$CHR == i,]$BP)
    dat$BPcum[dat$CHR == i] <- dat$BP[dat$CHR == i] + s
    s <- s + nbp[i]
  }
  axis.set <- dat %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ylim <- abs(floor(log10(min(dat$P)))) + 2 
  sig1 <- p.pwas
  
  sigline <- data.frame(sig=c(-log10(sig1)),val=c(paste0("P=",signif(sig1,2))))
  
  # layout(matrix(c(1,2),ncol = 2),widths = c(2,1))
  
  pdf(paste0(trait,"_UKB_WES_CV_Manhattan_Plot.pdf"), width=15, height=9)
  
  p1 <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                        color = as.factor(CHR), size = -log10(P))) +
    geom_point(alpha = 0.8, size=0.8) + 
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    geom_hline(data = sigline, aes(yintercept = sig), color= "red", linetype="dashed") +
    guides(color = F) + 
    labs(x = NULL, 
         y = "-log10(p)", 
         linetype = "",
         title = paste0(trait," for Europeans"))+
    theme_Publication()+
    theme(
      legend.position = "top",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5),
      axis.text.y = element_text(size = 16),
      plot.subtitle = element_text(size = 8)
    )
  
  print(p1)
  
  dev.off()
  
  pdf(paste0(trait,"_UKB_WES_CV_QQplot.pdf"), width=10, height=10)
  
  qqplotdata <- function(logpvector){
    o = sort(logpvector,decreasing=T)
    e = -log10(ppoints(length(o)))       
    qqdata <- data.frame(o,e)
    qqdata$o <- round(qqdata$o,3)
    qqdata$e <- round(qqdata$e,3)
    keepU <- which(!duplicated(qqdata))
    qqdata <- qqdata[keepU,]
    
    N <- length(logpvector) ## number of p-values
    ## create the confidence intervals
    qqdata$c975 <- NA
    qqdata$c025 <- NA
    
    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)
    
    for(i in 1:length(keepU)){
      j <- keepU[i]
      qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
      qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
    }
    return(qqdata)
  }
  yLine <- c(-log10(5E-8))
  colLine <- c("red")
  dat$log10P <- -log10(dat$P)
  gwas <- as.data.frame(dat)
  # Determine frequency bins and create variable for binned QQ plot
  
  minMAF <- min(gwas$MAF)
  
  freqbins <- c(c(0.5,0.05,0.01,0)[which(c(0.5,0.05,0.01,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
  gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
  freqtable <- table(gwas$freqbin)
  freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
  freqtable <- freqtable[freqtable > 0]
  
  ## Generate QQ plot data by frequency bin
  fbin <- character(0)
  fN <- integer(0)
  fx <- numeric(0)
  fy <- numeric(0)
  fcol <- character(0)
  legendcol <- character(0)
  conf <- list()
  allcols <- brewer.pal(4,"Set1")
  ycol <- "log10P"
  for(f in 1:length(freqtable)){
    fbin <- c(fbin,names(freqtable)[f])
    fsnps <- which(gwas$freqbin ==names(freqtable)[f])
    plotdata <- qqplotdata(gwas[[ycol]][fsnps])
    fN <- c(fN,freqtable[f])
    fx <- c(fx,plotdata$e)
    fy <- c(fy,plotdata$o)
    fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
    conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                            'y'=c(plotdata$c975,rev(plotdata$c025)))
    legendcol <- c(legendcol,allcols[f])
  }
  legendtext <- paste0("MAF=",fbin,"; N SNPs=",format(fN,big.mark=",",scientific=FALSE))
  opt <-  list(break.top = 15,
               top.size = 0.125)
  
  
  xlim <- c(0,max(fx,na.rm=T))
  ylim <- c(0,max(fy,na.rm=T))
  maxY <- max(fy,na.rm=T)
  print("okkkk2")
  par(mar=c(5.1,5.1,4.1,1.1))
  print("okkkk3")
  
  lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
  lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
  lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
  lab2 <- lab2[lab2 > max(lab1)]
  
  # resulting range of top scale in bottom scale units
  top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
  top.data = max(lab2)-opt$break.top
  
  # function to rescale the top part
  rescale <- function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
  rescaled.y <- rescale(fy[fy>opt$break.top])
  plot(0,0,
       ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
       xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)
  
  # Plot confidence intervals	
  for(p in 1:length(conf)){
    polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
            col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
            border = NA)
  }
  
  # add points
  points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)
  
  # identify line & add axis break
  lines(xlim,xlim,col="black",lty = 2)
  axis(1,cex.axis=1.5,cex.lab=1.5)
  par(las=1)
  axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
  par(las=0)
  box()
  par(las=0)
  points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
  par(las=1)
  axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
  axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
  axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
  lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
  abline(h=ifelse(yLine<opt$break.top,
                  yLine,
                  rescale(yLine)),
         col=colLine,lwd=1.5,lty=2)
  legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
  text(4,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
  text(4.5,1,paste(lambda_1000),cex = 1.5)
  title(paste0(trait," for European"))
  
  # p2 <- recordPlot()
  # 
  # title <- ggdraw() + 
  #   draw_label(
  #     paste0(trait," for Europeans"),
  #     fontface = 'bold',
  #     x = 0,
  #     hjust = 0,
  #     size = 30
  #   ) +
  #   theme(
  #     # add margin on the left of the drawing canvas,
  #     # so title is aligned with left edge of first plot
  #     plot.margin = margin(0, 0, 0, 7)
  #   )
  
  # plot_row <- plot_grid(p1, p2, align = "h", axis = "bt", rel_widths = c(6, 1))
  # 
  # p3 <- plot_grid(
  #   title, plot_row,
  #   ncol = 1,
  #   # rel_heights values control vertical title margins
  #   rel_heights = c(0.1, 1)
  # )
  # 
  # cowplot::save_plot(filename = paste0(trait,"_UKB_WES_CV_QQplot.pdf"),plot = p3,nrow = 2,ncol = 2,base_width = c(24))
  # 
  dev.off()
}

write.csv(lambda_dat,file = "UKB_WGS_lambda.csv",row.names = FALSE)