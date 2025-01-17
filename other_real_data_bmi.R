
sum_stats <- read.table("/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/BMI/bmi.txt", header=T)
head(sum_stats)
dim(sum_stats)


cslength2 <- numeric()
cslength <- numeric()



chrtotal <- as.character(1:22)



for(chrid in chrtotal){
  
  
  print(chrid)
  reference <- read.table("/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/UKBB_panel/hg38/reference.txt", header=T)
  referchr <- reference[reference$Chr==chrid,]
  
  directory_path <- paste("/gpfs/gibbs/pi/zhao/zy92/data/gtex/adjusted_expr/chr",chrid, sep="")
  # List all folders in the specified directory
  folder_names <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)
  genes <- sub("\\..*", "", folder_names)
  
  # Print the folder names
  #print(folder_names)
  
  library(biomaRt)
  
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  z <- getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
             "ensembl_gene_id", genes, mart = ensembl)
  index <- match(z$ensembl_gene_id, genes)
  
  #check
  
  table(genes[index] == z$ensembl_gene_id)
  
  z$folder_names <- folder_names[index]
  
  
  mt <- 2
  
  d <- 0.001
  d1 <- 0.01
  
  path <- paste('/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/GTEX/bychr/gtex_ukb_chr', chrid, sep="")
  library(genio)
  dat <- read_plink(path)
  
  prior <- numeric(nrow(dat$X)) + d
  eqtl <- cbind(prior, dat$bim$pos)
  
  colnames(eqtl) <- c("prior_wb","bp_hg38")
  
  count <- numeric(0)
  checkna <- numeric()
  leng_snp <- numeric()
  #[c(-164,-195,-434,-449,-472,-126,-147)]
  #c(-2,-43,-859,-866)
  for(j in (1:nrow(z))){
    gene_info <- z[j,]
    path <- paste("/gpfs/gibbs/pi/zhao/zy92/data/gtex/adjusted_expr/chr",chrid,"/",gene_info$folder_names,"/Adipose_Subcutaneous.adj_expr", sep="")
    if (!file.exists(path)) {
      #cat("File does not exist for iteration", j, "\n")
      next  # Skip this iteration and move to the next one
    }
    ge_j <- read.table(path)
    #length(intersect(ge_j$V1, dat$fam$fam))
    index <- match(ge_j$V1, dat$fam$fam) # 
    #check:
    #table(dat$fam$fam[index] == ge_j$V1)
    checkna <- append(checkna, length(which(is.na(index))))
    
    indexsnp <- which(dat$bim$pos>gene_info$start_position & dat$bim$pos<gene_info$end_position)
    
    leng_snp <- append(leng_snp, length(indexsnp))
    
    if(length(indexsnp==1)){
      eqtl[indexsnp,1]=1
    }

    if (length(indexsnp)<mt) {
      #cat("Gene contains less than 10 SNPs for iteration", j, "\n")
      next  # Skip this iteration and move to the next one
    }
    
    X_j <- dat$X[indexsnp,index]
    X_j <- bedNA(t(X_j))
    X_j <- scale(X_j)
    
    # check
    #table(rownames(X_j)==ge_j$V1)
    #table(colnames(X_j)==dat$bim[indexsnp,2])
    
    fit <- susie(X_j,scale(ge_j$V2),L=round(d1*length(indexsnp)+1),standardize = FALSE,
                 estimate_residual_variance = FALSE,
                 scaled_prior_variance = 0.5,
                 min_abs_corr = 0.5)
    priore1 <- numeric(ncol(fit$alpha))
    for(jj in 1:ncol(fit$alpha)){
      priore1[jj] <- 1-prod(1-fit$alpha[,jj])
    }
    priore1[priore1 == 0] <- 10^(-10)
    eqtl[indexsnp,1] <- priore1
    #  print(paste("na:",length(which(is.na(index)))))
    count <- append(count, j)
  }
  
  print(table(checkna))
  print(paste(chrid,':',length(count)))
  
  
  for (i in 1:1) {
    median_value <- median(eqtl[,i], na.rm = TRUE)
    eqtl[is.na(eqtl[,i]),i] <- median_value
  }
  
  
  referprior <- merge(referchr, eqtl, by = "bp_hg38", all.x = TRUE)
  referprior[, ncol(referprior)][is.na(referprior[, ncol(referprior)])] <- median(referprior[, ncol(referprior)], na.rm = TRUE)
  
  
  bmichr <- sum_stats[sum_stats$chr==chrid,]
  dim(bmichr)
  
  unique_indices <- which(!duplicated(bmichr$bp))
  bmichr <- bmichr[unique_indices,]
  dim(bmichr)
  
  path <- paste("/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/UKBB_panel/hg38/chr/ukb_gtex19_chr",chrid,".bim", sep="")
  ukbbim19chr <- read.table(path)
  
  print(paste('all(ukbbim19chr$V4 %in% bmichr$bp)',all(ukbbim19chr$V4 %in% bmichr$bp)))
  
  
  index <- match(ukbbim19chr$V4, bmichr$bp)
  table(is.na(index))
  bmichrcommon <- bmichr[index,]
  table(bmichrcommon$bp==ukbbim19chr$V4)
  
  #length(which(bmichrcommon$pval<0.000000005))
  
  table(bmichrcommon$a1==ukbbim19chr$V5 & bmichrcommon$a2==ukbbim19chr$V6)
  table(bmichrcommon$a1==ukbbim19chr$V6 & bmichrcommon$a2==ukbbim19chr$V5)
  
  inv <- which(bmichrcommon$a1==ukbbim19chr$V6 & bmichrcommon$a2==ukbbim19chr$V5)
  
  for(i in inv){
    bmichrcommon[i,10] <- -bmichrcommon[i,10]
    a <- bmichrcommon[i,14]
    bmichrcommon[i,14] <- bmichrcommon[i,15]
    bmichrcommon[i,15] <- a
  }
  # solve the problem of allele coding errors
  
  head(bmichrcommon)
  head(ukbbim19chr)
  
  print(table(bmichrcommon$a1==ukbbim19chr$V5 & bmichrcommon$a2==ukbbim19chr$V6))
  
  
  tail(bmichrcommon)
  tail(ukbbim19chr)
  
  # read in the reference panel:
  
  
  
  
  
  path <- paste("/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/UKBB_panel/hg38/chr/ukb_gtex19_chr",chrid, sep="")
  #path <- '/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/UKBB_panel/hg38/chr/ukb_gtex19_chr2'
  library(genio)
  dat <- read_plink(path)
  #X <- bedNA(t(dat$X))
  #X <- scale(X)
  
  print(paste("table(bmichrcommon$bp==dat$bim$pos):",table(bmichrcommon$bp==dat$bim$pos)))
  table(bmichrcommon$bp==dat$bim$pos) # check if this matches for all the SNPs
  
  
  
  
  risk <- which(bmichrcommon$pval<0.0000005)
  for(width in c(50000)){
   
    library(susieR)
    
    riskregion <- matrix(0,nrow=0,ncol=2)
    remainrisk <- risk
    snp_byregion <- list()
    
    risk_bp <- bmichrcommon$bp[risk]
    chr_start <- bmichrcommon$bp[1]
    chr_end <- bmichrcommon$bp[length(bmichrcommon$bp)]
    risk_diff <- diff(risk_bp)
    where_sep <- which(risk_diff>2*width)
    if(length(where_sep)>0){
      region1 <- c(max(chr_start, risk_bp[1]-width), min(chr_end,risk_bp[where_sep[1]]+width))
      riskregion <- rbind(riskregion, as.numeric(region1))
      if(length(where_sep)>1){
        for(sepidx in 2:length(where_sep)){
          region1 <- c(max(chr_start, risk_bp[where_sep[sepidx-1]+1]-width), min(chr_end,risk_bp[where_sep[sepidx]]+width))
          riskregion <- rbind(riskregion, as.numeric(region1))
        }
      }
    }
    
    
 
    
    print(dim(riskregion))
    # get the risk region and all the SNPs in this risk region
    
    
    
    result0 <- list()
    resulte <- list()
    bphg190 <- list()
    bphg19e <- list()
    purity0 <- numeric()
    puritye <- numeric()
    find_signal0 <- numeric()
    find_signale <- numeric()
    
    for(j in 1:nrow(riskregion)){
      region <- riskregion[j,]
      index <- which(bmichrcommon$bp>region[1] & bmichrcommon$bp<region[2])
      summary_j <- bmichrcommon[index, c(5,10,13)]
      X_j <- dat$X[index,]
      X_j <- bedNA(t(X_j))
      X_j <- scale(X_j)
      correlationgenet <- cor(X_j)
      bp_snp <- bmichrcommon$bp[index]
      index_prior <- match(bp_snp, referprior$bp_hg19)
      prior <- referprior$prior_wb[index_prior]
      print(estimate_s_rss(summary_j$tstat, correlationgenet, n=mean(summary_j$n_complete_samples)))
      fit <- susie_rss(z=summary_j$tstat, min_abs_corr = 0.90, 
                       R = correlationgenet, n = mean(summary_j$n_complete_samples), L=5, estimate_prior_variance = TRUE)
      fite <- susie_rss(z=summary_j$tstat, min_abs_corr = 0.90, prior_weights = prior,
                        R = correlationgenet, n = mean(summary_j$n_complete_samples), L=5, estimate_prior_variance = TRUE)
      result0 <- append(result0, fit$sets$cs)
      resulte <- append(resulte, fite$sets$cs)
      
      print(paste("iteration:",j))
      print(length(fit$sets$cs))
      print(length(fite$sets$cs))
      
      print(sapply(fit$sets$cs, length))
      print(sapply(fite$sets$cs, length))
      
      if(length(fit$sets$cs)>0){
        purity0 <- rbind(purity0, matrix(unlist(fit$sets$purity), ncol=3))
      }
      if(length(fite$sets$cs)>0){
        puritye <- rbind(puritye, matrix(unlist(fite$sets$purity), ncol=3))
      }
      
      bphg190 <- append(bphg190, lapply(fit$sets$cs, function(x) bp_snp[x]))
      bphg19e <- append(bphg19e, lapply(fite$sets$cs, function(x) bp_snp[x]))
      if(length(fit$sets$cs)>0){
        find_signal0 <- append(find_signal0, j)
      }
      if(length(fite$sets$cs)>0){
        find_signale <- append(find_signale, j)
      }
      
      #  print(paste(j,":",length(fite$sets$cs)))
    }
  
    
    
    print(colMeans(purity0))
    print(colMeans(puritye))
    
    cslength <- append(cslength,sapply(result0, length))
    cslength2 <- append(cslength2,sapply(resulte, length))
    
    susier <- c(length(result0),mean(sapply(result0, length)),length(which(sapply(result0, length)==1)),length(which(sapply(result0, length)<5)))
    susie2 <- c(length(resulte),mean(sapply(resulte, length)),length(which(sapply(resulte, length)==1)),length(which(sapply(resulte, length)<5)))
    
    total_result <- rbind(susier, susie2)
    colnames(total_result) <- c("total","ave","<=1","<5")
    print(total_result)
    
  }

  
  bphg190[which(sapply(result0, length)==1)]
  bphg19e[which(sapply(resulte, length)==1)]

  path <- paste("/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/BMI/result/chr",chrid,"/totalads_5w0.9check.txt", sep="")
  write.table(total_result, path, 
            quote=F, row.names = T, col.names = T)
  
  #path <- paste("/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/BMI/result/chr",chrid,"/riskregionads.txt", sep="")
  #write.table(riskregion, path,
  #            quote = F, row.names = F, col.names = F)
  
  #path <- paste("/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/BMI/result/chr",chrid,"/referpriorads.txt", sep="")
  #write.table(referprior, path,
  #            quote = F, row.names = F)
  
#  save(result0, resulte, find_signal0, find_signale, bphg190, bphg19e, file="/gpfs/gibbs/pi/zhao/xz527/susiereview/real_data/BMI/result/chr5/resultwb.RData")
}
