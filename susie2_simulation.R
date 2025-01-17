args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
args <- as.numeric(args)
# get the parameters from the command line

library(mvsusieR)
library(susieR)


bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm=T)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


path <- '/gpfs/gibbs/pi/zhao/xz527/susiereview/new_version_ukb/ukb'
# path <- '~/data/UKBB/chr1eg2q'
library(genio)
dat <- read_plink(path)

glen <- 5000
starting <- sample(1:(93246+1-glen), 1)
samplesnp <- starting:(starting+glen-1)

library(susieR)

fam <- dat$fam
bim <- dat$bim[samplesnp,]
zfile <- bim[,c(2,1,4,5,6)]
colnames(zfile) <- c("rsid","chromosome","position","allele1","allele2")
zfile$maf <- pmin(rowMeans(dat$X[samplesnp,], na.rm = T)/2,1-rowMeans(dat$X[samplesnp,], na.rm = T)/2)
X <- bedNA(t(dat$X[samplesnp,]))
X <- scale(X)

n <- nrow(X)
p <- ncol(X)
d <- 1
k <- 4
h2e <- args[1] # heritability of eQTL
h2p <- args[1] # heritability of phenotype
# h2e <- h2p <- 0.1
np <- 2
gamma <- 1
pre <- 0.5 # proportion of heritability explained by common SNPs
prp <- 0.5
ne <- 4/k
nc <- 4/k
iteration <- args[2]
dis <- 100
leng <- 100

h2 <- (gamma^2*h2e*k+h2p)/(gamma^2*k+1) # check the heritability
write.table(c(h2,10),file=paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/simulation/insample/result/herieg",".txt"),quote=F)


prop1 <- numeric(d) # proportion of casual variables detected
prop2 <- numeric(d) # proportion of credible sets that contains at least one casual variable
len <- numeric(d) # number of credible sets
ave <- numeric(d) # averaged length of credible sets
prop3 <- numeric(d)
prop1p <- numeric(d)
prop2p <- numeric(d)
lenp <- numeric(d)
avep <- numeric(d)
prop3p <- numeric(d)
prop1pa1 <- numeric(d)
prop2pa1 <- numeric(d)
lenpa1 <- numeric(d)
avepa1 <- numeric(d)
prop3pa1 <- numeric(d)
prop1mv <- numeric(d) # proportion of casual variables detected
prop2mv <- numeric(d) # proportion of credible sets that contains at least one casual variable
lenmv <- numeric(d) # number of credible sets
avemv <- numeric(d) # averaged length of credible sets
prop3mv <- numeric(d)
prop1mv1 <- numeric(d) # proportion of casual variables detected
prop2mv1 <- numeric(d) # proportion of credible sets that contains at least one casual variable
lenmv1 <- numeric(d) # number of credible sets
avemv1 <- numeric(d) # averaged length of credible sets
prop3mv1 <- numeric(d)
prop1fl <- numeric(d) # proportion of casual variables detected
prop2fl <- numeric(d) # proportion of credible sets that contains at least one casual variable
lenfl <- numeric(d) # number of credible sets
avefl <- numeric(d) # averaged length of credible sets
prop3fl <- numeric(d)





for (i in 1:d) {
  
  ## ----- simulate the phenotype with specified heritability -------
  
  
  s <- sample((dis+1):(p-leng*k-dis*(2*k+0)),k)
  s <- s[order(s)]
  s <- s+((1:k)-1)*(leng+dis*2)
  
  betae <- matrix(0, nrow = p, ncol = k)
  ye <- matrix(0, nrow = n, ncol = k)
  indp <- sample(p,np)
  inde <- indc <- numeric(0)
  for (j in 1:k) {
    inde1 <- sample(setdiff(s[j]:(s[j]+leng),indp),ne)
    indc1 <- sample(setdiff(s[j]:(s[j]+leng),c(inde1,indp)),nc)
    betae[inde1,j] <- rnorm(length(inde1),0,sqrt(h2e*(1-pre)/length(inde1)))
    betae[indc1,j] <- rnorm(length(indc1),0,sqrt(h2e*pre/length(indc1)))
    ye[,j] <- X%*%betae[,j]+rnorm(n,0,sqrt(1-h2e))
    inde <- append(inde,inde1)
    indc <- append(indc,indc1)
  }
  
  betap <- rep(0,p)
  #print(indp)
  betap[indp] <- rnorm(length(indp),0,sqrt(h2p*(1-prp)/length(indp)))
  betap[indc] <- rnorm(length(indc),0,sqrt(h2p*prp/length(indc)))
  yp <- X%*%betap+rnorm(n,0,sqrt(1-h2p))
  gammaall <- rep(gamma,k)
  y <- (yp+ye%*%gammaall)/sqrt(1+sum(gammaall^2))
  ## ---------------------------------------------------
  
  
  yeh <- matrix(0, nrow = n, ncol = k/2)
  
  for(j in 1:(k/2)){
    yeh[,j] <- ye[,2*j-1]+ye[,2*j]
    yeh[,j] <- scale(yeh[,j])
  }
  
  
  ## ----- generate the summary statistics  -------
  
  
  beta <- numeric(p)
  se <- numeric(p)
  tvalue <- numeric(p)
  for(l in 1:p){
    linearmodel <- lm(y ~ X[,l])
    beta[l] <- coef(summary(linearmodel))[2,1]
    se[l] <- coef(summary(linearmodel))[2,2]
    tvalue[l] <- coef(summary(linearmodel))[2,3]
  }
  
  datasetz <- vector("list",k/2+1)
  
  zlist <- zfile
  zlist$beta <- beta
  zlist$se <- se
  datasetz[[1]] <- zlist
  
  
  ## ---------------------------------------------------
  
  
  
  ## ----- conduct the SuSiE and SuSiE2  -----------------
  
  correlationgenet <- cor(X)

  lambda = estimate_s_rss(beta/se, correlationgenet, n=10000) # check the balance between summary statistics and LD matrix
  
  
  
  fit <- susie_rss(bhat = beta, shat = se, min_abs_corr = 0.5,
                   R = correlationgenet, n = n, L=10, estimate_prior_variance = TRUE)  # run the SuSiE
  
  cs <- fit$sets$cs
  out <- as.list(cs)
  
  prior3 <- rep(np/p, p)
  prior4 <- rep(np/p, p)
  
  
  ## -------- generate the prior used in SuSiE2 ---------------

  for (j in 1:k) {
    X1 <- X[,(s[j]-dis+1):(s[j]+leng+dis)]
    fite1 <- susie(X1,ye[,j],L=(ne+nc),standardize = FALSE,
                   estimate_residual_variance = FALSE,
                   scaled_prior_variance = 0.1,
                   min_abs_corr = 0.5)
    priore1 <- numeric(ncol(fite1$alpha))
    for(jj in 1:ncol(fite1$alpha)){
      priore1[jj] <- 1-prod(1-fite1$alpha[,jj])
    }
    #priore1 <- colSums(fite1$alpha)
    prior3[(s[j]-dis+1):(s[j]+leng+dis)] <- priore1
  }
  
  prior4[(s[k]-dis+1):(s[k]+leng+dis)] <- priore1
  
  ## ------------------------------------------------------------
  

  
  
  fitp1 <- susie_rss(bhat = beta, shat = se, prior_weights = prior4, min_abs_corr = 0.5,
                     R = correlationgenet, n = n, L=10, estimate_prior_variance = TRUE)
  csp1 <- fitp1$sets$cs
  outp1 <- as.list(csp1)
  # SuSiE2 only use a part of prior information

  
  fitpa <- susie_rss(bhat = beta, shat = se, prior_weights = prior3, min_abs_corr = 0.5,
                     R = correlationgenet, n = n, L=10, estimate_prior_variance = TRUE)
  cspa <- fitpa$sets$cs
  outpa <- as.list(cspa)
  # SuSiE2 use all the prior information
  
  
  uout <- unlist(out)
  uoutp1 <- unlist(outp1)
  uoutpa <- unlist(outpa)
  
  ## ------------------------------------------------------------
  
  
  
  
  ## ----- conduct the flashfm and mvSuSiE  ----------------
  
  betam <- matrix(0, nrow=p, ncol=k+1)
  sem <- matrix(0, nrow=p, ncol=k+1)
  tvaluem <- matrix(0, nrow=p, ncol=k+1)
  betam[,1] <- beta
  sem[,1] <- se
  tvaluem[,1] <- tvalue
  for(tt in 1:k){
    for(l in 1:p){
      linearmodel <- lm(ye[,tt] ~ X[,l])
      betam[l,tt+1] <- coef(summary(linearmodel))[2,1]
      sem[l,tt+1] <- coef(summary(linearmodel))[2,2]
      tvaluem[l,tt+1] <- coef(summary(linearmodel))[2,3]
    }
  }
  
  
  betam <- matrix(0, nrow=p, ncol=k/2+1)
  sem <- matrix(0, nrow=p, ncol=k/2+1)
  tvaluem <- matrix(0, nrow=p, ncol=k/2+1)
  betam[,1] <- beta
  sem[,1] <- se
  tvaluem[,1] <- tvalue
  for(tt in 1:(k/2)){
    for(l in 1:p){
      linearmodel <- lm(yeh[,tt] ~ X[,l])
      betam[l,tt+1] <- coef(summary(linearmodel))[2,1]
      sem[l,tt+1] <- coef(summary(linearmodel))[2,2]
      tvaluem[l,tt+1] <- coef(summary(linearmodel))[2,3]
    }
    zlist <- zfile
    zlist$beta <- betam[,tt+1]
    zlist$se <- sem[,tt+1]
    datasetz[[tt+1]] <- zlist
  }
  
  multi <- cbind(y, ye)
  z <- susieR:::calc_z(X, multi)
  
  
  ind <- append(append(inde,indp),indc)
  
  library(flashfm)
  raf <- numeric(p)+mean(zfile$maf)
  multi <- cbind(y, yeh)
  covY <- cor(multi)
  ybar <- numeric(k/2+1)
  N <- rep(n,k/2+1)
  names(raf) <- zfile$rsid
  path <- paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/finemap_v1.4.2_x86_64/result/iter_",iteration,"/region")
  fm2 <- FLASHFMwithFINEMAP(gwas.list=datasetz,corX=correlationgenet,raf=raf,ybar=ybar,N=N,fstub=path,
                            TOdds=1,covY,cpp=.95,NCORES=5,FMpath="/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64") 
  
  statement <- rownames(fm2$mpp.pp$PPg[[1]])[1]
  
  # Split the statement into individual component names
  component_names <- strsplit(statement, "%")[[1]]
  flashlist <- fm2$snpGroups$groups.flashfm[component_names]
  rsind <- zfile$rsid[ind]
  
  prop1fl[i] <- length(intersect(unlist(flashlist), rsind))/length(rsind)
  a <- lapply(flashlist, function(x) length(intersect(rsind,x)))
  prop2fl[i] <- sum(unlist(a)>0)/length(flashlist)
  # proportion of credible sets that contains at least one casual variable
  
  prop3fl[i] <- length(setdiff(unlist(flashlist),rsind))/(p-length(rsind))
  
  lenfl[i] <- length(as.numeric(lapply(flashlist, length)))
  # number of credible sets
  
  avefl[i] <- mean(as.numeric(lapply(flashlist, length)))
  
  ## ------------------------------------------------------------
  
  

  
  path <- paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/finemap_v1.4.2_x86_64/result/iter_",iteration,"/result",h2p,"_",iteration,".RData")
  save(fm2, file= path)
  
  
  
  
  ## ---------------save the input files for fastpaintor------------------------
  
  colnames(z) <- c("Z1","Z2","Z3","Z4","Z5")
  
  file=paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/simulation/insample/Inpu/Locus4k10",h2p,'_',iteration)
  write.table(z, file,
              row.names = F, col.names = T, quote = F)
  
  file=paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/simulation/insample/Inpu/Locus4k10",h2p,'_',iteration,".LD1")
  write.table(correlationgenet, file,
              row.names = F, col.names = F, quote = F)
  
  anno <- numeric(p)+1
  file=paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/simulation/insample/Inpu/Locus4k10",h2p,'_',iteration,".annotations")
  write.table(anno, file,
              row.names = F, col.names = T, quote = F)
  
  
  ## ------------------------------------------------------------
  

  
  ## --------------save the result------------------------------
  
  prior_var <- create_mixture_prior(R=k+1)
  mvss <- mvsusie_rss(z, R = correlationgenet, N = n, L=10, prior_variance = prior_var, 
                      estimate_prior_variance = TRUE, min_abs_corr = 0.5)
  
  csmv <- mvss$sets$cs
  outmv <- as.list(csmv)
  uoutmv <- unlist(outmv)
  
  csmv1 <- mvss$sets$cs
  outmv1 <- as.list(csmv1)
  lfsr <- mvss$single_effect_lfsr[,1]
  rem_id <- which(lfsr>0.05)
  l_id <- paste('L',rem_id, sep='')
  outmv1[l_id] <- NULL
  uoutmv1 <- unlist(outmv1)
  
  
  
  ind <- append(append(inde,indp),indc)
  
  
  file=paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/simulation/insample/result/ind4k10","_",h2p,"_",iteration,".txt")
  write.table(ind, file, quote=F, col.names = F, row.names = F)
  
  rsind <- zfile$rsid[ind]
  file=paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/simulation/insample/result/rsind4k10","_",h2p,"_",iteration,".txt")
  write.table(rsind, file, quote=F, col.names = F, row.names = F)
  
  
  prop1[i] <- length(intersect(as.numeric(ind),as.numeric(uout)))/length(ind)
  prop1p[i] <- length(intersect(as.numeric(ind),as.numeric(uoutp1)))/length(ind)
  prop1pa1[i] <- length(intersect(as.numeric(ind),as.numeric(uoutpa)))/length(ind)
  prop1mv[i] <- length(intersect(as.numeric(ind),as.numeric(uoutmv)))/length(ind)
  prop1mv1[i] <- length(intersect(as.numeric(ind),as.numeric(uoutmv1)))/length(ind)
  # proportion of casual variables detected
  
  a <- lapply(out, function(x) length(intersect(as.numeric(ind),as.numeric(x))))
  prop2[i] <- sum(unlist(a)>0)/length(out)
  a <- lapply(outp1, function(x) length(intersect(as.numeric(ind),as.numeric(x))))
  prop2p[i] <- sum(unlist(a)>0)/length(outp1)
  a <- lapply(outpa, function(x) length(intersect(as.numeric(ind),as.numeric(x))))
  prop2pa1[i] <- sum(unlist(a)>0)/length(outpa)
  a <- lapply(outmv, function(x) length(intersect(as.numeric(ind),as.numeric(x))))
  prop2mv[i] <- sum(unlist(a)>0)/length(outmv)
  a <- lapply(outmv1, function(x) length(intersect(as.numeric(ind),as.numeric(x))))
  prop2mv1[i] <- sum(unlist(a)>0)/length(outmv1)
  # proportion of credible sets that contains at least one casual variable
  
  prop3[i] <- length(setdiff(as.numeric(uout),as.numeric(ind)))/(p-length(ind))
  prop3p[i] <- length(setdiff(as.numeric(uoutp1),as.numeric(ind)))/(p-length(ind))
  prop3pa1[i] <- length(setdiff(as.numeric(uoutpa),as.numeric(ind)))/(p-length(ind))
  prop3mv[i] <- length(setdiff(as.numeric(uoutmv),as.numeric(ind)))/(p-length(ind))
  prop3mv1[i] <- length(setdiff(as.numeric(uoutmv1),as.numeric(ind)))/(p-length(ind))
  
  
  len[i] <- length(as.numeric(lapply(out, length)))
  lenp[i] <- length(as.numeric(lapply(outp1, length)))
  lenpa1[i] <- length(as.numeric(lapply(outpa, length)))
  lenmv[i] <- length(as.numeric(lapply(outmv, length)))
  lenmv1[i] <- length(as.numeric(lapply(outmv1, length)))
  # number of credible sets
  
  
  ave[i] <- mean(as.numeric(lapply(out, length)))
  avep[i] <- mean(as.numeric(lapply(outp1, length)))
  avepa1[i] <- mean(as.numeric(lapply(outpa, length)))
  avemv[i] <- mean(as.numeric(lapply(outmv, length)))
  avemv1[i] <- mean(as.numeric(lapply(outmv1, length)))
  # averaged length of credible sets
  
  
  # write.table(i,"~/data/UKBB/check.txt",quote = F)
}



result<-cbind(prop1,prop2,ave,len,prop3,prop1p,prop2p,avep,lenp,prop3p,prop1pa1,prop2pa1,avepa1,lenpa1,prop3pa1,
              prop1mv,prop2mv,avemv,lenmv,prop3mv,prop1mv1,prop2mv1,avemv1,lenmv1,prop3mv1,prop1fl,prop2fl,avefl,lenfl,prop3fl)

file=paste0("/gpfs/gibbs/pi/zhao/xz527/susiereview/fastpaintor/simulation/insample/result/result4k100.5","_",h2p,"_",iteration,".txt")
write.table(result,file, quote = F)

## ------------------------------------------------------------
