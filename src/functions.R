marginalEffSim_sixLoci <- function(data = data_151002, 
                                   snps = c("4944267_chrVIII_114337_G_A", "9716860_chrXIV_469224_A_G", "1242017_chrIII_198615_T_G", "2358650_chrIV_998628_A_T", "8733525_chrXIII_410320_T_C", "7890567_chrXII_645539_A_G"), 
                                   freqs = c(.2, .8), 
                                   trait = "trait_8.mean"){
  #This function simulates populations as described in the paper "Accounting for genetic interactions improves modeling of individual quantitative trait phenotypes in yeast"
  #The default arguments are just there as an example
  
  #Arguments:
  #data = An object of gwaa.data-class (implemented in the GenABEL package)
  #snps = Index, character or logical vector with subset of SNPs in data. The vector has to define 6 SNPs
  #freqs = The allele frequencies at which to simulate populations. WARNING: all length(freqs)^6 possible populations will be simulated, causing the running time to grow rapidly with length(freqs)
  #trait = A character string indicating a phenotype in data
  
  if(!require(GenABEL)){
    stop("You need GenABEL")
  }
  
  geno <- as.double.gwaa.data(data[, snps])
  geno[geno == 2] <- 1
  class.factor <- factor(paste(geno[,1], geno[,2], geno[,3], geno[,4], geno[,5], geno[,6], sep = "")) #, levels = c("00", "01", "10", "11"))
  classes <- sapply(levels(class.factor), function(x){as.numeric(strsplit(x, split = "")[[1]])}) #The possible genotype classes. 64 in this case
  y.split <- split(x = data@phdata[,trait], f = class.factor)
  y.split.means <- sapply(X = y.split, FUN = mean, na.rm = T)
  #y.split.sd <- sapply(X = y.split, FUN = sd, na.rm = T)
  y.split.sd <- sd(data@phdata[,trait], na.rm = T) #Simulate using the TOTAL variance when sampling from every genotype class
  N <- sum(!is.na(data@phdata[,trait]))
  #the elements in this output-list needs to be in the same order as the columns in geno
  marginalSim <- setNames(list(list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)),
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6))), 
                          snps)
  j <- 1
  k <- 1
  pb <- txtProgressBar(style=3)
  ptm <- proc.time()
  for(p1 in freqs){
    for(p2 in freqs){
      for(p3 in freqs){
        for(p4 in freqs){
          for(p5 in freqs){
            for(p6 in freqs){
              p <- c(p1, p2, p3, p4, p5, p6)
              q <- 1-p
              classes.p <- apply(classes, 2, FUN = function(x){prod(p[x == 0], q[x == 1])}) #The frequency of every genotype class
              x <- sample(x = 1:ncol(classes), size = N, prob = classes.p, replace = T) #number of observations in every class
              
              for(i in 1:6){
                a <- which(classes[i,] == 1) #the classes that makes up the first marginal class
                a.means <- y.split.means[a]
                a.sd <- y.split.sd[a]
                a.n <- table(x)[as.character(a)] #the number of observations in every class that makes up a
                a.n[is.na(a.n)] <- 0
                #a.pheno <- rnorm(n = sum(x %in% a), mean = rep(a.means, a.n), sd = rep(a.sd, a.n))
                a.pheno <- rnorm(n = sum(x %in% a), mean = rep(a.means, a.n), sd = y.split.sd)
                
                b <- which(classes[i,] == 0) #the classes that makes up the second marginal class
                b.means <- y.split.means[b]
                b.sd <- y.split.sd[b]
                b.n <- table(x)[as.character(b)] #the number of observations in every class that makes up b
                b.n[is.na(b.n)] <- 0
                #b.pheno <- rnorm(n = sum(x %in% b), mean = rep(b.means, b.n), sd = rep(b.sd, b.n))
                b.pheno <- rnorm(n = sum(x %in% b), mean = rep(b.means, b.n), sd = y.split.sd)
                
                pheno <- c(a.pheno, b.pheno)
                geno <- c(rep(1, length(a.pheno)), rep(0, length(b.pheno)))
                
                #Fit the model. This "manual" way is somewhat faster than using lm()
                #NULL
                X0 <- rep(1, length(geno))
                qr0 <- qr(X0)
                est0 <- qr.coef(qr0, pheno)
                res0 <- pheno - X0 * est0
                RSS.0 <- sum(res0^2)/N
                #Model with one SNP 
                X1 <- cbind(rep(1, length(geno)), geno)
                qr1 <- qr(X1)
                est1 <- qr.coef(qr1, pheno)
                res1 <- pheno - X1 %*% est1
                RSS.1 <- sum(res1^2)/N
                LRT <- -N*( log(RSS.1)-log(RSS.0) )
                #Model with all SNPs
                
                marginalSim[[i]]$effect[j] <- est1[2]
                marginalSim[[i]]$p[j] <- 1 - pchisq(LRT, df=1) #Needs to be changed if more than 1 parameter is tested
                marginalSim[[i]]$R2[j] <- 1 - RSS.1/RSS.0
                names(marginalSim[[i]]$effect)[j] <- names(marginalSim[[i]]$p)[j] <- names(marginalSim[[i]]$R2)[j] <- paste(p, collapse = "_")
              }
              j <- j+1
            }
          }
        }
      }
    }
    setTxtProgressBar(pb, k/length(freqs))
    k <- k+1
  }
  proc.time() - ptm
  close(pb)
  return(marginalSim)
}

marginalEffSim_sixLoci_add <- function(data = data_151002, 
                                        snps = c("4944267_chrVIII_114337_G_A", "9716860_chrXIV_469224_A_G", "1242017_chrIII_198615_T_G", "2358650_chrIV_998628_A_T", "8733525_chrXIII_410320_T_C", "7890567_chrXII_645539_A_G"), 
                                        freqs = c(.2, .8), 
                                        trait = "trait_8.mean"){
  #This function simulates populations as described in the paper "Accounting for genetic interactions improves modeling of individual quantitative trait phenotypes in yeast"
  #The phenotype of every individual is given AS IF the six loci acted completely additively
  #The default arguments are just there as an example
  
  #Arguments:
  #data = An object of gwaa.data-class (implemented in the GenABEL package)
  #snps = Index, character or logical vector with subset of SNPs in data. The vector has to define 6 SNPs
  #freqs = The allele frequencies at which to simulate populations. WARNING: all length(freqs)^6 possible populations will be simulated, causing the time running time to grow rapidly with length(freqs)
  #trait = A character string indicating a phenotype in data
  
  if(!require(GenABEL)){
    stop("You need GenABEL")
  }
  
  havePheno <- !is.na(data@phdata[, trait])
  geno <- as.double.gwaa.data(data[havePheno, snps])
  geno[geno == 2] <- 1
  addModel <- lm(data@phdata[havePheno, trait] ~ geno[,1] + geno[,2] + geno[,3] + geno[,4] + geno[,5] + geno[,6])
  
  class.factor <- factor(paste(geno[,1], geno[,2], geno[,3], geno[,4], geno[,5], geno[,6], sep = "")) #, levels = c("00", "01", "10", "11"))
  classes <- sapply(levels(class.factor), function(x){as.numeric(strsplit(x, split = "")[[1]])}) #The possible genotype classes. 64 in this case
  y.split <- split(x = addModel$fitted.values, f = class.factor) #The predicted value in every class
  y.split.means <- sapply(X = y.split, FUN = mean, na.rm = T) 
  #y.split.sd <- sapply(X = y.split, FUN = sd, na.rm = T)
  #y.split.sd <- sd(addModel$residuals) #Simulate using the RESIDUAL variance when sampling from every genotype class
  y.split.sd <- sd(data@phdata[,trait], na.rm = T) #Simulate using the TOTAL variance when sampling from every genotype class
  #y.split.sd <- sapply(X = split(x = data@phdata[havePheno, trait], f = class.factor), FUN = sd, na.rm = T) #Use the real sd within every class
  N <- sum(havePheno)
  #the elements in this output-list needs to be in the same order as the columns in geno
  marginalSim <- setNames(list(list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)),
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6)), 
                               list(effect = rep(NA, length(freqs)^6), R2 = rep(NA, length(freqs)^6), p = rep(NA, length(freqs)^6))), 
                          snps)
  j <- 1
  k <- 1
  pb <- txtProgressBar(style=3)
  ptm <- proc.time()
  for(p1 in freqs){
    for(p2 in freqs){
      for(p3 in freqs){
        for(p4 in freqs){
          for(p5 in freqs){
            for(p6 in freqs){
              p <- c(p1, p2, p3, p4, p5, p6)
              q <- 1-p
              classes.p <- apply(classes, 2, FUN = function(x){prod(p[x == 0], q[x == 1])}) #The frequency of every genotype class
              x <- sample(x = 1:ncol(classes), size = N, prob = classes.p, replace = T) #number of observations in every class
              
              for(i in 1:6){
                a <- which(classes[i,] == 1) #the classes that makes up the first marginal class
                a.means <- y.split.means[a]
                a.sd <- y.split.sd[a]
                a.n <- table(x)[as.character(a)] #the number of observations in every class that makes up a
                a.n[is.na(a.n)] <- 0
                a.pheno <- rnorm(n = sum(x %in% a), mean = rep(a.means, a.n), sd = y.split.sd)
                #a.pheno <- rnorm(n = sum(x %in% a), mean = rep(a.means, a.n), sd = rep(a.sd, a.n))
                
                b <- which(classes[i,] == 0) #the classes that makes up the second marginal class
                b.means <- y.split.means[b]
                b.sd <- y.split.sd[b]
                b.n <- table(x)[as.character(b)] #the number of observations in every class that makes up b
                b.n[is.na(b.n)] <- 0
                b.pheno <- rnorm(n = sum(x %in% b), mean = rep(b.means, b.n), sd = y.split.sd)
                #b.pheno <- rnorm(n = sum(x %in% b), mean = rep(b.means, b.n), sd = rep(b.sd, b.n))
                
                pheno <- c(a.pheno, b.pheno)
                geno <- c(rep(1, length(a.pheno)), rep(0, length(b.pheno)))
                
                #Fit the model
                #NULL
                X0 <- rep(1, length(geno))
                qr0 <- qr(X0)
                est0 <- qr.coef(qr0, pheno)
                res0 <- pheno - X0 * est0
                RSS.0 <- sum(res0^2)/N
                #Model with one SNP 
                X1 <- cbind(rep(1, length(geno)), geno)
                qr1 <- qr(X1)
                est1 <- qr.coef(qr1, pheno)
                res1 <- pheno - X1 %*% est1
                RSS.1 <- sum(res1^2)/N
                LRT <- -N*( log(RSS.1)-log(RSS.0) )
                #Model with all SNPs
                
                marginalSim[[i]]$effect[j] <- est1[2]
                marginalSim[[i]]$p[j] <- 1 - pchisq(LRT, df=1) #Needs to be changed if more than 1 parameter is tested
                marginalSim[[i]]$R2[j] <- 1 - RSS.1/RSS.0
                names(marginalSim[[i]]$effect)[j] <- names(marginalSim[[i]]$p)[j] <- names(marginalSim[[i]]$R2)[j] <- paste(p, collapse = "_")
              }
              j <- j+1
            }
          }
        }
      }
    }
    setTxtProgressBar(pb, k/length(freqs))
    k <- k+1
  }
  proc.time() - ptm
  close(pb)
  return(marginalSim)
}

vGWAS.dglm <- function(data, trait.name = NULL, pheno = NULL, coVars = NULL, returnModels = F, returnMap = F, ...){
  #This function performs a vGWAS analysis using a double generalized linear model. The model fitting is performed by the dglm function in the package dglm
  
  #Arguments:
  #data = An object of gwaa.data-class (implemented in the GenABEL package)
  #trait.name = A character string indicating a phenotype in data
  #pheno = A phenotype vector. Either trait.name or pheno must be specified
  #coVars = Index, character or logical vector with subset of SNPs in data. These SNPs will be included as fixed effects in the mean part of the model
  #returnModels = If true, a list will be returned containing the fitted dglm models for each SNP (Default: only p-values of the dispersion effect is returned)
  #returnMap = = If true, a list will be returned containing the position, chromosome, and p-value for the dispersion effect of every SNP (Default: only p-values of the dispersion effect is returned)
  if(!require(dglm)){
    stop("Could not load package dglm")
  }
  if(!is.null(trait.name)){
    pheno <- data@phdata[,trait.name]
  }
  else if(is.null(pheno)){
    stop("trait.name or pheno must be specified")
  }
  
  if(returnModels){
    output <- list(models = list(), p = c())
  }
  else if(returnMap){
    output <- list(map = c(), chr = c(), p = c())
  }
  else{
    output <- rep(x = NA, times = length(data@gtdata@snpnames))
  }
  i <- 1
  pb <- txtProgressBar(style=3)
  for(SNP in data@gtdata@snpnames){
    if(!is.null(coVars)){
      geno <- as.double.gwaa.data(data[,c(coVars, SNP)])
      data.dglm <- na.omit(as.data.frame(cbind(pheno, geno)))
      formula.mean <- formula(paste("pheno", "~", paste(coVars, collapse = "+"), sep = "")) #The formula for the mean model, with covariates
    }
    else{
      geno <- as.double.gwaa.data(data[, SNP])
      data.dglm <- na.omit(data.frame(pheno = pheno, geno = geno[,1]))
      #data.dglm <- na.omit(as.data.frame(cbind(pheno, geno)))
      formula.mean <- formula(paste("pheno", "~", "1", sep = "")) #The formula for the mean model, without covariates
    }
    
    formula.disp <- formula(paste("~", "geno")) #The formula for the dispersion model
    #formula.disp <- formula(paste("~", SNP)) #The formula for the dispersion model
    sink("/dev/null") #supress printing from the dglm-function
    model <- dglm(formula = formula.mean, dformula = formula.disp, data = data.dglm, ...)
    sink()
    
    #Wrap up the output
    if(returnModels){
      output$models[[i]] <- model 
      output$p[i] <- summary(model$dispersion.fit)$coefficients[SNP, 4] #the p-value of SNP in the dispersion model
    }
    else if(returnMap){
      output$map[i] <- data@gtdata@map[SNP]
      output$chr[i] <- as.numeric(as.character(data@gtdata@chromosome[SNP]))
      #output$p[i] <- summary(model$dispersion.fit)$coefficients[SNP, 4] #the p-value of SNP in the dispersion model
      output$p[i] <- summary(model$dispersion.fit)$coefficients["geno", 4] #the p-value of SNP in the dispersion model
    }
    else{
      output[i] <- summary(model$dispersion.fit)$coefficients[SNP, 4] #the p-value of SNP in the dispersion model
    }
    
    i <- i+1
    setTxtProgressBar(pb, i/data@gtdata@nsnps)
  }
  close(pb)
  return(output)
}

