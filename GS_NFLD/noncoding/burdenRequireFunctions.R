CNVBurdenTest <- function(cnv.matrix, geneset, label, covariates, correctGlobalBurden = T, standardizeCoefficient = T, 
                          permutation = T, nperm = 100, BiasedUrn = F){
  
  model = "lm"
  if(length(unique(cnv.matrix[, label])) == 2){
    message("dichotomous outcome variable detected. Logistic regression is being done ...")
    model = "glm"
  }else if(is.numeric(cnv.matrix[, label])){
    message("continuous outcome variable detected. Linear regression is being done ...")
    model = "lm"
  }else if(is.factor(cnv.matrix[, label])){
    message("ordinal outcome variable detected. Ordinal regression is being done ...")
    model = "clm"
  }else{
    stop("Non dichotomous or continuous outcome variable detected. The burden test cannot be run")
  }
  
  if(permutation){
    ref.term <- sprintf("%s ~ %s", label, paste(covariates, collapse = " + "))
    message(sprintf("Permuting sample labels for %s times", nperm))
    if(BiasedUrn){
      lm.odds <- glm(ref.term, cnv.matrix, family = binomial (link = "logit"))
      d.odds <- exp(lm.odds$linear.predictors)
      
      n.case <- sum(cnv.matrix[, label])
      n.all <- length(cnv.matrix[, label])
      
      perm.hg <- BiasedUrn::rMFNCHypergeo(nran = nperm, m = rep(1, n.all), n = n.case, odds = d.odds)
    }else{
      perm.hg <- data.frame(1:nrow(cnv.matrix), nrow(cnv.matrix):1)
      for(i in 1:nperm){
        permuted <- sample(cnv.matrix[, label])
        perm.hg <- cbind(perm.hg, permuted)
      }
      
      perm.hg <- perm.hg[, -c(1:2)]
    }
    
  }
  
  perm.test.pvalues <- data.frame()
  test.out <- data.frame()
  
  for(this.gs in geneset){
    out.message <- sprintf("Testing %s", this.gs)
    
    message(out.message)
    
    feature <- this.gs
    
    this.covariates <- covariates
    
    ref.term <- sprintf("%s ~ %s", label, paste(this.covariates, collapse = " + "))
    add.term <- sprintf("%s + %s", ref.term, feature)
    coeff.term <- gsub("+ TotalRare", "", add.term)
    if(standardizeCoefficient & mean(cnv.matrix[, feature]) != 0 & sd(cnv.matrix[, feature]) != 0){
      cnv.matrix[, feature] <- scale(cnv.matrix[, feature])
    }
    
    if(model == "lm"){
      ref.model <- lm(ref.term, cnv.matrix)
      add.model <- lm(add.term, cnv.matrix)
      coeff.model <- lm(coeff.term, cnv.matrix)
      ano <- anova(ref.model, add.model, test = "Chisq")
    }else if(model == "glm"){
      ref.model <- glm(ref.term, cnv.matrix, family = binomial(link = "logit"))
      add.model <- glm(add.term, cnv.matrix, family = binomial(link = "logit"))
      coeff.model <- glm(coeff.term, cnv.matrix, family = binomial(link = "logit"))
      ano <- anova(ref.model, add.model, test = "Chisq")
    }else{
      ref.model <- ordinal::clm(ref.term, data = cnv.matrix)
      add.model <- ordinal::clm(add.term, data = cnv.matrix)
      coeff.model <- ordinal::clm(coeff.term, data = cnv.matrix)
      ano <- anova(ref.model, add.model)
    }
    
    names(ano)[length(names(ano))] <- "pvalue"
    pvalue <- ano$pvalue[2]
    coefficient <- coeff.model$coefficients[feature]
    conf <- tryCatch({confint.default(add.model)}, 
                     error = function(e){return(NA)})
    
    if(is.na(conf)){
      coeff.l <- 0
      coeff.u <- 0
    }else{
      coeff.l <- conf[feature, 1]
      coeff.u <- conf[feature, 2]
    }
    
    
    temp.out <- data.frame("geneset" = this.gs, "coefficient" = coefficient, 
                           "coeff.upper" = coeff.u, "coeff.lower" = coeff.l, "pvalue" = pvalue)
    test.out <- rbind(test.out, temp.out)
    
    if(permutation){
      for(iperm in 1:nperm){
        cnv.matrix$outcome.perm <- perm.hg[, iperm]
        ref.perm.term <- sprintf("outcome.perm ~ %s", paste(this.covariates, collapse = " + "))
        add.perm.term <- sprintf("%s + %s", ref.perm.term, feature)
        
        if(model == "lm"){
          ref.perm.model <- lm(ref.perm.term, cnv.matrix)
          add.perm.model <- lm(add.perm.term, cnv.matrix)
          ano.perm <- anova(ref.perm.model, add.perm.model, test = "Chisq")
        }else if(model == "glm"){
          ref.perm.model <- glm(ref.perm.term, cnv.matrix, family = binomial(link = "logit"))
          add.perm.model <- glm(add.perm.term, cnv.matrix, family = binomial(link = "logit"))
          ano.perm <- anova(ref.perm.model, add.perm.model, test = "Chisq")
        }else{
          ref.perm.model <- ordinal::clm(ref.perm.term, data = cnv.matrix)
          add.perm.model <- ordinal::clm(add.perm.term, data = cnv.matrix)
          ano.perm <- anova(ref.perm.model, add.perm.model)
        }
        
        names(ano.perm)[length(names(ano.perm))] <- "pvalue"
        coeff <- add.perm.model$coefficients[feature]
        perm.test.pvalues <- rbind(perm.test.pvalues, data.frame("feat" = this.gs, "pvalue" = ano.perm$pvalue[2], coeff))
      }
    }
  }
  
  if(permutation){
    message("Calculating permutation-based FDR")
    test.out$permFDR <- 1
    for(i in 1:nrow(test.out)){
      rec <- test.out[i, ]
      this.perm <- perm.test.pvalues
      
      actual <- sum(test.out$pvalue <= rec$pvalue)/nrow(test.out)
      perm <- sum(this.perm$pvalue <= rec$pvalue)/nrow(this.perm)
      
      fdr <- ifelse(perm/actual > 1, 1, perm/actual)
      
      test.out$permFDR[i] <- fdr
    }
  }
  
  test.out <- test.out[order(test.out$pvalue), ]
  for(i in 1:nrow(test.out)){
    test.out$permFDR[i] <- min(test.out$permFDR[i:nrow(test.out)])
  }
  
  list.out <- list(test.out, perm.test.pvalues)
  names(list.out) <- c("Test", "Permutation.Test")
  return(list.out)
}

getPermPvalue <- function(target.pvalue, actual.pvalues, perm.pvalues){
  actual.times <- sum(actual.pvalues <= target.pvalue, na.rm = T)
  perm.times <- sum(perm.pvalues <= target.pvalue, na.rm = T)
  permFDR <- (perm.times/length(perm.pvalues))/(actual.times/length(actual.pvalues))
  permFDR <- ifelse(permFDR > 1, 1, permFDR)
  return(permFDR)
}

getPermFDR <- function(i, dt){
  pvalue <- dt$pvalue[i]
  return(min(dt$permFDR[dt$pvalue >= pvalue]))
}