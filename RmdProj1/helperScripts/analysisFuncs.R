##==================================================
##  Function for chi-square, logistic regression
##==================================================

AssoAnalysis <- function (sig.list, snpFiles, sheets, dt, column, row2use)
{
    dataIN <- loadWorkbook(snpFiles)
    sheets.names <-  getSheets(dataIN)  
    for (sheet2parse in sheets)
    {
      sheetIN  <- readWorksheet(dataIN, sheet=sheet2parse,   startRow = 7)
      sig.chisq.case.I <- NULL
      sig.chisq.case.II <- NULL
      sig.chisq.case.III <- NULL
    
      for (i in 2:dim(sheetIN)[2])
      {
        snp.name <- colnames(sheetIN)[i]
        geno.pheno <- cbind(sheetIN[,i], dt[,column])
      
        ##  ONLY use RSV positive
        geno.pheno <- geno.pheno [row2use ,]
      
        snp.geno.1 <- geno.pheno [which(geno.pheno [,2] == 1),1]
        snp.geno.0 <- geno.pheno [which(geno.pheno [,2] == 0),1]
      
        b.1 <- c(length(which(snp.geno.1 == 0)) , length(which(snp.geno.1 == 1)), length(which(snp.geno.1 == 2)))
        b.0 <- c(length(which(snp.geno.0 == 0)) , length(which(snp.geno.0 == 1)), length(which(snp.geno.0 == 2)))   
      

        # case 1
        dm <- as.data.frame (rbind (b.1, b.0))
        colnames(dm) <- c("AA", "AB", "BB")
        rownames(dm) <- c ("case-1", "case-0")
        if ( !is.na (chisq.test(dm)$p.value) & chisq.test(dm)$p.value < 0.05)
        {
          sig.chisq <- list ("SNP" = snp.name, pval =  chisq.test(dm)$p.value)
          if (is.null(sig.chisq.case.I))
          {
            sig.chisq.case.I <- as.data.frame(sig.chisq)
          }else{
            sig.chisq.case.I  = rbind (sig.chisq.case.I ,as.data.frame(sig.chisq))
          }
        }
      
        # case 2 
        dm <- as.data.frame (rbind (b.1, b.0))
        dm <- matrix (c(b.1[1], (b.1[2] + b.1[3]), b.0[1], (b.0[2]+b.0[3])), nrow = 2, ncol =2, byrow = TRUE)
        colnames(dm) <- c("AA", "AB/BB")
        rownames(dm) <- c ("case-1", "case-0")
      
        if ( !is.na (chisq.test(dm)$p.value) & chisq.test(dm)$p.value < 0.05)
        {
          sig.chisq <- list ("SNP" = snp.name, pval =  chisq.test(dm)$p.value)
          if (is.null(sig.chisq.case.II))
          {
            sig.chisq.case.II <- as.data.frame(sig.chisq)
          }else{
            sig.chisq.case.II  = rbind (sig.chisq.case.II ,as.data.frame(sig.chisq))
          }
        }
      
        # case 3
        dm <- as.data.frame (rbind (b.1, b.0))
        dm <- matrix (c((b.1[1] + b.1[2]),  b.1[3], (b.0[1] + b.0[2]), b.0[3]), nrow = 2, ncol =2, byrow = TRUE)
        colnames(dm) <- c("AA/AB", "BB")
        rownames(dm) <- c ("case-1", "case-0")
      
        if ( !is.na (chisq.test(dm)$p.value) & chisq.test(dm)$p.value < 0.05)
        {
          sig.chisq <- list ("SNP" = snp.name, pval =  chisq.test(dm)$p.value)
          if (is.null(sig.chisq.case.III))
          {
            sig.chisq.case.III <- as.data.frame(sig.chisq)
          }else{
            sig.chisq.case.III  = rbind (sig.chisq.case.III ,as.data.frame(sig.chisq))
          }
        }
      
      }
      sig.association <- list ( "Case-I" = sig.chisq.case.I , "Case-II" = sig.chisq.case.II , "Case-III" =  sig.chisq.case.III) 
      sig.list[[sheet2parse]] <- sig.association 
  #    names(sig.list[[sheet2parse]] ) <- sheets.names[sheet2parse]
    }
    names(sig.list) <- sheets.names
    return (sig.list)
}

##  logistic regression: SAT on SNP-genotype
logisticAnalysis <- function (sig.list, snpFiles, sheets, dt, column, row2use)
{
  dataIN <- loadWorkbook(snpFiles)
  sheets.names <-  getSheets(dataIN)  
  for (sheet2parse in sheets)
  {
    sheetIN  <- readWorksheet(dataIN, sheet=sheet2parse,   startRow = 7)
    sig.chisq.case.I <- NULL
    sig.chisq.case.II <- NULL
    sig.chisq.case.III <- NULL
    
    for (i in 2:dim(sheetIN)[2])
    {
      snp.name <- colnames(sheetIN)[i]
      geno.pheno <- cbind(sheetIN[,i], dt[,column])
      
      ##  ONLY use RSV positive
      geno.pheno <- geno.pheno [row2use ,]
      
      ##  Exclude individual with missing SNP information
      if (length(which(is.na(geno.pheno[,1]))) > 0){
        geno.pheno <-  geno.pheno[-which(is.na(geno.pheno[,1])),]
      }

      # case 1   
      ##  |   AA  |   AB    |  BB   |
      ##  Keep all three phenotypes
      
      m <- as.data.frame(geno.pheno)
      m2 <- data.frame(sapply(m, function(x) as.numeric(as.character(x))))
      tmp.logistic= glm(m2[,2] ~ m2[,1], family=binomial("logit"))
      
      if ( !is.na(as.numeric(coef(summary(tmp.logistic))[, 4][2]))  & as.numeric(coef(summary(tmp.logistic))[, 4][2]) < 0.05 )
        
      {
        sig.chisq <- list ("SNP" = snp.name, pval =   as.numeric(coef(summary(tmp.logistic))[,4][2]))
        if (is.null(sig.chisq.case.I))
        {
          sig.chisq.case.I <- as.data.frame(sig.chisq)
        }else{
          sig.chisq.case.I  = rbind (sig.chisq.case.I ,as.data.frame(sig.chisq))
        }
      }
      
      # case 2
      ##     AA  |   AB/BB 
      ##    Convert "1" to "2"
      
      m <- as.data.frame(geno.pheno)
      m[which(m[,1] == 1),1] = 2 
      m2 <- data.frame(sapply(m, function(x) as.numeric(as.character(x))))
      tmp.logistic= glm(m2[,2] ~ m2[,1], family=binomial("logit"))
      if ( !is.na(as.numeric(coef(summary(tmp.logistic))[, 4][2]))  & as.numeric(coef(summary(tmp.logistic))[, 4][2]) < 0.05 )
      {
        sig.chisq <- list ("SNP" = snp.name, pval =   as.numeric(coef(summary(tmp.logistic))[,4][2]))
        if (is.null(sig.chisq.case.II))
        {
          sig.chisq.case.II <- as.data.frame(sig.chisq)
        }else{
          sig.chisq.case.II  = rbind (sig.chisq.case.II ,as.data.frame(sig.chisq))
        }
      }
      
      # case 3
      ##    AA/AB |   BB   |
      ##    Convert "1" to "0"
      
      m <- as.data.frame(geno.pheno)
      m[which(m[,1] == 1),1] = 0
      m2 <- data.frame(sapply(m, function(x) as.numeric(as.character(x))))
      tmp.logistic= glm(m2[,2] ~ m2[,1], family=binomial("logit"))
      if ( !is.na(as.numeric(coef(summary(tmp.logistic))[, 4][2]))  & as.numeric(coef(summary(tmp.logistic))[, 4][2]) < 0.05 )
      {
        sig.chisq <- list ("SNP" = snp.name, pval =   as.numeric(coef(summary(tmp.logistic))[,4][2]))
        if (is.null(sig.chisq.case.III))
        {
          sig.chisq.case.III <- as.data.frame(sig.chisq)
        }else{
          sig.chisq.case.III  = rbind (sig.chisq.case.III ,as.data.frame(sig.chisq))
        }
      }
    }
    sig.association <- list ( "Case-I" = sig.chisq.case.I , "Case-II" = sig.chisq.case.II , "Case-III" =  sig.chisq.case.III) 
    sig.list[[sheet2parse]] <- sig.association 
   # names(sig.list[[sheet2parse]] ) <- sheets.names[sheet2parse]
  }
  names(sig.list) <- sheets.names
  return (sig.list)
}


##  logistic regression: SAT on SNP-genotype with covariate information: corvCol
logisticAnalysisWCorrV <- function (sig.list, snpFiles, sheets, dt, column, row2use, corvCol)
{
  dataIN <- loadWorkbook(snpFiles)
  sheets.names <-  getSheets(dataIN)  
  for (sheet2parse in sheets)
  {
    sheetIN  <- readWorksheet(dataIN, sheet=sheet2parse,   startRow = 7)
    sig.chisq.case.I <- NULL
    sig.chisq.case.II <- NULL
    sig.chisq.case.III <- NULL
    
    for (i in 2:dim(sheetIN)[2])
    {
      
      snp.name <- colnames(sheetIN)[i]

      ## take genotype and covariate and construct datamatrix
      geno.pheno <- cbind(sheetIN[,i], dt[,c(column,corvCol)])     
      ##  ONLY use RSV positive
      geno.pheno <- geno.pheno [row2use ,]
      
      ##  Exclude individual with missing SNP information
      if (length(which(is.na(geno.pheno[,1]))) > 0){
        geno.pheno <-  geno.pheno[-which(is.na(geno.pheno[,1])),]
      }
      
      # case 1   
      ##  |   AA  |   AB    |  BB   |
      ##  Keep all three phenotypes
      
      m <- as.data.frame(geno.pheno)
      m2 <- data.frame(sapply(m, function(x) as.numeric(as.character(x))))
      tmp.logistic= glm(m2[,2] ~ m2[,1] * m2[,3], family=binomial("logit"))

      if ( !is.na(as.numeric(coef(summary(tmp.logistic))[, 4][4]))  & as.numeric(coef(summary(tmp.logistic))[, 4][4]) < 0.05 )
      {
        sig.chisq <- list ("SNP" = snp.name, pval =   as.numeric(coef(summary(tmp.logistic))[,4][4]))
        if (is.null(sig.chisq.case.I))
        {
          sig.chisq.case.I <- as.data.frame(sig.chisq)
        }else{
          sig.chisq.case.I  = rbind (sig.chisq.case.I ,as.data.frame(sig.chisq))
        }
      }
      
      # case 2
      ##     AA  |   AB/BB 
      ##    Convert "1" to "2"
      
      m <- as.data.frame(geno.pheno)
      m[which(m[,1] == 1),1] = 2 
      m2 <- data.frame(sapply(m, function(x) as.numeric(as.character(x))))
      tmp.logistic= glm(m2[,2] ~ m2[,1] * m2[,3], family=binomial("logit"))
      
      if ( !is.na(as.numeric(coef(summary(tmp.logistic))[, 4][4]))  & as.numeric(coef(summary(tmp.logistic))[, 4][4]) < 0.05 )
      {
        sig.chisq <- list ("SNP" = snp.name, pval =   as.numeric(coef(summary(tmp.logistic))[,4][4]))
        if (is.null(sig.chisq.case.II))
        {
          sig.chisq.case.II <- as.data.frame(sig.chisq)
        }else{
          sig.chisq.case.II  = rbind (sig.chisq.case.II ,as.data.frame(sig.chisq))
        }
      }
      
      # case 3
      ##    AA/AB |   BB   |
      ##    Convert "1" to "0"
      
      m <- as.data.frame(geno.pheno)
      m[which(m[,1] == 1),1] = 0
      m2 <- data.frame(sapply(m, function(x) as.numeric(as.character(x))))
      tmp.logistic= glm(m2[,2] ~ m2[,1] * m2[,3], family=binomial("logit"))
      
      if ( !is.na(as.numeric(coef(summary(tmp.logistic))[, 4][4]))  & as.numeric(coef(summary(tmp.logistic))[, 4][4]) < 0.05 )
      {
        sig.chisq <- list ("SNP" = snp.name, pval =   as.numeric(coef(summary(tmp.logistic))[,4][4]))
        if (is.null(sig.chisq.case.III))
        {
          sig.chisq.case.III <- as.data.frame(sig.chisq)
        }else{
          sig.chisq.case.III  = rbind (sig.chisq.case.III ,as.data.frame(sig.chisq))
        }
      }
    }
    sig.association <- list ( "Case-I" = sig.chisq.case.I , "Case-II" = sig.chisq.case.II , "Case-III" =  sig.chisq.case.III) 
    sig.list[[sheet2parse]] <- sig.association 
 #   names(sig.list[[sheet2parse]] ) <- sheets.names[sheet2parse]
  }
  names(sig.list) <- sheets.names
  return (sig.list)
}
