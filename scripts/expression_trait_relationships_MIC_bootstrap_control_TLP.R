
# define what to run

# treatment
treatment <- 'control1'
#treatment <- 'drought1'

# physiology measure
#phys_measure <- 'leaf_WP_avg'
phys_measure <- 'TLP_avg'

# setup
library(minerva)
library(parallel) # to run on multiple cores

# load data
# physiology data - named 'phys'
load('phys_data_cleaned_dataframe.Rdata')

# expression data
expr <- read.csv('gene_expression_normalized_voom_batch_effects_removed.csv', row.names = 1, check.names = F)

# setup functions

# function to get confidence interval
ci.get <- function(distribution, confidence.level = 0.95, ...){
  ci.range <- c(0, confidence.level)+(1-confidence.level)/2
  return(quantile(distribution, probs = ci.range, names = F, ...))
}

# p value calculator
# not: have not tested two-sided option extensively, but it 
# should be fine for doing one-sided test whether a value (MIC)
# is greater than expected from a null distribution of MIC
pcalc <- function(nullDist, observedStatistic, sided = "two"){
  switch (sided,two = {
    # two-sided test: actual value is different than null 
    # (either bigger or smaller)
    # calculate for most likely scenario, then multiply by 2
    # to account for using both tails
      pValue <- (min(sum(nullDist > observedStatistic), sum(nullDist < observedStatistic))/length(nullDist))*2
    },
    one = {
      # one-sided test: actual is bigger or smaller than null
      # get null observations that are less than actual
      # return pvalue for probability that observed is 
      # greater than null
      # (would have to change for prob observed is less than null)
      pValue <- sum(nullDist > observedStatistic)/length(nullDist)
    },
    
    {
      stop("Unknown sidedness")
    }
  )
  return(pValue)
}



# bootstrap function
mic.boot <- function(xvar, yvar, nBootstraps = 10000, ci.level = 0.95, alpha.level = 0.05, sided = "one", ncores = 1, ...){
  
  # get null distribution of MIC values
  MIC.null <- vector()
  for(b in 1:nBootstraps){
    nsamples <- length(xvar)
    resampledX <- xvar[sample(1:nsamples, replace = F)]
    mic <- mine(resampledX, yvar, n.cores = ncores)
    MIC.null[b] <- mic$MIC
  }
  # get actual MIC
  MIC.actual <- mine(xvar, yvar)$MIC
  
  return(
    list(
      MIC.actual = MIC.actual,
      null.confidence.interval = ci.get(MIC.null, ci.level),
      p.value = pcalc(nullDist = MIC.null, observedStatistic = MIC.actual, sided = sided)
    )
  )
  
}

#################################################
# run bootstrap

# get subset of expression data
samples <- phys$Tree[which(phys$Treatment == treatment)]
expr <- expr[,which(colnames(expr) %in% samples)]
phys <- phys[which(phys$Tree %in% colnames(expr)),] # use expression samples because not all in treatment were sequenced

# remove 069-15 outlier
expr <- expr[, which(colnames(expr) != '069-15')]
phys <- phys[which(phys$Tree != '069-15'),]

# get measurement vector
column <- which(colnames(phys) == phys_measure)
measure <- phys[,column]

# run bootstrap for all genes
# this is calculating pvals for MIC of expression and TLP for control individuals
# haven't been able to get it to run with multiple cores
top <- nrow(expr)
mic.boot.result <- as.data.frame(matrix(nrow = top, ncol = 2))
colnames(mic.boot.result) <- c('MIC', 'p.value')
rownames(mic.boot.result) <- rownames(expr)[1:top]

for(n in 1:top){
  mic <- mic.boot(measure, as.numeric(expr[n,]), nBootstraps = 5000, ncores = 1)
  mic.boot.result[n,1] <- mic$MIC.actual
  mic.boot.result[n,2] <- mic$p.value
  if (n %% 100 == 0){
    print(c('done with gene', n), quote = F)
  }
}

save(mic.boot.result, file = paste('MIC_bootstrap_results', treatment, colnames(phys)[column], '.Rdata', sep = '_'))
