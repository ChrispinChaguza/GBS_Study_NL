##R code for the rank-based inverse normal transformation of quantitative phenotypes/traits in GWAS (slightly modified)
##https://yuxuanstat.com/posts/2020/06/rank-based-inverse-normal-transformation/

irnt <- function(pheno) {
  set.seed(1234)
  
  lm.fit = lm(pheno ~ 1) 
  r = residuals(lm.fit)
  
  numPhenos = length(which(!is.na(r)))
  quantilePheno = (rank(r, na.last="keep", ties.method="random")-0.5)/numPhenos
  phenoIRNT = qnorm(quantilePheno)	
  
  return(unname(phenoIRNT));
}
