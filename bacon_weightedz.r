BACON_adjusted_WeightedZ <- function(Pvalue,Dir,SampleCase,SampleControl,MultipleCorrection=c("BH","holm","hochberg","bonferroni","BY"),seed=111){
if(nrow(Pvalue)<2 | nrow(Dir)<2){
   stop("Results for at least two studies needed")}
if(nrow(Pvalue)!=nrow(Dir)){
   stop("No. of studies in Pvalue data and Direction data do not match")}
if(ncol(Pvalue)!=ncol(Dir)){
   stop("No. of genes in Pvalue data and Direction data do not match")}
if(sum(rownames(Pvalue)!=rownames(Dir))>0){
   stop("Study names in Pvalue data and Direction data are not in same order")}
if(sum(colnames(Pvalue)!=colnames(Dir))>0){
   stop("Gene names in Pvalue data and Direction data are not in same order")}
if(nrow(Pvalue)!=length(SampleCase)){
   stop("No. of studies in Pvalue data and Sample sizes do not match")}
if(length(SampleCase)!=length(SampleControl)){
   stop("No. of sample sizes for two groups do not match")}
if(as.numeric(MultipleCorrection%in%c("BH","holm","hochberg","bonferroni","BY"))==0){
   stop("Multiple testing correction method must be 'BH', 'holm', 'hochberg', 'bonferroni', or 'BY'")}

else{
raw_pval <- Pvalue
raw_dir <- Dir
Ztrans_BacWeightedZ <- qnorm(raw_pval/2,lower.tail=FALSE)*raw_dir
Zvec_BacWeightedZ <- as.vector(Ztrans_BacWeightedZ)

library(bacon)
set.seed(seed)
BacWeightedZ <- bacon(Zvec_BacWeightedZ)
sig.hat_BacWeightedZ <- round(inflation(BacWeightedZ),2) 
mu.hat_BacWeightedZ <- round(bias(BacWeightedZ),2)
Mod_Ztrans_BacWeightedZ <- (Ztrans_BacWeightedZ-mu.hat_BacWeightedZ)/sig.hat_BacWeightedZ

N1 <- SampleCase
N2 <- SampleControl
Neff <- 4/((1/N1)+(1/N2))
w <- sqrt(Neff)	
BacWeightedZ_meta <- apply(Mod_Ztrans_BacWeightedZ,2,FUN=function(x) sum(x*w)/sqrt(sum(w^2)))
pval_BacWeightedZ <- 2*(1-pnorm(abs(-BacWeightedZ_meta)))
BacWeightedZ_adj_pval <- p.adjust(pval_BacWeightedZ,method=MultipleCorrection)
return(list(pvalue=data.frame(pval_BacWeightedZ, BacWeightedZ_adj_pval),
            EmpiricalEstimatedMean=mu.hat_BacWeightedZ,
            EmpiricalEstimatedSD=sig.hat_BacWeightedZ,
            BacWeightedZ_meta = BacWeightedZ_meta))
}
}