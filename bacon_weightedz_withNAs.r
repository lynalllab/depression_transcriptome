BACON_adjusted_WeightedZwithNA <- function(Pvalue,Dir,SampleCase,SampleControl,MultipleCorrection=c("BH","holm","hochberg","bonferroni","BY"),seed=111){
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
if(as.numeric(MultipleCorrection%in%c("BH","holm","hochberg","bonferroni","BY"))==0){
   stop("Multiple testing correction method must be 'BH', 'holm', 'hochberg', 'bonferroni', or 'BY'")}

else{
raw_pval <- Pvalue
raw_dir <- Dir
Ztrans_BacWeightedZ <- apply(raw_pval, 2, function(x) ifelse(is.na(x) == F, qnorm(x/2, lower.tail = F), NA))
Ztrans_BacWeightedZ <- Ztrans_BacWeightedZ*raw_dir
Zvec_BacWeightedZ <- as.vector(Ztrans_BacWeightedZ)

library(bacon)
set.seed(seed)
BacWeightedZ <- bacon(Zvec_BacWeightedZ, na.exclude = T)
sig.hat_BacWeightedZ <- round(inflation(BacWeightedZ),2) 
mu.hat_BacWeightedZ <- round(bias(BacWeightedZ),2)
Mod_Ztrans_BacWeightedZ <- (Ztrans_BacWeightedZ-mu.hat_BacWeightedZ)/sig.hat_BacWeightedZ

BacWeightedZ_meta <- c()
for(i in 1:ncol(raw_pval)){
N1 <- as.vector(na.omit(SampleCase[i,]))
N2 <- as.vector(na.omit(SampleControl[i,]))
Neff <- 4/((1/N1)+(1/N2))
w <- sqrt(Neff)	
BacWeightedZ_meta_temp <- apply(as.matrix(na.omit(Mod_Ztrans_BacWeightedZ[,i])),2,FUN=function(x) sum(x*w)/sqrt(sum(w^2)))
BacWeightedZ_meta<- c(BacWeightedZ_meta, BacWeightedZ_meta_temp)
}

pval_BacWeightedZ <- 2*(1-pnorm(abs(-BacWeightedZ_meta)))
BacWeightedZ_adj_pval <- p.adjust(pval_BacWeightedZ,method=MultipleCorrection)
return(list(pvalue=data.frame(pval_BacWeightedZ, BacWeightedZ_adj_pval),
            EmpiricalEstimatedMean=mu.hat_BacWeightedZ,
            EmpiricalEstimatedSD=sig.hat_BacWeightedZ,
            BacWeightedZ_meta = BacWeightedZ_meta))
}
}

