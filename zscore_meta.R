zscore_meta <- function(dat="",size=""){
  
  for(i in 1:length(dat)){
    #read data
    sum_stat <- read.table(dat[i],header=T)
    #quality filtering can be done here e.g exclude snps with maf of zero
    #select snp,beta and p-value
    z_dat <- sum_stat
    z_dat <- z_dat[,c("SNP","b","p")]
    names(z_dat)[-1] <- paste0(names(z_dat)[-1],i)
    
    if(i==1){
      sumdat <- z_dat
    }else{
      sumdat <- merge(sumdat,z_dat,by="SNP")
      }
  }
  
  ##sample size based meta  analysis adapted from Goutam Sahana Summer school lecture note
  row.names(sumdat) <- sumdat$SNP
  b <- sumdat[,paste0("b",1:length(dat))]
  P <- sumdat[,paste0("p",1:length(dat))]
  
  delta <- b/abs(b) ## direction of effect for study i
  
  #function to convert p values to Zscore
  p2zscore <- function(x){
    qnorm(1-(x/2))
  }
  
  Zi <-apply(P,c(1,2),p2zscore)
  Zi <- Zi*delta
  Wi <- sqrt(size)
  Zw <- t(t(Zi)*Wi)
  Z <- apply(Zw,1,sum)/sqrt(sum(Wi^2))
  Z <- as.data.frame(Z)
  Z$Pmeta <- 2*pnorm(-abs(Z$Z))
  Z$N <- sum(size)
  
  Zmeta <- merge(sum_stat[,c("Chr","SNP","bp")],Z,by.x="SNP",by.y = 0)
  
  return(Zmeta)
}


