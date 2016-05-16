library(plyr)
orgDetails <- data.frame(org=c("Sorghum","Maize","Soybean","Arabidopsis"),nChrs=c(10,10,20,5),nSNPs=c(126,640,93,943))

Maize = c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
Arabidopsis = c(30427671,19698289,23459830,18585056,26975502)
Sorghum = c(73727935,77694824,74408397,67966759,62243505,62192017,64263908,55354556,59454246,61085274)
Soybean = c(56831624,48577505,45779781,52389146,42234498,51416486,44630646,47837940,50189764,51566898,34766867,40091314,45874162,49042192,51756343,37887014,41641366,58018742,50746916,47904181)

chrLengths <- list(Maize,Arabidopsis,Sorghum,Soybean)
names(chrLengths) <- c("Maize","Arabidopsis","Sorghum","Soybean")
#rice = c(43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856)

# rMaize <- 1e6
# rSoy <- 1e6
# rCorn <- 1e6
# rArab <- 1e6
#   

nDatasets <- 100

#####Select nSNP chromsomes, from 
for(i in 1:nDatasets){
  dir.create(paste0("RandomDatasets/Random",i), showWarnings = FALSE)
  chrList <- alply(orgDetails,1,function(x){
      thisChr <- sort(sample(1:x$nChrs,x$nSNPs,replace=TRUE))
      print(x)
      data.frame(org=as.character(x$org),chr=thisChr,bp=unlist(lapply(1:x$nChrs,function(y) sample(1:get(as.character(x$org))[y],length(which(thisChr==y)),replace = FALSE))))
    })
  lapply(chrList, function(x){write.table(x,paste0("RandomDatasets/Random",i,"/Zn.Random.Permutation",i,".",unique(x$org),".csv"),sep=",",col.names=TRUE,row.names=FALSE)})
}
