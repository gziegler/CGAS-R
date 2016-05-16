library(data.table)

ranges <- data.frame(orgs=c("Arabidopsis","Maize","Sorghum","Soybean"),range=c(25000,1e6,1e5,1e6),stringsAsFactors = FALSE)
randomDatasets <- list.dirs("./RandomDatasets",recursive = FALSE,full.names = TRUE)

#compFiles <- list.files("./static/data/")
arabCorn <- fread("./static/data/A.-thaliana-TAIR9-TAIR10_Z.-mays-AGPv3-5b+.csv",sep=",",stringsAsFactors = FALSE)
arabSorg <- fread("./static/data/A.-thaliana-TAIR9-TAIR10_S.-bicolor-v2.0-v2.1.csv",sep=",",stringsAsFactors = FALSE)
soyCorn  <- fread("./static/data/G.-max-v2.0-Wm82.a2.v1_Z.-mays-AGPv3-5b+.csv",sep=",",stringsAsFactors = FALSE)
soySorg  <- fread("./static/data/G.-max-v2.0-Wm82.a2.v1_S.-bicolor-v2.0-v2.1.csv",sep=",",stringsAsFactors = FALSE)
arabCorn[,V14:=NULL,]
arabSorg[,V14:=NULL,]
soyCorn[,V14:=NULL,]
soySorg[,V14:=NULL,]
setnames(arabCorn,c("Gene1id","Gene1org","Gene2id","Gene2org","Relationship","Gene2end","Gene2start","Gene1end","Gene1start","Gene1chr","Gene2chr","Gene1Defline","Gene2Defline"))
setnames(arabSorg,c("Gene1id","Gene1org","Gene2id","Gene2org","Relationship","Gene2end","Gene2start","Gene1end","Gene1start","Gene1chr","Gene2chr","Gene1Defline","Gene2Defline"))
setnames(soyCorn,c("Gene1id","Gene1org","Gene2id","Gene2org","Relationship","Gene2end","Gene2start","Gene1end","Gene1start","Gene1chr","Gene2chr","Gene1Defline","Gene2Defline"))
setnames(soySorg,c("Gene1id","Gene1org","Gene2id","Gene2org","Relationship","Gene2end","Gene2start","Gene1end","Gene1start","Gene1chr","Gene2chr","Gene1Defline","Gene2Defline"))
arabCorn[,Gene1chr := gsub("Chr","",Gene1chr),]
arabSorg[,Gene1chr := gsub("Chr","",Gene1chr),]
arabSorg[,Gene2chr := gsub("Chr","",Gene2chr),]
soyCorn[,Gene1chr := gsub("Chr","",Gene1chr),]
soySorg[,Gene1chr := gsub("Chr","",Gene1chr),]
soySorg[,Gene2chr := gsub("Chr","",Gene2chr),]

#make chromosome and base pair columns numeric (this will convert some non-numeric chromsome IDs to NA (with warnings), they aren't the major chromosomes anyways)
for(col in c("Gene1chr","Gene2chr","Gene1start","Gene2start","Gene1end","Gene2end")) set(arabCorn, j=col, value=as.numeric(arabCorn[[col]]))
for(col in c("Gene1chr","Gene2chr","Gene1start","Gene2start","Gene1end","Gene2end")) set(arabSorg, j=col, value=as.numeric(arabSorg[[col]]))
for(col in c("Gene1chr","Gene2chr","Gene1start","Gene2start","Gene1end","Gene2end")) set(soyCorn, j=col, value=as.numeric(soyCorn[[col]]))
for(col in c("Gene1chr","Gene2chr","Gene1start","Gene2start","Gene1end","Gene2end")) set(soySorg, j=col, value=as.numeric(soySorg[[col]]))

for(thisDataset in randomDatasets){
  snpFiles <- list.files(thisDataset,full.names = TRUE)
  snps <- data.frame()
  for(i in snpFiles){
    snps <- rbind(snps,read.table(i, header=TRUE, sep=",",stringsAsFactors = FALSE))
  }
  orgs <- unique(snps$org)
  for(i in orgs){
    snps$start[snps$org==i] <- snps$bp[snps$org==i]-ranges$range[ranges$orgs==i]
    snps$end[snps$org==i] <- snps$bp[snps$org==i]+ranges$range[ranges$orgs==i]
  }
  setDT(snps)
  setkey(snps,chr,start,end)
  
  getOverlapIndxs <- function(org1,org2,snpTable,orthoTable){
    #####Find indxs of genes that are in arabidopsis snp ranges####
    org1Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org1,],by.x=c("Gene1chr","Gene1start","Gene1end"),by.y=c("chr","start","end"),mult="first",type="within",which=TRUE))))
    #####Find indxs of genes that are in corn snp ranges####
    org2Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org2,],by.x=c("Gene2chr","Gene2start","Gene2end"),by.y=c("chr","start","end"),mult="first",type="within",which=TRUE))))
    intersect(org1Indxs,org2Indxs)
  }
  
  ####Get the overlap genes for the four comparisons
  arabCornIndxs <- getOverlapIndxs("Arabidopsis","Maize",snps,arabCorn)
  thisArabCornOverlaps <- arabCorn[arabCornIndxs,]
  arabSorgIndxs <- getOverlapIndxs("Arabidopsis","Sorghum",snps,arabSorg)
  thisArabSorgOverlaps <- arabSorg[arabSorgIndxs,]
  
  soyCornIndxs <- getOverlapIndxs("Soybean","Maize",snps,soyCorn)
  thisSoyCornOverlaps <- soyCorn[soyCornIndxs,]
  soySorgIndxs <- getOverlapIndxs("Soybean","Sorghum",snps,soySorg)
  thisSoySorgOverlaps <- soySorg[soySorgIndxs,]
  
  #####Use this code to spot check that there are snps in the range for genes found
    #checkRow <- 2
    #arabSorg[arabSorgIndxs[checkRow],]
    #snps[(snps$org=="Arabidopsis") & (snps$chr==arabSorg[arabSorgIndxs[checkRow],Gene1chr]) & (snps$start <= arabSorg[arabSorgIndxs[checkRow],Gene1start]) & (snps$end >= arabSorg[arabSorgIndxs[checkRow],Gene1end]),]
    #snps[(snps$org=="Sorghum") & (snps$chr==arabSorg[arabSorgIndxs[checkRow],Gene2chr]) & (snps$start <= arabSorg[arabSorgIndxs[checkRow],Gene2start]) & (snps$end >= arabSorg[arabSorgIndxs[checkRow],Gene2end]),]
    #snps[(snps$org=="Soybean") & (snps$chr==7) & (snps$start <= 12550508) & (snps$end >= 12555061)]
  
  #####Find which arabidopsis IDs are in common between the comparisons#####
  #####Find which soybean IDs are in common between the comparisons#####
  commonArabIds <- intersect(thisArabCornOverlaps$Gene1id,thisArabSorgOverlaps$Gene1id)
  commonSoyIds  <- intersect(thisSoyCornOverlaps$Gene1id,thisSoySorgOverlaps$Gene1id)
  
  #####Narrow down the ortholog lists to only those genes that are in common between the two sets#####
  commonArabCorn <- thisArabCornOverlaps[Gene1id %in% commonArabIds,]
  commonArabSorg <- thisArabSorgOverlaps[Gene1id %in% commonArabIds,]
  commonSoyCorn  <- thisSoyCornOverlaps[Gene1id %in% commonSoyIds,]
  commonSoySorg  <- thisSoySorgOverlaps[Gene1id %in% commonSoyIds,]
  
  
  ###Find which corn IDs are in common between the soy-corn comp and the arab-corn comp###
  commonCornIds  <- intersect(commonArabCorn$Gene2id,commonSoyCorn$Gene2id)
  
  ###Subset arab-sorg and soy-sorg, based on arabIDs and soyIDs present in rows containign corn IDs###
  arabSorgAlsoinCorn <- commonArabSorg[Gene1id %in% commonArabCorn$Gene1id[commonArabCorn$Gene2id %in% commonCornIds]]
  soySorgAlsoinCorn  <- commonSoySorg[Gene1id %in% commonSoyCorn$Gene1id[commonSoyCorn$Gene2id %in% commonCornIds]]
  
  commonSorgIds  <- intersect(arabSorgAlsoinCorn$Gene2id,soySorgAlsoinCorn$Gene2id)
  
  finalArabinAllSorg <- arabSorgAlsoinCorn[Gene2id %in% commonSorgIds]
  finalArabinAllCorn <- commonArabCorn[Gene1id %in% finalArabinAllSorg$Gene1id]
  finalSoyinAllSorg  <- soySorgAlsoinCorn[Gene2id %in% commonSorgIds]
  finalSoyinAllCorn <- commonSoyCorn[Gene1id %in% finalSoyinAllSorg$Gene1id]
  
  setnames(finalArabinAllSorg,c("ArabID","ArabOrg","SorghumID","SorghumOrg","RelationShipArabSorghum","SorghumGeneEnd","SorghumGeneStart","ArabGeneEnd","ArabGeneStart","ArabGeneChr","SorghumGeneChr","ArabDefline","SorghumDefline"))
  setnames(finalArabinAllCorn,c("ArabID","ArabOrg","MaizeID","MaizeOrg","RelationShipArabMaize","MaizeGeneEnd","MaizeGeneStart","ArabGeneEnd","ArabGeneStart","ArabGeneChr","MaizeGeneChr","ArabDefline","MaizeDefline"))

  setnames(finalSoyinAllCorn,c("SoyID","SoyOrg","MaizeID","MaizeOrg","RelationshipSoyMaize","MaizeGeneEnd","MaizeGeneStart","SoyGeneEnd","SoyGeneStart","SoyGeneChr","MaizeGeneChr","SoyDefline","MaizeDefline"))
  setnames(finalSoyinAllSorg,c("SoyID","SoyOrg","SorghumID","SorghumOrg","RelationshipSoySorghum","SorghumGeneEnd","SorghumGeneStart","SoyGeneEnd","SoyGeneStart","SoyGeneChr","SorghumGeneChr","SoyDefline","SorghumDefline"))
  
  arabSorgCorn <- merge(finalArabinAllSorg,finalArabinAllCorn[,.(ArabID,MaizeID,MaizeOrg,MaizeGeneEnd,MaizeGeneStart,MaizeDefline)],by="ArabID",all=TRUE,allow.cartesian=TRUE)
  soySorgCorn  <- merge(finalSoyinAllSorg,finalSoyinAllCorn[,.(SoyID,MaizeID,MaizeOrg,MaizeGeneEnd,MaizeGeneStart,MaizeDefline)],by="SoyID",all=TRUE,allow.cartesian=TRUE)
  
#  setkey(arabSorgCorn,SorghumID,MaizeID)
#  setkey(soySorgCorn,SorghumID,MaizeID)
#  merged <- arabSorgCorn[soySorgCorn,allow.cartesian=TRUE]
  
  arabSorgCornSoy <- merge(arabSorgCorn,soySorgCorn[,.(SoyID,SoyOrg,SoyGeneStart,SoyGeneEnd,SoyGeneChr,SoyDefline,SorghumID,MaizeID)],by=c("SorghumID","MaizeID"),all=FALSE,allow.cartesian = TRUE)

  ###########Write it all out#######
  randomFolder <- strsplit(thisDataset,"/")[[1]][3]
  dir.create(paste0("./randomDatasetResults/",randomFolder), showWarnings = FALSE)
  write.table(arabSorgCornSoy,file=paste0("./randomDatasetResults/",randomFolder,"/fourwayArabSoyCornSorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")

  write.table(thisArabCornOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/ArabCornGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")  
  write.table(thisArabSorgOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/ArabSorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")  
  write.table(thisSoyCornOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/SoyCornGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
  write.table(thisSoySorgOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/SoySorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
  
  rm(snps,soySorgCorn,arabSorgCorn,finalArabinAllSorg,finalArabinAllCorn,finalSoyinAllCorn,finalSoyinAllSorg,commonSorgIds,commonCornIds,commonSoySorg,commonSoyCorn,commonArabSorg,commonArabCorn,soySorgAlsoinCorn,commonArabIds,commonSoyIds,
     arabSorgAlsoinCorn,soySorgIndxs,soyCornIndxs,arabSorgIndxs,arabCornIndxs,thisSoySorgOverlaps,thisSoyCornOverlaps,thisArabSorgOverlaps,thisArabCornOverlaps,
     arabSorgCornSoy,arabSorgAlsoinCorn,finalSoyinAllCorn,finalSoyinAllSorg)  
  
} #end for each dataset


