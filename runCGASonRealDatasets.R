library(data.table)
library(ggplot2)
source("makePermuteDatasets.R")

#ranges <- data.frame(orgs=c("Arabidopsis","Maize","Sorghum","Soybean"),range=c(25000,1e6,1e5,1e6),stringsAsFactors = FALSE)
#ranges <- data.frame(orgs=c("Arabidopsis","Maize","Sorghum","Soybean"),range=c(25000,500000,1e5,500000),stringsAsFactors = FALSE)
ranges <- data.frame(orgs=c("Arabidopsis","Maize","Sorghum","Soybean"),range=c(25000,100000,1e5,100000),stringsAsFactors = FALSE)
sorg <- fread("./realDatasets/SAP.AllDatasetsSNPs.onlyCofactors.max.csv",stringsAsFactors = FALSE,sep=",",header=TRUE)
maize <- fread("./realDatasets/sigGWASsnpsCombinedIterations.allLoc.v3locs.csv",stringsAsFactors = FALSE,sep=",",header=TRUE)
arab <- fread("./realDatasets/AtIonome.csv",stringsAsFactors = FALSE,sep="\t",header=TRUE)
soy <- fread("./realDatasets/MLMM.maf0.05.allCof.SoyLM.onlyCofactors.csv",stringsAsFactors = FALSE,sep=",",header=TRUE)

arab$chrom <- as.numeric(gsub("Chr","",arab$chrom))
arab <- arab[arab$GWAS=="Seed",]

maize <- maize[maize$v3chr != "None",]

maize$org <- "Maize"
sorg$org <- "Sorghum"
arab$org <- "Arabidopsis"
soy$org <- "Soybean"

maize[,chr:=NULL,]
maize[,bp:=NULL,]

setnames(maize,old=c("org","v3chr","v3bp","el"),new=c("org","chr","bp","trait"))
setnames(arab,old=c("org","chrom","pos","Trait"),new=c("org","chr","bp","trait"))
soy$trait <- gsub("\\d+","",soy$trait)
maize$trait <- gsub("\\d+","",maize$trait)
arab <- arab[arab$trait %in% soy$trait]
sorg$trait <- gsub("\\d+","",sorg$phenotype)

traits <- Reduce(intersect,list(sorg$trait,maize$trait,arab$trait,soy$trait))
sorg <- sorg[sorg$trait %in% traits,]
maize <- maize[maize$trait %in% traits,]
arab <- arab[arab$trait %in% traits,]
soy <- soy[soy$trait %in% traits,]

allSnps <- rbind(maize[,c("org","chr","bp","trait"),with=FALSE],
              sorg[,c("org","chr","bp","trait"),with=FALSE],
              arab[,c("org","chr","bp","trait"),with=FALSE],
              soy[,c("org","chr","bp","trait"),with=FALSE])
allSnps[,chr:=as.integer(chr)]
####Make random datasets###
Maize = c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204)
Arabidopsis = c(30427671,19698289,23459830,18585056,26975502)
Sorghum = c(73727935,77694824,74408397,67966759,62243505,62192017,64263908,55354556,59454246,61085274)
Soybean = c(56831624,48577505,45779781,52389146,42234498,51416486,44630646,47837940,50189764,51566898,34766867,40091314,45874162,49042192,51756343,37887014,41641366,58018742,50746916,47904181)
chrLengths <- list(Maize,Arabidopsis,Sorghum,Soybean)
names(chrLengths) <- c("Maize","Arabidopsis","Sorghum","Soybean")

nPermute <- 100

for(i in unique(allSnps$org)){
  permSNPtable <- generateSNPs(chrLengths = chrLengths[[i]],nrow(allSnps[allSnps$org==i,])*100)
  for(j in traits){
    for(k in 1:nPermute){
      thisSampIndxs <- sample(nrow(permSNPtable),nrow(allSnps[allSnps$org==i & allSnps$trait==j,]),replace=FALSE)
      dir.create(paste0("RandomDatasets/Random",k), showWarnings = FALSE)
      write.table(permSNPtable[thisSampIndxs,],paste0("RandomDatasets/Random",k,"/",j,".Random.Permutation",k,".",i,".csv"),sep=",",col.names=TRUE,row.names=FALSE)
      permSNPtable <- permSNPtable[-thisSampIndxs,]
    }
  }
}

#for(i in traits){
#  orgDetails <- data.frame(org=c("Arabidopsis","Maize","Sorghum","Soybean"),nChrs=c(5,10,10,20),nSNPs=as.vector(table(allSnps$org[allSnps$trait==i])))
#  permuteDatasets(orgDetails,chrLengths,i,100)
#}

######Takes names of two organisms and a table with org,chr,start,end and returns indices from overlap table were both ranges match
getOverlapIndxs <- function(org1,org2,snpTable,orthoTable){
  setkey(snpTable,chr,start,end)
  #####Find indxs of genes that are in arabidopsis snp ranges####
  org1Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org1,],by.x=c("Gene1chr","Gene1start","Gene1end"),by.y=c("chr","start","end"),mult="first",type="within",which=TRUE))))
  #####Find indxs of genes that are in corn snp ranges####
  org2Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org2,],by.x=c("Gene2chr","Gene2start","Gene2end"),by.y=c("chr","start","end"),mult="first",type="within",which=TRUE))))
  intersect(org1Indxs,org2Indxs)
}



randomDatasets <- list.dirs("./RandomDatasets",recursive = FALSE,full.names = TRUE)

#compFiles <- list.files("./static/data/")
arabCorn <- fread("../OrthologTools/static/data/A.-thaliana-TAIR9-TAIR10_Z.-mays-AGPv3-5b+.csv",sep=",",stringsAsFactors = FALSE)
arabSorg <- fread("../OrthologTools/static/data/A.-thaliana-TAIR9-TAIR10_S.-bicolor-v2.0-v2.1.csv",sep=",",stringsAsFactors = FALSE)
soyCorn  <- fread("../OrthologTools/static/data/G.-max-v2.0-Wm82.a2.v1_Z.-mays-AGPv3-5b+.csv",sep=",",stringsAsFactors = FALSE)
soySorg  <- fread("../OrthologTools/static/data/G.-max-v2.0-Wm82.a2.v1_S.-bicolor-v2.0-v2.1.csv",sep=",",stringsAsFactors = FALSE)
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
dir.create("./randomDatasetResults/", showWarnings = FALSE)
dir.create("./realDatasetResults/", showWarnings = FALSE)
summaryTable <- data.frame()

orgs <- unique(allSnps$org)
for(i in orgs){
  allSnps$start[allSnps$org==i] <- allSnps$bp[allSnps$org==i]-ranges$range[ranges$orgs==i]
  allSnps$end[allSnps$org==i] <- allSnps$bp[allSnps$org==i]+ranges$range[ranges$orgs==i]
}

####Get the overlap genes for the four comparisons
for(thisTrait in traits){
  arabCornIndxs <- getOverlapIndxs("Arabidopsis","Maize",allSnps[allSnps$trait==thisTrait,],arabCorn)
  thisArabCornOverlaps <- arabCorn[arabCornIndxs,]
  arabSorgIndxs <- getOverlapIndxs("Arabidopsis","Sorghum",allSnps[allSnps$trait==thisTrait,],arabSorg)
  thisArabSorgOverlaps <- arabSorg[arabSorgIndxs,]
  
  soyCornIndxs <- getOverlapIndxs("Soybean","Maize",allSnps[allSnps$trait==thisTrait,],soyCorn)
  thisSoyCornOverlaps <- soyCorn[soyCornIndxs,]
  soySorgIndxs <- getOverlapIndxs("Soybean","Sorghum",allSnps[allSnps$trait==thisTrait,],soySorg)
  thisSoySorgOverlaps <- soySorg[soySorgIndxs,]
  
  mergeArab <- merge(thisArabCornOverlaps,thisArabSorgOverlaps,by="Gene1id",all=FALSE)
  mergeSoy <- merge(thisSoyCornOverlaps,thisSoySorgOverlaps,by="Gene1id",all=FALSE)
  fourway <- merge(mergeArab,mergeSoy,by=c("Gene2id.x","Gene2id.y"),all=FALSE)
  fourway[,c("Gene1org.y.x","Gene1Defline.y.x","Gene1Defline.x.y","Gene1end.y.x","Gene1start.y.x","Gene2Defline.x.y","Gene1chr.y.x",
             "Gene2org.x.y","Gene2end.x.y","Gene2start.x.y","Gene2end.y.y","Gene2start.y.y","Gene2Defline.y.y","Gene2chr.y.y","Gene2chr.x.y",
             "Gene1org.y.y","Gene2org.y.y","Gene1end.y.y","Gene1start.y.y","Gene1chr.y.y"):=NULL]
  setnames(fourway,c("MaizeID","SorghumID","ArabID","ArabOrg","MaizeOrg","RelationShipArabMaize","MaizeGeneEnd","MaizeGeneStart","ArabGeneEnd","ArabGeneStart",
                     "ArabGeneChr","MaizeGeneChr","ArabDefline","MaizeDefline","SorghumOrg","RelationShipArabSorghum","SorghumGeneEnd","SorghumGeneStart",
                     "SorghumGeneChr","SorghumDefline","SoyID","SoyOrg","RelationshipSoyMaize","SoyGeneEnd","SoyGeneStart","SoyGeneChr","RelationshipSoySorghum",
                     "SoyDefline"))
  
  write.table(fourway,file=paste0("./realDatasetResults/",thisTrait,".fourwayArabSoyCornSorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
  
  write.table(thisArabCornOverlaps,file=paste0("./realDatasetResults/",thisTrait,".ArabCornGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")  
  write.table(thisArabSorgOverlaps,file=paste0("./realDatasetResults/",thisTrait,".ArabSorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")  
  write.table(thisSoyCornOverlaps,file=paste0("./realDatasetResults/",thisTrait,".SoyCornGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
  write.table(thisSoySorgOverlaps,file=paste0("./realDatasetResults/",thisTrait,".SoySorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
  summaryTable <- rbind(summaryTable,data.frame(permutation="RealData",trait=thisTrait,Arabidopsis=length(unique(fourway$ArabID)),Soybean=length(unique(fourway$SoyID)),
                                                Sorghum=length(unique(fourway$SorghumID)),Maize=length(unique(fourway$MaizeID))))
}


for(thisDataset in randomDatasets){
  print(thisDataset)
  snpFiles <- list.files(thisDataset,full.names = TRUE)
  for(trait in traits){
    traitFiles <- grep(paste0("\\/",trait),snpFiles,value=TRUE)
    snps <- data.frame()
    for(i in traitFiles){
      thisTable <- read.table(i, header=TRUE, sep=",",stringsAsFactors = FALSE)
      thisTable$org <- strsplit(i,"\\.")[[1]][5]
      snps <- rbind(snps,thisTable)
    }
    orgs <- unique(snps$org)
    for(i in orgs){
      snps$start[snps$org==i] <- snps$bp[snps$org==i]-ranges$range[ranges$orgs==i]
      snps$end[snps$org==i] <- snps$bp[snps$org==i]+ranges$range[ranges$orgs==i]
    }
    setDT(snps)
    ####Get the overlap genes for the four comparisons
    arabCornIndxs <- getOverlapIndxs("Arabidopsis","Maize",snps,arabCorn)
    thisArabCornOverlaps <- arabCorn[arabCornIndxs,]
    arabSorgIndxs <- getOverlapIndxs("Arabidopsis","Sorghum",snps,arabSorg)
    thisArabSorgOverlaps <- arabSorg[arabSorgIndxs,]
    
    soyCornIndxs <- getOverlapIndxs("Soybean","Maize",snps,soyCorn)
    thisSoyCornOverlaps <- soyCorn[soyCornIndxs,]
    soySorgIndxs <- getOverlapIndxs("Soybean","Sorghum",snps,soySorg)
    thisSoySorgOverlaps <- soySorg[soySorgIndxs,]
    
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
    write.table(arabSorgCornSoy,file=paste0("./randomDatasetResults/",randomFolder,"/",trait,".fourwayArabSoyCornSorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
    
    write.table(thisArabCornOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/",trait,".ArabCornGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")  
    write.table(thisArabSorgOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/",trait,".ArabSorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")  
    write.table(thisSoyCornOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/",trait,".SoyCornGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(thisSoySorgOverlaps,file=paste0("./randomDatasetResults/",randomFolder,"/",trait,".SoySorghumGenes.csv"),row.names=FALSE,col.names=TRUE,sep=",")
    
    summaryTable <- rbind(summaryTable,data.frame(permutation=randomFolder,trait=trait,Arabidopsis=length(unique(arabSorgCornSoy$ArabID)),Soybean=length(unique(arabSorgCornSoy$SoyID)),
                                  Sorghum=length(unique(arabSorgCornSoy$SorghumID)),Maize=length(unique(arabSorgCornSoy$MaizeID))))
    rm(snps,soySorgCorn,arabSorgCorn,finalArabinAllSorg,finalArabinAllCorn,finalSoyinAllSorg,commonSorgIds,commonCornIds,
       commonSoySorg,commonSoyCorn,commonArabSorg,commonArabCorn,soySorgAlsoinCorn,commonArabIds,commonSoyIds,
       arabSorgAlsoinCorn,soySorgIndxs,soyCornIndxs,arabSorgIndxs,arabCornIndxs,thisSoySorgOverlaps,
       thisSoyCornOverlaps,thisArabSorgOverlaps,thisArabCornOverlaps,
       arabSorgCornSoy,arabSorgAlsoinCorn,finalSoyinAllCorn,finalSoyinAllSorg)
  }#end foreach trait
} #end for each dataset
write.table(summaryTable,"./randomDatasetResults/summaryTable.csv",sep=",",row.names=FALSE,col.names=TRUE)
write.table(ranges,"./randomDatasetResults/ranges.csv",sep=",",row.names=FALSE,col.names=TRUE)

meltedSummary <- melt(summaryTable)
meltedSummary$type <- "Bootstrap"
meltedSummary$type[meltedSummary$permutation=="RealData"] <- "Actual Dataset"
meltedSummary$color <- "Blue"
meltedSummary$color[meltedSummary$permutation=="RealData"] <- "Red"
meltedSummary$size <- 0.5
meltedSummary$size[meltedSummary$permutation=="RealData"] <- 3

head(meltedSummary)
pdf("permutationDensitiesvsActual.shortRanges.pdf")
for(trait in traits){
  #ggplot(meltedSummary[meltedSummary$trait==trait,],aes(y=value,x=variable,fill=type))+geom_boxplot()
  print(ggplot(meltedSummary[meltedSummary$trait==trait,],aes(value))+geom_density()+facet_wrap(~variable)+
    geom_rug(color=meltedSummary$color[meltedSummary$trait==trait],size=meltedSummary$size[meltedSummary$trait==trait])+
    ggtitle(paste0("Four-way overlap results for ",trait)))
}
dev.off()
  #geom_rug(meltedSummary[meltedSummary$trait==trait & meltedSummary$type=="Actual Dataset",],color="red",size="2")


#######FDR calc########
#proportion of permutation test exceeding actual
pptea <- data.frame()
for(trait in traits){
  thisTrait <- summaryTable[summaryTable$trait==trait,]
  propPermExceedsActualArab <- (length(which(thisTrait$Arabidopsis>=thisTrait$Arabidopsis[thisTrait$permutation=="RealData"]))-1)/100
  propPermExceedsActualSorg <- (length(which(thisTrait$Sorghum>=thisTrait$Sorghum[thisTrait$permutation=="RealData"]))-1)/100
  propPermExceedsActualSoy <- (length(which(thisTrait$Soybean>=thisTrait$Soybean[thisTrait$permutation=="RealData"]))-1)/100
  propPermExceedsActualCorn <- (length(which(thisTrait$Maize>=thisTrait$Maize[thisTrait$permutation=="RealData"]))-1)/100
  pptea <- rbind(pptea,data.frame(trait=trait,Arabidopsis=propPermExceedsActualArab,Sorghum=propPermExceedsActualSorg,Soybean=propPermExceedsActualSoy,
             Maize=propPermExceedsActualCorn))
  #thisTrait <- thisTrait[order(thisTrait$Arabidopsis,decreasing = TRUE),]
  #which(thisTrait$permutation=="RealData")[1]
}
write.table(pptea,"./realDatasetResults/ProportionofPermutationsExceedingActual.mediumRanges.csv",row.names=FALSE,col.names=TRUE,sep=",")
#####Use this code to spot check that there are snps in the range for genes found
#checkRow <- 2
#arabSorg[arabSorgIndxs[checkRow],]
#snps[(snps$org=="Arabidopsis") & (snps$chr==arabSorg[arabSorgIndxs[checkRow],Gene1chr]) & (snps$start <= arabSorg[arabSorgIndxs[checkRow],Gene1start]) & (snps$end >= arabSorg[arabSorgIndxs[checkRow],Gene1end]),]
#snps[(snps$org=="Sorghum") & (snps$chr==arabSorg[arabSorgIndxs[checkRow],Gene2chr]) & (snps$start <= arabSorg[arabSorgIndxs[checkRow],Gene2start]) & (snps$end >= arabSorg[arabSorgIndxs[checkRow],Gene2end]),]
#snps[(snps$org=="Soybean") & (snps$chr==7) & (snps$start <= 12550508) & (snps$end >= 12555061)]
