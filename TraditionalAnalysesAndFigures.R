#TeamEpidermis - Figure Generation Scripts
#load libraries
library(Momocs)
#library(plyr)
library(ggplot2)

#read in RV data/organize in table
dataDir = "G:/My\ Drive/TeamEpidermis/DataAnalysis/Rozi_CellsAll"
list.files(path = dataDir, full = T, pattern = "txt") %>% import_txt() %>% Out()-> l
measures <- matrix(data = NA, nrow = length(l), ncol = 4)
for(i in 1:length(l))
{
	measures[i,] <- c(coo_area(l[i]), coo_eccentricityboundingbox(l[i]), coo_solidity(l[i]), coo_circularity(l[i]) )
}

RVmeasuresDF <- data.frame(measures)
colnames(RVmeasuresDF) <- c('Area', 'AspectRatio', 'Solidity', 'Circularity')
RVmeasuresDF$Sample_Number <- names(l)
RVmeasuresDF <- RVmeasuresDF[,c(5,3,2,1,4)]

#add factors for calculations
rv.fac <- read.csv("G:\\My Drive\\TeamEpidermis\\Manuscript\\Revision\\rv_fac.txt")
RVmeasuresDF <- merge(RVmeasuresDF,rv.fac[,c(2:6)],by.x="Sample_Number",by.y="sample_number")

#read in GP data/organize table
gp_list <- list.files(path = "G:\\My Drive\\TeamEpidermis\\DataAnalysis\\GraceCells", full.names = TRUE, pattern = "jpg")
gp_cells <- import_jpg(gp_list)
gp_cells <- Out(gp_cells)

measures <- matrix(data = NA, nrow = length(gp_cells), ncol = 4)
for(i in 1:length(gp_cells))
{
	measures[i,] <- c(coo_area(gp_cells[i]), coo_eccentricityboundingbox(gp_cells[i]), coo_solidity(gp_cells[i]), coo_circularity(gp_cells[i]) )
}
GPmeasuresDF <- data.frame(measures)
colnames(GPmeasuresDF) <- c('Area', 'AspectRatio', 'Solidity', 'Circularity')
GPmeasuresDF$Sample_Number <- names(gp_cells)
GPmeasuresDF <- GPmeasuresDF[,c(5,3,2,1,4)]

gp_fac <- read.csv("gp_fac.txt")
GPmeasuresDF <- merge(GPmeasuresDF,gp_fac[,c(2:5)],by.x="Sample_Number",by.y="sample_number")

#cut RVmeasuresDF$sample_ID, not important and incongruous with GPmeasuresDF
RVmeasuresDF <- RVmeasuresDF[,-6]

#bind two data frames
AllmeasuresDF <- rbind(RVmeasuresDF,GPmeasuresDF)
#set clades so that ferns match
AllmeasuresDF$clade[AllmeasuresDF$clade == "fern"] = "ferns"
colnames(AllmeasuresDF)[6:8] <- c("Species_Code","Side","Clade")
write.table(AllmeasuresDF, file = 'all.measures.csv', append = F, quote = F, sep = ",", col.names = T, row.names = F)

#calculate means for each species/leaf side
Species <- unique(AllmeasuresDF$Species_Code)
meanTable <- matrix(data = NA, nrow = length(Species), ncol = 13)
for(i in 1:length(Species)){
	tempDF <- AllmeasuresDF[AllmeasuresDF$Species_Code %in% Species[i],]
	#I know this isn't the best way to do it, but I think it is more clear what is happening
	solAllMean <- mean(tempDF$Solidity)
	solAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),2])
	solAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),2])
	solMeanDiff <- solAdMean - solAbMean 
	arAllMean <- mean(tempDF$AspectRatio)
	arAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),3])
	arAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),3])
	areaAllMean <- mean(tempDF$Area)
	areaAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),4])
	areaAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),4])
	cirAllMean <- mean(tempDF$Circularity)
	cirAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),5])
	cirAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),5])
	
	meanTable[i,] = c(solAllMean, solAdMean, solAbMean, solMeanDiff, arAllMean, arAdMean, arAbMean,
				areaAllMean, areaAdMean, areaAbMean, cirAllMean, cirAdMean, cirAbMean)
}

meanTableDF <- data.frame(meanTable)
meanTableDF$Species_Code <- Species
colnames(meanTableDF) <- c("Sol_Mean","Sol_Ad_Mean","Sol_Ab_Mean","Sol_Mean_Diff","AR_Mean","AR_Ad_Mean","AR_Ab_Mean",
				   "Area_Mean","Area_Ad_Mean","Area_Ab_Mean","Cir_Mean","Cir_Ad_Mean","Cir_Ab_Mean","Species_Code")
speciesKey <- AllmeasuresDF[,c("Species_Code","Clade")]
head(meanTableDF <- merge(meanTableDF,unique(speciesKey),by="Species_Code",all.x=TRUE,all.y=FALSE))
write.table(meanTableDF, file = 'allData.means.csv', append = F, quote = F, sep = ",", col.names = T, row.names = F)

#some code to figure out which sample has the median value
whichmedian <- function(x) which.min(abs(x - median(x)))

#set NAs to unknown for side
levels(gp_fac$side) <-  c("abaxial","adaxial","unknown")
gp_fac$side[is.na(gp_fac$side)] = "unknown"

#print all median cell shapes as representatives
gp_fac.species <- unique(gp_fac$sample)
for(i in 1:length(gp_fac.species)){
	tempCoo <- slice(gp_cells,grep(gp_fac.species[i],names(gp_cells)))
	temp_fac <- gp_fac[gp_fac$sample %in% gp_fac.species[i],]
	for(j in unique(temp_fac$side)){
		tempCoo2 <- slice(tempCoo,temp_fac$side == j)
		tempMeas <- vector()
		for(k in 1:length(tempCoo2)){
			tempMeas[k] <- coo_solidity(tempCoo2[k])
		}
		med <- whichmedian(tempMeas)
		pdf(file=paste0(gp_fac.species[i],"_",j,".medianShape.pdf"))
		coo_plot(tempCoo2[med],col="blue")
		dev.off()
	}
}


#make violin plot for all AR data
#add an extra column to make plotting easy
AllmeasuresDF$Plot <- "Yes"
q <- ggplot(AllmeasuresDF[AllmeasuresDF$AspectRatio <= 1,],aes(x=Plot,y=AspectRatio)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = FALSE) + theme_minimal()

r <- ggplot(AllmeasuresDF[AllmeasuresDF$AspectRatio <= 1,],aes(x=Clade,y=AspectRatio, colour=Clade, fill=Clade)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = FALSE,alpha=0.5) + 
scale_x_discrete(limits=c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots")) +
labs(x = "", y = "Aspect Ratio") + ylim(0,1) + theme_minimal()
r

#do area with scales
scale <- read.table("CellMagnifications_GraceAndJeff.txt",header=F)
colnames(scale) = c("sample","Magn")
scale$Scale <- 2700
scale$Scale[scale$Magn == 40] = 108
scale$Scale[scale$Magn == 100] = 270
scale$Scale[scale$Magn == 200] = 540


GPscaledMeasuresDF <- merge(GPmeasuresDF,scale)
head(GPscaledMeasuresDF)
GPscaledMeasuresDF$Scaled_Area <- GPscaledMeasuresDF$Area*(100/GPscaledMeasuresDF$Scale)^2

RVscaledMeasuresDF <- merge(RVmeasuresDF,scale,all.x=TRUE)
RVscaledMeasuresDF$Scaled_Area <- RVscaledMeasuresDF$Area

#remake measures table with scaled area
AllScaledMeasuresDF <- rbind(RVscaledMeasuresDF,GPscaledMeasuresDF)

names(AllScaledMeasuresDF)[1] <- "Species_Code"
names(AllScaledMeasuresDF)[7] <- "Side"
names(AllScaledMeasuresDF)[8] <- "Clade"

write.table(AllScaledMeasuresDF, file = 'all.measures.scaled.csv', append = F, quote = F, sep = ",", col.names = T, row.names = F)
write.table(allData.measures.scaled, file = 'all.measures.scaled.csv', append = F, quote = F, sep = ",", col.names = T, row.names = F)


#calculate means for each species/leaf side
Species <- unique(AllScaledMeasuresDF$Species_Code)
meanTable <- matrix(data = NA, nrow = length(Species), ncol = 13)
for(i in 1:length(Species)){
	tempDF <- AllScaledMeasuresDF[AllScaledMeasuresDF$Species_Code %in% Species[i],]
	#I know this isn't the best way to do it, but I think it is more clear what is happening
	solAllMean <- mean(tempDF$Solidity)
	solAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),3])
	solAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),3])
	solMeanDiff <- solAdMean - solAbMean 
	arAllMean <- mean(tempDF$AspectRatio)
	arAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),4])
	arAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),4])
	areaAllMean <- mean(tempDF$Scaled_Area)
	areaAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),11])
	areaAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),11])
	cirAllMean <- mean(tempDF$Circularity)
	cirAdMean <- mean(tempDF[grep("adaxial",tempDF$Side),6])
	cirAbMean <- mean(tempDF[grep("abaxial",tempDF$Side),6])
	
	meanTable[i,] = c(solAllMean, solAdMean, solAbMean, solMeanDiff, arAllMean, arAdMean, arAbMean,
				areaAllMean, areaAdMean, areaAbMean, cirAllMean, cirAdMean, cirAbMean)
}

meanScaledTableDF <- data.frame(meanTable)
meanScaledTableDF$Species_Code <- Species
colnames(meanScaledTableDF) <- c("Sol_Mean","Sol_Ad_Mean","Sol_Ab_Mean","Sol_Mean_Diff","AR_Mean","AR_Ad_Mean","AR_Ab_Mean",
				   "Scaled_Area_Mean","Scaled_Area_Ad_Mean","Scaled_Area_Ab_Mean","Cir_Mean","Cir_Ad_Mean","Cir_Ab_Mean","Species_Code")
speciesKey <- AllScaledMeasuresDF[,c("Species_Code","Clade")]
head(meanScaledTableDF <- merge(meanScaledTableDF,unique(speciesKey),by="Species_Code",all.x=TRUE,all.y=FALSE))
write.table(meanScaledTableDF, file = 'allData.means.scaled.csv', append = F, quote = F, sep = ",", col.names = T, row.names = F)
#Need to read in allData.means.scaled.csv after adding 9.01-9.2
allData.means.scaled <- read.csv("allData.means.scaled.csv",header=T)
allData.means.scaled$Clade[allData.means.scaled$Clade == "fern"] = "ferns"
write.table(allData.means.scaled, file = 'allData.means.scaled.csv', append = F, quote = F, sep = ",", col.names = T, row.names = F)

##ok, making correlations plots, remember to save FREQUENTLY
#Correlation plots to make:
#1. cell vs leaf AR
#2. cell area vs solidity
#3. leaf area vs solidity
#4. cell area vs AR
#5. cell AR vs solidity

#cell AR vs Solidity
pdf("Fig4_CellARvsCellSolidity.pdf",h=15,w=20)
ggplot(allData.means.scaled[!is.na(allData.means.scaled$Clade),],aes(x=AR_Mean,y=Sol_Mean)) + geom_jitter(aes(color = Clade, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Cell Aspect Ratio Mean",y="Cell Solidity Mean") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()

model <- lm(Sol_Mean ~ AR_Mean, allData.means.scaled[!is.na(allData.means.scaled$Clade),])
summary(model)

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Sol_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"AR_Mean"],use="complete.obs",method="pearson")

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Sol_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"AR_Mean"],use="complete.obs",method="spearman")


#cell area vs solidity
pdf("Fig4_CellAreavsCellSolidity.pdf",h=15,w=20)
ggplot(allData.means.scaled[!is.na(allData.means.scaled$Clade),],aes(x=Scaled_Area_Mean,y=Sol_Mean)) + geom_jitter(aes(color = Clade, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Cell Area Mean",y="Cell Solidity Mean") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()

model <- lm(Sol_Mean ~ Scaled_Area_Mean, allData.means.scaled[!is.na(allData.means.scaled$Clade),])
summary(model)

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Sol_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Scaled_Area_Mean"],use="complete.obs",method="pearson")

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Sol_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Scaled_Area_Mean"],use="complete.obs",method="spearman")

#cell area vs aspect ratio
pdf("Fig4_CellAreavsCellAR.pdf",h=15,w=20)
ggplot(allData.means.scaled[!is.na(allData.means.scaled$Clade),],aes(x=Scaled_Area_Mean,y=AR_Mean)) + geom_jitter(aes(color = Clade, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Cell Area Mean",y="Cell Aspect Ratio Mean") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()

model <- lm(AR_Mean ~ Scaled_Area_Mean, allData.means.scaled[!is.na(allData.means.scaled$Clade),])
summary(model)

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"AR_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Scaled_Area_Mean"],use="complete.obs",method="pearson")

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"AR_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Scaled_Area_Mean"],use="complete.obs",method="spearman")

##read in leaf data and calculate leaf area, solidity, and AR

dataDir = "G:/My\ Drive/TeamEpidermis/DataAnalysis/LEAVES/TXT"
list.files(path = dataDir, full = T, pattern = "txt") %>% import_txt() %>% Out()-> rv_leaf

dataDir = "G:/My\ Drive/TeamEpidermis/Manuscript/Ferns/Leaf\ Outlines/"
list.files(path = dataDir, full = T, pattern = "jpg") %>% import_jpg() %>% Out()-> gp_leaf

allLeafCoo <- combine(rv_leaf,gp_leaf)
measures <- matrix(data = NA, nrow = length(allLeafCoo), ncol = 4)
for(i in 1:length(allLeafCoo))
{
	measures[i,] <- c(coo_area(allLeafCoo[i]), coo_eccentricityboundingbox(allLeafCoo[i]), coo_solidity(allLeafCoo[i]), coo_circularitynorm(allLeafCoo[i]) )
}

LeafMeasuresDF <- data.frame(measures)
colnames(LeafMeasuresDF) <- c('Area', 'AspectRatio', 'Solidity', 'Circularity')
LeafMeasuresDF$Sample_Number <- names(allLeafCoo)

LeafMeasuresDF$Sample_Number <- sapply(strsplit(LeafMeasuresDF$Sample_Number,split=" "),"[",1)
LeafMeasuresDF$Sample_Number <- gsub(pattern="-leaf|-ab-leaf|-ad-leaf",rep="",LeafMeasuresDF$Sample_Number)

#Calcuate means
Species <- unique(LeafMeasuresDF$Sample_Number)
meanTable <- matrix(data = NA, nrow = length(Species), ncol = 3)
for(i in 1:length(Species)){
	tempDF <- LeafMeasuresDF[LeafMeasuresDF$Sample_Number %in% Species[i],]
	#I know this isn't the best way to do it, but I think it is more clear what is happening
	solAllMean <- mean(tempDF$Solidity)
	arAllMean <- mean(tempDF$AspectRatio)
	areaAllMean <- mean(tempDF$Area)
	meanTable[i,] = c(solAllMean, arAllMean,areaAllMean)
}
meanLeafTableDF <- data.frame(meanTable)
meanLeafTableDF$Species_Code <- Species
colnames(meanLeafTableDF) <- c("Sol_Mean","AR_Mean","Area_Mean","Species_Code")
speciesKey <- AllmeasuresDF[,c("Species_Code","Clade")]
head(meanLeafTableDF <- merge(meanLeafTableDF,unique(speciesKey),by="Species_Code",all.x=TRUE,all.y=FALSE))
write.table(meanLeafTableDF, file = 'allLeaf.means.scaled.csv', append = F, quote = F, sep = ",", col.names = T, row.names = F)
colnames(meanLeafTableDF) <- c("Species_Code","Leaf_Sol_Mean","Leaf_AR_Mean","Leaf_Area_Mean","Clade")

head(allData.means.scaled);head(meanLeafTableDF)
head(cellAndLeafMeansDF<-merge(allData.means.scaled,meanLeafTableDF,by="Species_Code"))

#Make cell AR vs leaf AR plot
pdf("Fig4_CellARvsLeafAR.pdf",h=15,w=20)
ggplot(cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Leaf_AR_Mean < 1,],aes(x=AR_Mean,y=Leaf_AR_Mean)) + geom_jitter(aes(color = Clade.x, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Cell Aspect Ratio Mean",y="Leaf Aspect Ratio") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()
model <- lm(Leaf_AR_Mean ~ AR_Mean, cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Leaf_AR_Mean < 1,])
summary(model)

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Leaf_AR_Mean < 1,"AR_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Leaf_AR_Mean < 1,"Leaf_AR_Mean"],use="complete.obs",method="pearson")

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Leaf_AR_Mean < 1,"AR_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Leaf_AR_Mean < 1,"Leaf_AR_Mean"],use="complete.obs",method="spearman")

#Make leaf area vs solidity plot
pdf("Fig4_LeafAreavsLeafSolidity.pdf",h=15,w=20)
ggplot(cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),],aes(x=Leaf_Area_Mean,y=Leaf_Sol_Mean)) + geom_jitter(aes(color = Clade.x, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Leaf Area",y="Leaf Solidity") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()

model <- lm(Leaf_Sol_Mean ~ Leaf_Area_Mean, cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),])
summary(model)

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Leaf_Area_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Leaf_Sol_Mean"],use="complete.obs",method="pearson")

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Leaf_Area_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Leaf_Sol_Mean"],use="complete.obs",method="spearman")

#cell area vs aspect ratio and cell area vs solidity with log area
head(allData.means.scaled)
allData.means.scaled$Log_10_Area_Mean <- log10(allData.means.scaled$Scaled_Area_Mean)

#plot log cell area vs aspect ratio
pdf("Fig4_LogCellAreavsCellAR.pdf",h=15,w=20)
ggplot(allData.means.scaled[!is.na(allData.means.scaled$Clade),],aes(x=Log_10_Area_Mean,y=AR_Mean)) + geom_jitter(aes(color = Clade, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Log Cell Area Mean",y="Cell Aspect Ratio Mean") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()

model <- lm(AR_Mean ~ Log_10_Area_Mean, allData.means.scaled[!is.na(allData.means.scaled$Clade),])
summary(model)

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"AR_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Log_10_Area_Mean"],use="complete.obs",method="pearson")

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"AR_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Log_10_Area_Mean"],use="complete.obs",method="spearman")

#plot log cell area vs solidity
pdf("Fig4_LogCellAreavsCellSolidity.pdf",h=15,w=20)
ggplot(allData.means.scaled[!is.na(allData.means.scaled$Clade),],aes(x=Log_10_Area_Mean,y=Sol_Mean)) + geom_jitter(aes(color = Clade, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Log10 Cell Area Mean",y="Cell Solidity Mean") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()

model <- lm(Sol_Mean ~ Log_10_Area_Mean, allData.means.scaled[!is.na(allData.means.scaled$Clade),])
summary(model)

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Sol_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Log_10_Area_Mean"],use="complete.obs",method="pearson")

cor.test(x=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Sol_Mean"],y=allData.means.scaled[!is.na(allData.means.scaled$Clade),"Log_10_Area_Mean"],use="complete.obs",method="spearman")

#need to do log leaf area
head(cellAndLeafMeansDF)
cellAndLeafMeansDF$Log_10_Leaf_Area_Mean <- log10(cellAndLeafMeansDF$Leaf_Area_Mean)

#Make log10 leaf area vs solidity plot
pdf("Fig4_LogLeafAreavsLeafSolidity.pdf",h=15,w=20)
ggplot(cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),],aes(x=Log_10_Leaf_Area_Mean,y=Leaf_Sol_Mean)) + geom_jitter(aes(color = Clade.x, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Log10 Leaf Area",y="Leaf Solidity") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#77A340", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE,shape = FALSE, size = FALSE)
dev.off()

model <- lm(Leaf_Sol_Mean ~ Log_10_Leaf_Area_Mean, cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),])
summary(model)

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Log_10_Leaf_Area_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Leaf_Sol_Mean"],use="complete.obs",method="pearson")

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Log_10_Leaf_Area_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x),"Leaf_Sol_Mean"],use="complete.obs",method="spearman")

#Make cell AR vs leaf AR plot w/o eudicots
pdf("Fig4_CellARvsLeafAR_noEudicots.pdf",h=15,w=20)
ggplot(cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Leaf_AR_Mean < 1 & cellAndLeafMeansDF$Clade.x != "eudicots",],aes(x=AR_Mean,y=Leaf_AR_Mean)) + geom_jitter(aes(color = Clade.x, shape="circle",size=10)) + geom_smooth(method="lm",se=FALSE, colour="black") + 
theme_minimal() + labs(x="Cell Aspect Ratio Mean",y="Leaf Aspect Ratio") + 
scale_color_manual(breaks = c("ferns", "gymnosperms", "early_diverging_angiosperms", "eudicots", "monocots"), values=c("#369B9A", "#304EA1", "#562D51", "#D93B95")) +
theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30)) + guides(colour=FALSE, shape = FALSE, size = FALSE)
dev.off()
model <- lm(Leaf_AR_Mean ~ AR_Mean, cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Clade.x != "eudicots" & cellAndLeafMeansDF$Leaf_AR_Mean < 1,])
summary(model)

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Clade.x != "eudicots" & cellAndLeafMeansDF$Leaf_AR_Mean < 1,"AR_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Clade.x != "eudicots","Leaf_AR_Mean"],use="complete.obs",method="pearson")

cor.test(x=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Clade.x != "eudicots" & cellAndLeafMeansDF$Leaf_AR_Mean < 1,"AR_Mean"],y=cellAndLeafMeansDF[!is.na(cellAndLeafMeansDF$Clade.x) & cellAndLeafMeansDF$Clade.x != "eudicots" & cellAndLeafMeansDF$Leaf_AR_Mean < 1,"Leaf_AR_Mean"],use="complete.obs",method="spearman")

#Generate cell PCAs (do this using Momocs code to get pretty figures)
#need to get circularitynorm() results 
l_circnorm<-measure(l,coo_circularitynorm)
gpcells_circnorm<-measure(gp_cells,coo_circularitynorm)
rvCircNorm <- data.frame(l_circnorm$coe)
rvCircNorm$Sample_Number <- rownames(rvCircNorm)
gpCircNorm <- data.frame(gpcells_circnorm$coe)
gpCircNorm$Sample_Number <- rownames(gpCircNorm)
CircNormDF <- rbind(rvCircNorm,gpCircNorm)

#Ok, need to get all the data together into a TraCoe
allData.measures.scaled <- read.csv("all.measures.scaled.csv")
allData.measures.scaled <- merge(allData.measures.scaled,CircNormDF)
allData.measures.scaled$Clade[allData.measures.scaled$Clade == "fern"] = "ferns"
cellCoeInputDF <- allData.measures.scaled[!is.na(allData.measures.scaled$Clade),]

cellCoe <- TraCoe(cellCoeInputDF[,c(3,4,11,12)],fac=cellCoeInputDF[,c(2,8)])
#test PCA building
cellCoe %>% PCA() %>% plot(ellipsesax=F,contour=T,lev.contour=10,eigen=F,rug=F,labelsgroups=F,alpha=0.1,col="grey",col.contour="black")

#build PCA file for ggplot graph
cellCoe.pc <- PCA(cellCoe)
PCmeansBySpecies <- data.frame()
species <- unique(cellCoe.pc$fac[,1])
for(i in 1:length(species)){
	tempDF <- cellCoe.pc$x[grep(pattern = species[i], cellCoe.pc$fac[,1]),]
	tempRow <- colMeans(tempDF)
	PCmeansBySpecies <- rbind(PCmeansBySpecies,tempRow)
}
colnames(PCmeansBySpecies) <- c("PC1","PC2","PC3","PC4")
PCmeansBySpecies$Species_Code = species
PCmeansBySpecies <- merge(PCmeansBySpecies,unique(allData.measures.scaled[,c(2,8)]),all.x=T,all.y=F)

pdf("PCA_alltaxa.pdf",w=20,h=20)
cellCoe %>% PCA() %>% as_df() %>% ggplot +
aes(x=PC1,y=PC2) + coord_equal() +
geom_point(colour="darkgrey",alpha=0.25,size=5) + 
geom_density2d(colour="#000000",bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + geom_point(data=PCmeansBySpecies, aes(x=PC1,y=PC2),shape=18,fill="darkgrey",size=5)
dev.off()
#need to do a separate figure with loadings
pdf("cellPCAloadings.pdf",width=12,height=12)
cellCoe %>% PCA() %>% plot(loadings=T,points=F,eigen=F,rug=F)
dev.off()
#get the PC percentages (easiest on the graph)
cellCoe %>% PCA() %>% plot()

#redo with color for each clade, and then only plotting eudicots and monocots
pdf("PCA_alltaxa_means_lightColorMatch.pdf",w=20,h=20)
cellCoe %>% PCA() %>% as_df() %>% ggplot +
aes(x=PC1,y=PC2,colour=Clade) + coord_equal() +
geom_point(alpha=0.25,size=5) + geom_density2d(bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + geom_point(data=PCmeansBySpecies, aes(x=PC1,y=PC2),shape=18,fill="darkgrey",size=5) +
scale_color_manual(breaks = c("ferns", "eudicots", "gymnosperms", "early_diverging_angiosperms", "monocots"), values=c("#97C8C7", "#B0C86C", "#8F9BCD", "#9B8498", "#E184B7")) +
guides(colour=FALSE)
dev.off()

#eudicots and monocots
pdf("PCA_alltaxa_eudicots_and_monocots_cladeColors.pdf",w=20,h=20)
cellCoe.pc.DF <- (cellCoe %>% PCA() %>% as_df())
ggplot(data = subset(cellCoe.pc.DF,Clade %in% c("eudicots","monocots"))) + aes(x=PC1,y=PC2,colour=Clade) + coord_equal() +
geom_point(alpha=0.25,size=5) + geom_density2d(bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + geom_point(data=subset(PCmeansBySpecies,Clade %in% c("eudicots","monocots")), aes(x=PC1,y=PC2),shape=18,fill="darkgrey",size=5) +
scale_color_manual(breaks = c("eudicots", "monocots"), values=c("green", "red")) +
guides(colour=FALSE) + scale_x_continuous(limits = c(-3.5, 7.5))
dev.off()

#redo with correct colors
#add variable for shape and color match
cellCoe.pc.DF$Shape = 1
PCmeansBySpecies$Shape = 2
cellCoe.pc.DF$Color = paste(cellCoe.pc.DF$Clade,cellCoe.pc.DF$Shape,sep="_")
PCmeansBySpecies$Color = paste(PCmeansBySpecies$Clade,PCmeansBySpecies$Shape,sep="_")
cellCoe.pc.DF <- cellCoe.pc.DF[,c(2,4,5,6,7,3,8,9)]
PCA.df <- rbind(cellCoe.pc.DF,PCmeansBySpecies)

pdf("PCA_alltaxa_brightColorMatch.pdf",w=20,h=20)
ggplot(data = subset(PCA.df,!is.na(Clade))) +
aes(x=PC1,y=PC2,group = Clade,colour=Color,shape=Shape) + coord_equal() +
geom_jitter(size=5) + geom_density2d(aes(col = Clade),bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + xlim(-3.5,8) + ylim(-3,5) +
scale_color_manual(values=c("#369B9A", "#97C8C7", "#369B9A", "#77A340", "#A9C065", "#77A340", "#304EA1", "#8E9ACB", "#304EA1", "#562D51", "#9C8599", "#562D51", "#D93B95", "#E184B7", "#D93B95")) +
scale_shape_manual(values = c(16,18)) + guides(colour=FALSE,shape=FALSE)
dev.off()

#make graphs for each clade
pdf("PCA_eudicots_brightColorMatch.pdf",w=20,h=20)
ggplot(data = subset(PCA.df,Clade == "eudicots")) +
aes(x=PC1,y=PC2,group = Clade,colour=Color,shape=Shape) + coord_equal() +
geom_jitter(size=5) + geom_density2d(aes(col = Clade),bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + xlim(-3.5,8) + ylim(-3,5) +
scale_color_manual(values=c("#77A340", "#A9C065", "#77A340")) + 
scale_shape_manual(values = c(16,18)) + guides(colour=FALSE,shape=FALSE)
dev.off()

pdf("PCA_monocots_brightColorMatch.pdf",w=20,h=20)
ggplot(data = subset(PCA.df,Clade == "monocots")) +
aes(x=PC1,y=PC2,group = Clade,colour=Color,shape=Shape) + coord_equal() +
geom_jitter(size=5) + geom_density2d(aes(col = Clade),bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + xlim(-3.5,8) + ylim(-3,5) +
scale_color_manual(values=c("#D93B95","#E184B7","#D93B95")) + 
scale_shape_manual(values = c(16,18)) + guides(colour=FALSE,shape=FALSE)
dev.off()

pdf("PCA_earlydivergingangiospers_brightColorMatch.pdf",w=20,h=20)
ggplot(data = subset(PCA.df,Clade == "early_diverging_angiosperms")) +
aes(x=PC1,y=PC2,group = Clade,colour=Color,shape=Shape) + coord_equal() +
geom_jitter(size=5) + geom_density2d(aes(col = Clade),bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + xlim(-3.5,8) + ylim(-3,5) +
scale_color_manual(values=c("#369B9A","#97C8C7","#369B9A")) + 
scale_shape_manual(values = c(16,18)) + guides(colour=FALSE,shape=FALSE)
dev.off()

pdf("PCA_gymnosperms_brightColorMatch.pdf",w=20,h=20)
ggplot(data = subset(PCA.df,Clade == "gymnosperms")) +
aes(x=PC1,y=PC2,group = Clade,colour=Color,shape=Shape) + coord_equal() +
geom_jitter(size=5) + geom_density2d(aes(col = Clade),bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + xlim(-3.5,8) + ylim(-3,5) +
scale_color_manual(values=c("#562D51","#9C8599","#562D51")) + 
scale_shape_manual(values = c(16,18)) + guides(colour=FALSE,shape=FALSE)
dev.off()

pdf("PCA_ferns_brightColorMatch.pdf",w=20,h=20)
ggplot(data = subset(PCA.df,Clade == "ferns")) +
aes(x=PC1,y=PC2,group = Clade,colour=Color,shape=Shape) + coord_equal() +
geom_jitter(size=5) + geom_density2d(aes(col = Clade),bins=10) + 
theme_minimal() + geom_hline(yintercept=0,col="black") +
geom_vline(xintercept=0,col="black") + xlim(-3.5,8) + ylim(-3,5) +
scale_color_manual(values=c("#304EA1","#8E9ACB","#304EA1")) + 
scale_shape_manual(values = c(16,18)) + guides(colour=FALSE,shape=FALSE)
dev.off()

##Remake table S1 to match desired format
##start with allData.means.scaled
View(allData.means.scaled)
adData.means <- allData.means.scaled[!is.na(allData.means.scaled$Sol_Ad_Mean),c(1,15,3,7,10,13)]
abData.means <- allData.means.scaled[!is.na(allData.means.scaled$Sol_Ab_Mean),c(1,15,4,8,11,14)]
naData.means <- allData.means.scaled[is.na(allData.means.scaled$Sol_Ad_Mean) & 
	is.na(allData.means.scaled$Sol_Ab_Mean),c(1,15,2,6,9,12)]
colnames(adData.means) <- c("Species Code","Major Group", "Mean Solidity", "Mean Aspect Ratio", "Mean Area", "Mean Circularity") 
colnames(abData.means) <- c("Species Code","Major Group", "Mean Solidity", "Mean Aspect Ratio", "Mean Area", "Mean Circularity") 
colnames(naData.means) <- c("Species Code","Major Group", "Mean Solidity", "Mean Aspect Ratio", "Mean Area", "Mean Circularity")
adData.means$Side = "ad"
abData.means$Side = "ab"
naData.means$Side = "NA" 
#make table with all data
meansBySide <- rbind(adData.means,abData.means,naData.means)

#change to Log10 Measures to match data from Siobhan
head(meansBySide.log)
meansBySide$"Log Mean Area" = log10(meansBySide$"Mean Area")
meansBySide$"Log Mean Circularity" = log10(meansBySide$"Mean Circularity")
meansBySide.log <- meansBySide[,c(1,2,7,8,4,9,3)]
#Now add Leaf Data by Species Code Column
#I am just going to attach the first leaf measure, and then (in the rare case they were different
#leaves) change the leaf measures manually
#Also don't worry about 9.* right now, will add manually later
LeafMeasuresDF$LogArea = log10(LeafMeasuresDF$Area)
LeafMeasuresDF$LogCirc = log10(LeafMeasuresDF$Circularity)
LeafMeasuresDF.log = LeafMeasuresDF[,c(5,9,2,7,3)]
colnames(LeafMeasuresDF.log) = c("Species Code","Leaf Log Area","Leaf Aspect Ratio","Leaf Log Circularity","Leaf Solidity")
suppTable1.precursor <- merge(meansBySide.log, LeafMeasuresDF.log, by="Species Code", all = TRUE)
View(suppTable1.precursor)

#Forgot to scale the leaves
LeafMeasuresDF$Scaled_Area = LeafMeasuresDF$Area * (1/118)^2
LeafMeasuresDF$LogScaledArea = log10(LeafMeasuresDF$Scaled_Area)
head(LeafMeasuresDF)

#redid merge, now save before continuing
write.csv(suppTable1.precursor,file  = "SupplTable1.precursor.csv")

#Redo leaf area vs leaf solidity
#Luckily everything else should be the same
#And since area is scaled by a factor, shouldn't affect the p-values for correlations hopefully
