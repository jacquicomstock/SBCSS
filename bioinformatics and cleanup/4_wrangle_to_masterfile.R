#load data
RA <- read.csv("SBCSS_seqtab-nochimtaxa_rar13000_NoSingle.csv", header=TRUE,sep=",",row.names=1)
meta <- read.csv("SBCSS__meta_rar13000.csv", header=TRUE,sep=",", row.names = 1)
taxa <- read.csv("SBCSS_taxa_rar13000_NoSingle.csv", header=TRUE,sep=",", row.names = 1, na.strings = "NOPE")

#Create column in taxa file with ascending number (for unique ASV identification)
num <- data.frame(x = c(1:nrow(taxa)))
taxa_num <- as.data.frame(cbind(taxa, num))

#concatinate taxa names with unique number identifier
taxa$fullname <- paste(taxa$Kingdom, taxa$Phylum, taxa$Class, taxa$Order, taxa$Family, taxa$Genus, taxa$Species, taxa_num$x, sep = "_")

#use these concatinated names as row names for RA table & add metadata
row.names(RA) <- taxa$fullname
tRA <- as.data.frame(t(RA))
fulldf <- as.data.frame(cbind(meta, tRA))

#write file
write.csv(fulldf, file="SBCSS_EnvASV_rar13000.csv")
