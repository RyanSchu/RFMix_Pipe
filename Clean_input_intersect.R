##Take the intersection of snps in each file

library(data.table)
library(argparse)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--ref", help="reference snp list")
parser$add_argument("--query", help="query snp list")
parser$add_argument("--map", help="genetic map")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

names1<-c("ID","chr","bp")
names2<-c("ID","bp","cm")

refSNPS<-fread(args$ref, skip = 1, header = F)
querySNPS<-fread(args$query, skip = 1, header = F)
mapSNPS<-fread(args$map, skip = 1, header = F)

colnames(refSNPS)<-names1
colnames(querySNPS)<-names1
colnames(mapSNPS)<-names2

refSNPS<-refSNPS %>% distinct(ID, .keep_all = TRUE) %>% distinct(bp, .keep_all = TRUE)
querySNPS<-querySNPS  %>% distinct(ID, .keep_all = TRUE) %>% distinct(bp, .keep_all = TRUE)
mapSNPS<-mapSNPS %>% distinct(ID, .keep_all = TRUE) %>% distinct(bp, .keep_all = TRUE)

ID_intersect<-base::intersect(refSNPS$ID,querySNPS$ID) %>% base::intersect(mapSNPS$ID)
bp_intersect<-base::intersect(refSNPS$bp,querySNPS$bp) %>% base::intersect(mapSNPS$bp)
nIDs<-length(ID_intersect)
nBPs<-length(bp_intersect)

if(nIDs != nBPs){
  cat("Warning: Number of overlapping SNP IDs != to number of overlapping base pair positions\n")
  cat("Number of overlapping SNP IDs: ", nIDs,"\n")
  cat("Number of overlapping base pair positions: ", nBPs,"\n")
  cat("Files may be of different builds, contain multiple chromosomes, or have different SNP ID types.\n")
  if (nIDs ==0 | nBPs ==0) {
    cat("One or more list is empty. Exiting.\n")
  }
  if(nIDs/nBPs <= 0.90 | nIDs/nBPs >= 1/(0.90)){
    cat("Intersections less than 90% concordant in length. Please resolve above before continuing.\nExiting.\n")
    quit(1)
  } else{
    cat("Majority of SNPS overlap. Continuing with smallest number of snps\n")
  }  
}

final_list<-inner_join(refSNPS,querySNPS,by=c("ID","chr","bp")) %>% inner_join(mapSNPS,by=c("ID","bp"))
cat("Length of final overlap: ", dim(final_list)[1],"\n")
fwrite(final_list,args$out %&% "_intersect_snp.list.txt",sep="\t",col.names=F)
