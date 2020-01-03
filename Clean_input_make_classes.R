library(data.table)
library(argparse)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--samples", help="sample list")
parser$add_argument("--pop", help="population codes for references")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

samples<-fread(args$samples,stringsAsFactors = F,header=F)
pop<-fread(args$pop,stringsAsFactors = F,header = F)
populations<-unique(pop$V2)

joined<-left_join(samples,pop,by=c("V1"))
joined[is.na(joined)]<-"adm"

joined$V2<-as.numeric(factor(joined$V2,levels=c("adm",populations)))
joined$V2<-(joined$V2 - 1)
joined$V3<-joined$V2
joined$V1<-NULL
classes<-as.vector(unlist(t(joined),use.names = F))
cat(classes, file = args$out)

