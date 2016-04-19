## Microbial co-occurrence from temporally resolved 16S data

### Data
The data for this project are from David *et al.* 2014. Genome Biol. paper.
Downloaded from NCBI using fastq-dump from sratoolkit. BioProject PRJEB6518

``` {r}
setwd("/Users/okamneva/Projects/David14") ## CHANGE THIS to your working directory 
samples=read.table("SraRunTable.txt", sep="\t", head=T, stringsAsFactors=F) ## read in table of samples
dim(samples)
samples=samples[samples$description_s=="DonorA Stool",]
dim(samples)
samples=samples[!is.na(as.numeric(samples$collection_day_s)),]
dim(samples)
samples=samples[order(as.numeric(samples$collection_day_s)),]
sratk="~/Downloads/sratoolkit.2.5.4-mac64/" ## CHANGE THIS to the location of sra tool kit on your machine
commands=paste(sratk, "bin/fastq-dump --split-files ", samples$Run_s, sep="")
if(!dir.exists("fastq")){dir.create("fastq")}
write(commands, "fastq/commands_get_reads")
```

``` shell
cd /Users/okamneva/Projects/David14/fastq
source commands_get_reads
```


### Checking read quality using ShortRead R package





```
dataDir="/Users/okamneva/Projects/David14" ### 
fastqDir <- file.path(dataDir, "fastq")
fls <- list.files(fastqDir, "fastq$", full=TRUE)


postscript(paste(dataDir, "/fastQC.eps", sep=""))
summary=c()
for(i in 1:length(fls))
{
	qaSummary <- qa(fls[i], type="fastq")
	perCycle <- qaSummary[["perCycle"]]
	ShortRead:::.plotCycleQuality(perCycle$quality)
	summary=rbind(summary, c(qaSummary[["readCounts"]]$read,paste(qaSummary[["adapterContamination"]][1], collapse=", ") ))
}
dev.off()
```



* item1
* item2
* item3