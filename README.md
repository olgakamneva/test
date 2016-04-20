## Microbial co-occurrence from temporally resolved 16S data

### Getting data
The data for this project are from David *et al.* 2014. Genome Biol. paper. BioProject PRJEB6518.
File SraRunTable.txt was downloaded from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP006059 and contains summary information about each run in the study.
First I look at the list of samples in R and create command file to download fastq files from NCBI using fastq-dump from sratoolkit. 
There are a lot of samples in this study, I will only include stool samples from one donor for now.


``` {r}
## In R
setwd("/Users/okamneva/Projects/David14") ## CHANGE THIS to your working directory 
samples=read.table("SraRunTable.txt", sep="\t", head=T, stringsAsFactors=F) ## read in table of samples
names(samples)
table(samples$description_s)
table(samples$LibraryLayout_s)
samples=samples[samples$description_s=="DonorA Stool",]
dim(samples)
samples=samples[!is.na(as.numeric(samples$collection_day_s)),]
dim(samples)
samples=samples[order(as.numeric(samples$collection_day_s)),]
sratk="~/Downloads/sratoolkit.2.5.4-mac64/" ## CHANGE THIS to the location of sra tool kit on your machine
commands=paste(sratk, "bin/fastq-dump ", samples$Run_s, sep="")
if(!dir.exists("fastq")){dir.create("fastq")}
write(commands, "fastq/commands_get_reads")
write.table(samples,"SraRunTable_filterred.txt" ,quote =F, sep ="\t", row.names = F, col.names = TRUE)
```

``` shell
## In terminal
cd /Users/okamneva/Projects/David14/fastq
source commands_get_reads
head ERR531447.fastq
```
It it

### Checking read quality using ShortRead R package


```{r}
## In R
setwd("/Users/okamneva/Projects/David14") ## CHANGE THIS to your working directory 
fastqDir <- file.path("fastq")
fls <- list.files(fastqDir, "fastq$", full=TRUE)

library(ShortRead)
summary=c()
for(i in 1:length(fls))
{
	cat(i, "\n")
	qaSummary <- qa(fls[i], type="fastq")
	perCycle <- qaSummary[["perCycle"]]
	perCycle=split(perCycle$quality, as.factor(perCycle$quality[,1]))
	
	means=unlist(lapply(perCycle, function(x) mean(unlist(mapply(rep, x[,3], x[,4]))) ))
	medians=unlist(lapply(perCycle, function(x) median(unlist(mapply(rep, x[,3], x[,4]))) ))
	q25=unlist(lapply(perCycle, function(x) quantile(unlist(mapply(rep, x[,3], x[,4])), probs = 0.25) ))
	summary=rbind(summary, 
		c(
		qaSummary[["readCounts"]]$read,
		paste(qaSummary[["adapterContamination"]][1], collapse=",") ,
		paste(which(means<28), collapse=",") ,
		paste(which(medians<28), collapse=","),
		paste(which(q25<28), collapse=","),
		paste(which(means<25), collapse=",") ,
		paste(which(medians<25), collapse=","),
		paste(which(q25<25), collapse=","),
		paste(which(means<20), collapse=",") ,
		paste(which(medians<20), collapse=","),
		paste(which(q25<20), collapse=",")
		))
}
names(summary)=c("readCount","AdapterContamination","means28", "medians28", "lq28" ,"means25", "medians25", "lq25" ,"means20", "medians20", "lq20"  )
```


