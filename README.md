### Checking read quality using ShortRead R package. 

### Data
The data for this project are from David *et al.* 2014. Genome Biol. paper.
Downloaded from NCBI using fastq-dump from sratoolkit. BioProject PRJEB6518

```{r}
dataDir="/Users/okamneva/Projects/David14"
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