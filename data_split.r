library("caret")

args <- commandArgs(trailingOnly=TRUE)
indir <- file.path(getwd(), args[1])
outdir <- file.path(getwd(), args[2])

setwd(indir)
input_files <- list.files(path = ".")

for(i in 1:length(input_files)){
	data <- read.csv(input_files[i], header=TRUE, sep=",")
	name <- strsplit(input_files[i], "\\.")[[1]][1]
	
	train <- createDataPartition(y=data$value, p=.80, list=FALSE)
	test <- data[-train,]
	train <- data[train,]
	write.table(test$id, paste(outdir, name, "_test.csv", sep=""), row.names=FALSE, col.names=FALSE)
	write.table(train$id, paste(outdir, name, "_train.csv", sep=""), row.names=FALSE, col.names=FALSE)
}