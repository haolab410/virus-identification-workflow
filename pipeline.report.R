library(tools)
library(brew)
library(optparse)

args <- commandArgs(TRUE)

folder<-args[1]
sample <- args[2]
c <- as.integer(args[3])

path<-paste(folder,sample,"",sep="")
setwd(folder)

qcs <- ''
for (i in 1:c){
    qc <- system(paste("sed 's/count/", i, "/' qc.brew", sep=''),intern=T)
    writeLines(qc, paste(i, ".qcs.brew", sep=''))
    }

system(paste("cat ",folder,"/title.brew ", folder, "/*.qcs.brew ", folder, "/section2.brew ", folder, "/end.brew >", folder, "/temp.brew",sep=""), intern=F)
label<-paste('<% cat(sample<-"',sample,'") %>',sep="")

file.brew<-system(paste("sed 's/<% cat(sample<-\".*\") %>/",label,"/' temp.brew",sep=""),intern=T)

writeLines(file.brew,"reportofvirus.brew")
brew("reportofvirus.brew", "report.tex",text=sample)
texi2dvi("report.tex", pdf = TRUE)
 