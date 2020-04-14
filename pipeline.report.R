library(tools)
library(brew)
library(optparse)

args <- commandArgs(TRUE)

folder<-args[1]
sample <- args[2]
  path.cur <- getwd()
  print(path.cur)
  path<-paste(folder,sample,"",sep="")
  setwd(folder)
  print(path)
  system(paste("cp ",folder,"/reportofvirus.brew temp.brew",sep=""), intern=F)
     
    label<-paste('<% cat(sample<-"',sample,'") %>',sep="")
    print(label)
file.brew<-system(paste("sed 's/<% cat(sample<-\".*\") %>/",label,"/' temp.brew",sep=""),intern=T)
    print(file.brew)
    writeLines(file.brew,"reportofvirus.brew")
    
brew("reportofvirus.brew", "report.tex",text=sample)
texi2dvi("report.tex", pdf = TRUE)
 