#' This is MakeCSV4Plot
#' @param stuff and more
#' @keywords multi
#' @export
#' @examples
#' #not run


MakeCSV4Plot<-function( .df,  FileName ="test", PlotNames=OfficialNames){
.m<-merge( OfficialNames,.df,all=TRUE,sort = FALSE) #
.me<-na.omit(data.frame("Filename"=.m[,2], "Outcome"=.m[,3]))
rownames(.me)<-.me[,1]
.me2<-as.data.frame(t(.me[,2]))
colnames(.me2)<-rownames(.me)
.me3<-cbind(data.frame("Filename"="Outcome"), .me2)
write.csv( .me3,file=paste(FileName, ".csv", sep=""),row.names=FALSE)
}
