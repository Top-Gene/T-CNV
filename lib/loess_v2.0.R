#!Rscript
library('getopt')

spec = matrix(c(
          'help', 'h',  0,   "logical",
        'infile', 'i',  1, "character",
       'prefix', 'p',  1, "character"
        ), byrow=TRUE, ncol=4);
opt = getopt(spec)

#------------------------------------------------------
# define usage function				      
#------------------------------------------------------
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE))
  cat("
Usage example:
        Rscript *.R --infile tumor_VS_normal.log2 --outfile tumor_VS_normal.adj.xls
Options:
        --help          -h      NULL             get the help
        --infile        -i      character        the input file [forced]
        --prefix        -p      character        the output filename of picture  [forced]
\n")
        q(status=1)
}

dealData <- function(log2score,gcScore){

	pdf(paste(opt$prefix,".pdf",sep=""),width=12,height=6)
	par(mar = c(5,5,4,5)+0.1,mfrow=c(1,2))

	l = loess(log2score ~ gcScore)
	mc <- predict(l, gcScore)
	adj = log2score-mc

	plot(gcScore,log2score,ylim=c(min(log2score)-1,max(log2score)+1),xlab="gc content",ylab="log2",cex.lab=1.5,type="p",pch=18,lwd=2,col="purple4")
	plot(gcScore,adj,ylim=c(min(adj)-1,max(adj)+1),xlab="gc content",ylab="adj_log2",cex.lab=1.5,type="p",pch=18,lwd=2,col="purple4")

	dev.off()

	return(adj)
}

#--------------------FUNCTION END-----------------------

## --help或参数错误 打印帮助文档
if ( !is.null(opt$help) || is.null(opt$infile) || is.null(opt$prefix) ) { 
	print_usage(spec) 
}

content = read.table(opt$infile,header=F,sep="\t")

## 处理样本log2
SampleAdj = dealData(content[,9],content[,6])

result <- data.frame(content,SampleAdj)
finalFile <- paste(opt$prefix,".result.txt",sep="")
write.table (result,file = finalFile, row.names = FALSE, col.names =FALSE, quote =FALSE,sep="\t")

# 在文件中添加表头
#system(paste("sed -i '1iInterval\tGene\tExon\ttotal_tumor\tmean_tumor\ttotal_normal\tmean_normal\tgc_content\tlog2_sample\tlog2_gene\tadj_sam\tadj_gene\n' ",finalFile))
