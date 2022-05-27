file = "~/Google_Drive/Projects/Bashford-Rogers Lab/Lab reports/Lauren Overend/All_raw_values_SEPSIS_BCR1_1500_PRODUCTIVE (1).txt"
p <- as.matrix(read.csv(file, head=T, sep="\t"))
colnames(p)

p1 = p[,which(colnames(p)!="ReadDepth" )]
samples = p1[,1]
p1 = p1[,c(2:length(p1[1,]))]
headers = colnames(p1)
mat = matrix(data = -1, nrow = length(samples), ncol = length(headers), dimnames = c(list(samples), list(headers)))
for(i in c(1:length(headers))){
	mat[,headers[i]] = as.numeric(p1[, headers[i]])
}

## remove all incomplete rows then columns
nz_cols = apply(mat, 2, function(x){length(which(x!=-1))})
matr1 = mat[,which(nz_cols>= max(nz_cols))]
nz_rows = apply(mat1, 1, function(x){length(which(x!=-1))})
matr2 = mat1[which(matr1>= max(nz_rows)),]


###########
fileout1='~/Google_Drive/Projects/Bashford-Rogers Lab/Lab reports/Lauren Overend/Metrics_anlaysis_BCR_1.pdf'
w=3.35
pdf(file=fileout1, height=w*2, width=w*2)
par(mfrow= c(2,2), mar = c(5,5,3,3))

nz_cols = apply(mat, 2, function(x){length(which(x!=-1))})
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main ="run1")
threshold = quantile(nz_cols, 0.25)
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat1 = mat[,which(nz_cols>= threshold)]


nz_rows = apply(mat1, 1, function(x){length(which(x!=-1))})
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main ="run2")
threshold = quantile(nz_rows, 0.1)
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)

mat2 = mat1[which(nz_rows>= threshold),]

nz_cols = apply(mat2, 2, function(x){length(which(x!=-1))})
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main ="run3")
threshold = quantile(nz_cols, 0.25)
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat3 = mat2[,which(nz_cols>= threshold)]


nz_rows = apply(mat3, 1, function(x){length(which(x!=-1))})
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main ="run4")
threshold = quantile(nz_rows, 0.1)
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)

mat4 = mat3[which(nz_rows>= threshold),]
dev.off()
dim(mat)
dim(mat1)
dim(mat2)
dim(mat3)
dim(mat4)
dim(matr2)


cbind(colnames(mat4))

#################################
file = "~/Google_Drive/Projects/Bashford-Rogers Lab/Lab reports/Lauren Overend/isotyper_metrics_1500_PRODUCTIVE (1).txt"
p <- as.matrix(read.csv(file, head=T, sep="\t"))
colnames(p)
p=p[which(p[,"Include.metric...1..use."]=="1"),]
metric = p[,"Metric"]
mat_filtered = mat[,metric]

###########
fileout1='~/Google_Drive/Projects/Bashford-Rogers Lab/Lab reports/Lauren Overend/Metrics_anlaysis_BCR_2.pdf'
w=3.35
pdf(file=fileout1, height=w*2, width=w*2)
par(mfrow= c(2,2), mar = c(5,5,3,3))

nz_rows = apply(mat_filtered, 1, function(x){length(which(x!=-1))})
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main ="run2.1")
threshold = quantile(nz_rows, 0.1)
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat_filtered1 = mat_filtered[which(nz_rows>= threshold),]


nz_cols = apply(mat_filtered1, 2, function(x){length(which(x!=-1))})
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main ="run2.2")
# threshold = quantile(nz_cols, 0.25)
# segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
# mat1 = mat[,which(nz_cols>= threshold)]

dev.off()


dim(mat_filtered)
dim(mat_filtered1)





