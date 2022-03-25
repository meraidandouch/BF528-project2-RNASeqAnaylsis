#Reading in data
FPKM<-read.table("/projectnb/bf528/users/hedgehog_2022/Project2/Programmer-Project2/P0_1_tophat/P0_1_cufflinks/genes.fpkm_tracking")
#Data formatting
colnames(FPKM)<-FPKM[1,]
FPKM<-FPKM[-1,]
Num<-as.numeric(FPKM$FPKM)
#filtering fpkm values lower than 1
fpkm<-Num[Num>=1]
#histogram
hist(log10(fpkm),xlim=c(0,7),col = "#FF9933", border = 'black',
     main = paste("Histogram of log10 FPKM values greater than 0"),
     xlab = "fpkm values")
hist(log10(Num [Num != 0]),col = "#FF9933", border = 'black',
     main = paste("Histogram of log10 FPKM values other than 0"),
     xlab = "fpkm values")



