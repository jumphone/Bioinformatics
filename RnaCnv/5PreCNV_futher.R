


REFERENCE_CELL="reference_cell.txt"
WINDOW=100

DATA='center_ALLCNV.txt'

scale_center=function(x){return(scale(x,center=T,scale=F))}
raw_data=read.delim(DATA,header=T,sep='\t')
center_ALLCNV=as.matrix(raw_data[,c(5:length(raw_data[1,]))])

reference_cell_name=as.character(read.delim(REFERENCE_CELL,header=F,sep='\t')[,1])
reference_cell= which(colnames(center_ALLCNV) %in% reference_cell_name )

#max_base_line = apply(center_ALLCNV[,reference_cell],1,max)
#min_base_line = apply(center_ALLCNV[,reference_cell],1,min)

max_base_line = apply(center_ALLCNV[,reference_cell],1,quantile,0.95)
min_base_line = apply(center_ALLCNV[,reference_cell],1,quantile,0.05)

#mean_base_line=apply(center_ALLCNV[,reference_cell],1,mean)
#sd_base_line=apply(center_ALLCNV[,reference_cell],1,sd)

#max_base_line =qnorm(0.975, mean=mean_base_line, sd=sd_base_line)      #apply(center_ALLCNV[,reference_cell],1,mean)
#min_base_line =qnorm(0.025, mean=mean_base_line, sd=sd_base_line) #apply(center_ALLCNV[,reference_cell],1,mean)


relative_center_ALLCNV = center_ALLCNV
i=1
while(i <= length(center_ALLCNV[,1])){
j=1
while(j <= length(center_ALLCNV[1,])){
#relative_center_ALLCNV[i,j] = center_ALLCNV[i,j] - max_base_line[i]	
	
if(center_ALLCNV[i,j] > max_base_line[i]  ){
relative_center_ALLCNV[i,j] = (center_ALLCNV[i,j] - max_base_line[i]) }
else if(center_ALLCNV[i,j] < min_base_line[i]  ){
relative_center_ALLCNV[i,j] = (center_ALLCNV[i,j] - min_base_line[i]) }
else {relative_center_ALLCNV[i,j] = 0}

j=j+1}
i=i+1
print(i)}

rownames(relative_center_ALLCNV)=rownames(raw_data[,4])
colnames(relative_center_ALLCNV)=colnames(center_ALLCNV)



CELL_TAG=read.delim('cell_tag.txt',header=F,sep='\t')
GENE_TAG=read.delim('gene_tag.txt',header=F,sep='\t')
CELL_COLOR= as.character(CELL_TAG[,2])
CHROM_COLOR= as.character(GENE_TAG[,2][(WINDOW/2+1):(WINDOW/2+length(center_ALLCNV[,1]))])
CHROM_TAG= GENE_TAG[(WINDOW/2+1):(WINDOW/2+length(center_ALLCNV[,1])),]

write.table(t(relative_center_ALLCNV),file='relative_ALLCNV.txt',row.names=T,col.names=T,quote=F,sep='\t')
write.table(CELL_COLOR, file='CELL_COLOR.txt',row.names=F,col.names=F,quote=F,sep='\t')
write.table(CHROM_COLOR, file='CHROM_COLOR.txt',row.names=F,col.names=F,quote=F,sep='\t')
write.table(CHROM_TAG, file='CHROM_TAG.txt',row.names=F,col.names=F,quote=F,sep='\t')


