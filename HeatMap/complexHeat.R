
library('ComplexHeatmap')
library('circlize')
library('seriation')
set.seed(1)

mat=COM
METHOD='ARSA'
o1 = get_order(seriate(dist(mat), method = METHOD))

o2=order(SCORE)

#o2 = get_order(seriate(dist(t(mat)), method = METHOD))

#SUM=apply(log(NTAB+1,10),2,sum)
#o2=order(SUM)
o.mat=mat[o1,o2]
#col_fun =colorRamp2(c(0,10,50,100 ), c('royalblue3','yellow','red3','red3'))
col_fun =colorRamp2(c(0,50,100 ), c('grey98','red3','red3'))





ha = HeatmapAnnotation(
	Pred.PC2 = anno_lines(SCORE[o2], add_points = TRUE) ,
    #foo = 1:length(LABEL), 
    Type = TYPE[o2],
    H3K27M.mut=GPV.ALL[o2],
    TP53.mut=TP53V.ALL[o2],
    Survival=SURV.ALL[o2],
    Response=RESP.ALL[o2],
    
	#PC2.bulk=PC2V.ALL[o2],
    col = list(
	       Type = c("Autopsy"='gold','Biopsy'='royalblue3','CellLine'='grey90'),
          Survival=c('Long'='green','Short'='red','NA'='grey90'),
	    Response=c('Good'='green','Bad'='red','NA'='grey90'),
	    H3K27M.mut=c('H3.3'='hotpink1','H3.1'='purple','WT'='grey90'),
	    TP53.mut=c('TP53'='purple','WT'='hotpink1','NA'='grey90'),
	    BMI1.exp=c('Low'='green','High'='red','NA'='grey90')
	    #PRC2.Target.exp=c('Low'='green','High'='red','NA'='grey90')
	    #PC2.bulk=c('Low'='green3','High'='red','NA'='grey80')
          #Response=c('Good'='olivedrab2','Bad'='red','NA'='grey80')
	      ),
    gp = gpar(col = "black"),
    
    BMI1.exp=EXPV.ALL[o2]
    #PRC2.Target.exp=EXPV1.ALL[o2]
    #Real.PC2.bulk=anno_lines(Real.PC2[o2], add_points = TRUE)
    
)


Heatmap(o.mat,row_title='',name="%",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=TRUE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	bottom_annotation = ha
	#cell_fun = function(j, i, x, y, width, height, fill) {
        #grid.text(sprintf("%.0f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
        #}
	)


