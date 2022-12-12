install.packages("ggplot2")
library(ggplot2)
install.packages("ggsci")
library(ggsci)


econ = economics
dia = diamonds
mpg = mpg 

chromosomes = read.table("../bioinfo. visualization/Table1.txt")

matrix = read.table("../bioinfo. visualization/GSE46224_Yang_et_al_human_heart_RNASeq.txt")

# give me information about my data
str(dia)
# to display all columns of diamond data
colnames(dia)
# ggplot 
# igv pallet have 50 color 
ggplot(data = chromosomes , mapping = aes( x = chrom , fill = type)) + geom_bar()+
  theme_light() + scale_fill_manual(values = pal_igv()(15))

#================================== histogram ================================#

# xlab change the x axis name 
# y lab change the y axis name
ggplot(data = dia , mapping = aes(x=price)) + 
  geom_histogram(bins = 10 ,fill="red" , color="black") +
  theme_light()+ 
  xlab("Price of diamond")+
  ylab("count of diamond")+
  ggtitle("this is histogram")

#==================================scatter plot===============================#

library(ggplot2)
library(ggsci)

str(dia)

p1= ggplot(data = dia , aes(x=carat , y=price , color=cut)) +
  geom_point()+
  theme_bw()+
  xlab("Carat")+
  ylab("Price")+
  ggtitle("relationship")

p2= ggplot(data = dia , aes(x=carat , y=price , color=cut)) +
  geom_point()+
  theme_bw()+
  xlab("Carat")+
  ylab("Price")+
  ggtitle("relationship")+
  facet_grid(~cut) # split each cut value

#=====================================box plot=================================#

mpg = mpg 
str(mpg)

ggplot(data = mpg , aes(drv,cty, fill=trans)) + 
  geom_boxplot()+
  theme_bw()+
  theme_light()


#=======================================violin plot============================#

library(ggplot2)
library(ggsci)

mpg=mpg

head(mpg)

ggplot(data = mpg , aes(drv,cty)) + geom_violin(color = "yellow" , fill="light blue") +
  theme_light()+ theme_bw() + geom_boxplot( width = 0.2 , color="red" , fill= "blue")+
  coord_flip() # flip the row with col 


#==================================heat map=====================================#  

install.packages("pheatmap")
install.packages("RColorBrewer")

library(pheatmap)
library(RColorBrewer)

data = read.table("../bioinfo. visualization/SRP029880.raw_counts.tsv")

ann = read.table("../bioinfo. visualization/SRP029880.colData.tsv" , sep = "\t")

# preprocessing of data 

datacut = data[,-11]  # remove last col
head(datacut)

# convert data as matrix
class(datacut)
matrix_data = as.matrix(datacut) 
class(matrix_data)

# margin = 1 -> indicate rows
# margin = 2 -> indicate columns
dim(matrix_data)
# select 100 variant genes from  19719

most_var_genes = apply( X = matrix_data , MARGIN = 1 ,FUN = var)
selected_genes = names(most_var_genes[order(most_var_genes , decreasing = T)][1:100]) 
final_data = matrix_data[selected_genes,]

class(final_data)

# visualization of heatmap
# main() -> used for title
pheatmap(mat = final_data , scale = "row" , cluster_rows = T  , cluster_cols = T ,
         cutree_rows = 2  , cutree_cols = 2 , show_rownames = T , fontsize_row = 5,
         show_colnames = T , main = "this is the first heatmap" , 
         color =  brewer.pal(n=9 , name = "PuBu")  , annotation_col = ann )
# if we use ggsci instead of brewer paler ->  ggsci::pal_aaas()(10)

# to display the palet of colors of RColorBrewer
RColorBrewer::display.brewer.all()
RColorBrewer::brewer.pal.info # all information about the colorpalet
# show the colors of single palet  # n -> show the number of colors
RColorBrewer::display.brewer.pal(n=11 , name = "BrBG")
# show the names of all colors in this palet
RColorBrewer::brewer.pal(n=12 , name = "Paired")

#==================================== PCA =====================================#
# Dimensionality reduction -> principle component analysis

paste("lama","medhat") # put spaces between names
paste0("lama" , "medhat")  # don't put spaces between names

for(i in 1:5){
  print(paste("num. is" , i))
}





data_matrix = matrix(nrow=100 , ncol=10)
colnames(data_matrix) = c(
  paste("normal" , 1:5 , sep=""),
  paste0("ubnormal" , 1:5 , sep="")
)
rownames(data_matrix) = paste("gene" , 1:100 , sep="")
head(data_matrix)

for(i in 1:100){
  normal.values = rpois(5 , lambda = sample(x=10:1000 , size = 1))
  ubnormal.values = rpois(5 , lambda = sample(x=10:1000 , size = 1))
  data_matrix[i,] = c( normal.values , ubnormal.values)
}

head(data_matrix)
data_matrix
dim(data_matrix)


# t -> mean transpose
# i make transpose to see the variation in samples (genes)
# prcomp this is a function of PCA

my_pca = prcomp(t(data_matrix), scale = T)
my_pca

# plot PCA 1 -> most variation in data
# plot PCA 2 -> second most variation in data

# number of PCAs = number of sample in your data

plot(my_pca$x[,1] , my_pca$x[,2] )

# i will compute the variation based on standard deviation of data

pca.var = my_pca$sdev^2

# compute the percentage of variation

pca.var.per = round(pca.var/sum(pca.var)*100 , 1)
barplot(pca.var.per, main = "screen plot" , xlab = "principle component" , 
        ylab = "percent variation")





































