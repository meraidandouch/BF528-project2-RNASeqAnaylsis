#-----------
#title:"BF528 Project_2 Biologist"
#author:"Qinrui Wu"
#date:"03/01/2021"
#-----------

#read library
library(ggplot2)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

#----------------------7.1--------------------------------

#read files
P0_1 <- read.csv('/projectnb/bf528/users/hedgehog_2022/project2/Programmer-Project2/cuffdiff_out/genes.fpkm_tracking', sep="\t")
fpkm <- read.csv('/projectnb/bf528/project_2/data/fpkm_matrix.csv', sep="\t")

#read all tables as tibble and select two cols; tracking_id, gene_short_name
Ad <- as_tibble(read.csv("/projectnb/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking",sep="\t")) %>% 
  select(tracking_id, gene_short_name)

#P0_1 with appropriate values
P0_1 <- P0_1 %>%
  select(tracking_id,gene_short_name, P0_FPKM) %>%
  rename(P0_1_FPKM = P0_FPKM)
P0_1 <- aggregate(P0_1_FPKM ~ gene_short_name, P0_1, mean)

#fpkm_all contains all 8 samples with values
fpkm_all <- merge(fpkm, Ad, by="tracking_id") %>%
  select(!tracking_id)
fpkm_all <- aggregate(. ~ gene_short_name, fpkm_all, mean)
fpkm_all <- merge(fpkm_all, P0_1, by="gene_short_name") %>% relocate(P0_1_FPKM, .after = Ad_2_FPKM)


#make a new tibble that contains tracking_id and corresponding mean value for each type of gene, then merge Ad table 
#with short gene name, then pivot_long the tibble to plot it.
plot_fpkm <- tibble(gene_short_name = fpkm_all$gene_short_name,
          Ad = (fpkm_all$Ad_1_FPKM+fpkm_all$Ad_2_FPKM)/2,
          P4 = (fpkm_all$P4_1_FPKM+fpkm_all$P4_2_FPKM)/2,
          P7 = (fpkm_all$P7_1_FPKM+fpkm_all$P7_2_FPKM)/2,
          P0 = (fpkm_all$P0_1_FPKM+fpkm_all$P0_2_FPKM)/2
          ) %>%
          pivot_longer(cols=c("P4", "P7", "Ad","P0"), names_to = 'type', values_to='value')

#gene lists
Sarcomere_genes <- c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
Mitochondria_genes <- c('Mpc1', 'Prdx3', 'Acat1', 'Echs1', 'Slc25a11', 'Phyh')
CellCycle_genes <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")
level_order = c("P0", "P4", "P7", "Ad")

#create a function to help plot praphs
plot_genes <- function(title, gene_list){
  temp <- filter(plot_fpkm, gene_short_name %in% gene_list)
  ggplot(temp, aes(x=factor(type, level=level_order),y=value,colour=factor(gene_short_name,level=gene_list),group=gene_short_name))+
    geom_line()+
    geom_point()+
    labs(title = title, x='Samples', y="FPKM", colour="Genes")
}

#plot three graphs
sarc <- plot_genes('Sarcomere', Sarcomere_genes)
mito <- plot_genes('Mitochondria', Mitochondria_genes)
cc <- plot_genes('CellCycle', CellCycle_genes)

#combine all plots
grid.arrange(grobs = list(sarc,mito, cc))


#----------------------7.3--------------------------------
#read file
gene_exp <- as_tibble(read.csv('/projectnb/bf528/users/hedgehog_2022/project2/Programmer-Project2/cuffdiff_out/gene_exp.diff', sep="\t"))

#get significant genes
subset <- filter(gene_exp, significant == 'yes') %>%
  arrange(desc(q_value)) %>%
  slice_head(n=1000) %>%
  select(gene)

#filtered fpkm matrix
filter_fpkm <- filter(fpkm_all, gene_short_name %in% subset$gene) %>%
  rename_all(funs(str_replace_all(., "_FPKM", " ")))


#colors for heatmap
my_colors <- colorRampPalette(c("blue","black", "yellow"), interpolate = "spline")
rownames(filter_fpkm)<-filter_fpkm[,1]

#plot heatmap 
my_heatmap <- heatmap(as.matrix(filter_fpkm[, 2:9]),
                        scale = "row", 
                        col = my_colors(10),
                        xlab = "Sample name", 
                        ylab = "Gene ID"
                        ) 




