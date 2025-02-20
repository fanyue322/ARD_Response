library(Seurat)
library(scRNAtoolVis)
library(ggplot2)
setwd("Specify your work directory")
MarkerPlot <- function(object, marker_gene, group) {
  plot <- DotPlot(
    object, group.by = group,
    features = marker_gene
  ) +
    scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3") +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
    theme(axis.text.y = element_text(size = 10, face = "italic")) +
    theme(legend.title = element_text(face = "bold", size = 12)) +
    coord_flip() +
    xlab("") +
    ylab("")
  return(plot)
}

## our markers
marker_gene<-c("CD74","GPR183","CXCL8","IFI27","IFI16","COL1A1","COL1A2","MMP2","WWP2","PLCG2","COL27A1","COL10A1","IBSP","SPP1",
               "CRTAC1","ABI3BP","PRG4","HSPA1B","HSPA1A","HSPA6","DDIT3","JUN","OGN","CILP","CILP2","CHI3L2","CHI3L1","CYTL1","FRZB","CHRDL2","HMGA1","BMP2","C11orf96")

object <- readRDS("/data/public/OA/hannah/all.dim_30_res_0.9.harmony.rds")
object$group2 <- seurat$group
levels(object$group2) <- c("this study","Hannah et al." )
DefaultAssay(object) <- "RNA"
idx <- which(object$seurat_clusters=="20")
seurat <- seurat[,-idx]
cell.levels <- c("InfC","preInfC","FC","preFC","HTC","preHTC","HomC","RepC","RegC","EC","ProC")
color<-c("#717200","#763931","#2561A0","#F8C7B4","#DF5406","#2EB9CC","#C657B3","#ADD187",
         "#159931","#CB2B1F","#E99407")# "#EF856A","#1864AA"
cell.color <- setNames(color, cell.levels)
object<-RenameIdents(object,
                     "0"="preHTC",
                     "1"="FC",
                     "2"="EC",
                     "3"="RegC",
                     "4"="RegC",
                     "5"="RepC",
                     "6"="HomC",
                     "7"="EC",
                     "8"="preHTC",
                     "9"="ProC",
                     "10"="preHTC",
                     "11"="HTC", 
                     "12"="EC",
                     "13"="preFC",
                     "14"="RegC",
                     "15"="InfC",
                     "16"="HomC",
                     "17"="FC",
                     "18"="preInfC",
                     "19"="FC")
object$celltype<-object@active.ident
levels <- c("ProC","EC","RegC","RepC","HomC","preHTC","HTC","preFC","FC","preInfC","InfC")
object$celltype <- factor(object$celltype,levels = levels)
###Figure 1A & B
p1<-DimPlot(object,pt.size = 0.5,label = F,label.size = 5,group.by = "celltype",raster = F,cols = rev(cell.color) )
p2<-DimPlot(object,pt.size = 0.5,label = F,label.size = 5,group.by = "group2",raster = F,cols = c("#DFC86B","#6176B5"))
###Figure 1C###
p3<-MarkerPlot(seurat,marker_gene,group = "celltype")

#################Marker Visualization####################
GeneVlnPlot <- function(object, marker_gene, group) {
 plot<-VlnPlot(object,features=marker_gene,group.by=group,assay='RNA',slot='data',pt.size = 0)+
          geom_boxplot(width=0,fill="white")+
          theme(axis.text = element_text(size = 25,face="bold"), 
                axis.title = element_text(size = 25,face="bold"),
                legend.text = element_text(size = 25,face="bold"),
                plot.title = element_text(size=25,face="bold"),
                legend.position='none')+
          xlab('Cell subtype')
  return(plot)
}

gene.to.plot=c("IFI27","PECAM1","CD34","CDH5")
obj_ard<-object[,which(object$group=='ARD')]
for(gene in gene.to.plot)
{
  p<-GeneVlnPlot(obj_ard,gene,'celltype')
  print(p)
}
