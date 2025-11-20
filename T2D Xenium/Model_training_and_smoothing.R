##############################################################
rm(list=ls())
library(Seurat)
library(hdf5r)
library(DEGAS)
library(dplyr)
library(ggplot2)
library(progressr)
library(progress)
library(ggrepel)
source("~/DEGAS_functions_for_smoothing.R")

################################# DEGAS model training##############

#Options
save_input = TRUE
save_output = TRUE
retrain = TRUE

# Loading bulk T2D data
Xpat = read.csv("~/bulkRNA_seq_Diabetes.csv",row.names=1)
Ypat = read.csv("~/patLab.csv",row.names=1)
Ypat = Ypat[colnames(Xpat),]


# Loading Xenium data
#######************ The Xenium T2D data is available in https://www.synapse.org/Synapse:syn68699752/files/ ***********
S0017155_ND = LoadXenium("~/output-XETG00126__0017155__ND16091_ND__20231010__204554")
S0017155_T2D = LoadXenium("~/output-XETG00126__0017155__ND16154_T2D__20231010__204554")
S0017308_T2D = LoadXenium("~/output-XETG00126__0017308__ND16090_T2D__20231010__204554")
S0017308_ND = LoadXenium("~/output-XETG00126__0017308__ND16091_ND__20231010__204554")
# Formatting Xenium data
nd1 = S0017155_ND@assays$Xenium$counts
rownames(nd1) = rownames(S0017155_ND)
colnames(nd1) = sub("-1","-ND1",colnames(S0017155_ND))
t2d1 = S0017155_T2D@assays$Xenium$counts
rownames(t2d1) = rownames(S0017155_T2D)
colnames(t2d1) = sub("-1","-T2D1",colnames(S0017155_T2D))
t2d2 = S0017308_T2D@assays$Xenium$counts
rownames(t2d2) = rownames(S0017308_T2D)
colnames(t2d2) = sub("-1","-T2D2",colnames(S0017308_T2D))
nd2 = S0017308_ND@assays$Xenium$counts
rownames(nd2) = rownames(S0017308_ND)
colnames(nd2) = sub("-1","-ND2",colnames(S0017308_ND))
Xst = cbind(nd1,t2d1,nd2,t2d2)

# Defining final matrices for training
feats1 = intersect(rownames(Xpat),rownames(Xst))
stDat = Xst[feats1,]
stDat = preprocessCounts(as.matrix(stDat))
stDat = as.matrix(stDat)
patDat = Xpat[feats1,]
patDat = preprocessCounts(as.matrix(patDat))
patLab = as.matrix(Ypat[rownames(patDat),])
rm(S0017155_ND,S0017155_T2D,S0017308_ND,S0017308_T2D,nd1,nd2,t2d1,t2d2,Xst)
if(save_input){
  saveRDS(stDat,file="~/Desktop/Xenium/stDat.rds")
}

# Training DEGAS model

model1 = runDEGASatlas(stDat,NULL,patDat,patLab,"~/Desktop/Xenium/tmp/","BlankClass","DenseNet",3,5,20000,123)
saveRDS(model1,file="~/Desktop/Xenium/model1.rds") #saving the DEGAS model
preds = predClassBag(model1,stDat,"pat")
saveRDS(preds,file="~/Desktop/Xenium/preds1.rds") #saving the predicted values


group = sub(".*[-]","",rownames(stDat))
mapper = c("ND","ND","T2D","T2D")
names(mapper) = c("ND1","ND2","T2D1","T2D2")
group2 = mapper[group]

# Loading UMAP coordinates
umap_nd1 = read.csv("~/output-XETG00126__0017155__ND16091_ND__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_nd1$Barcode = sub("-1","-ND1",umap_nd1$Barcode)
umap_t2d1 = read.csv("~/output-XETG00126__0017155__ND16154_T2D__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_t2d1$Barcode = sub("-1","-T2D1",umap_t2d1$Barcode)
umap_nd2 = read.csv("~/output-XETG00126__0017308__ND16090_T2D__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_nd2$Barcode = sub("-1","-T2D2",umap_nd2$Barcode)
umap_t2d2 = read.csv("~/output-XETG00126__0017308__ND16091_ND__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_t2d2$Barcode = sub("-1","-ND2",umap_t2d2$Barcode)
umap_coords = rbind(umap_nd1,umap_t2d1,umap_nd2,umap_t2d2)
rownames(umap_coords) = umap_coords$Barcode
umap_coords = umap_coords[rownames(stDat),]

# Loading cellular locations
locs_nd1 = read.csv("~/output-XETG00126__0017155__ND16091_ND__20231010__204554/cells.csv.gz")
locs_nd1$cell_id = sub("-1","-ND1",locs_nd1$cell_id)
locs_t2d1 = read.csv("~/output-XETG00126__0017155__ND16154_T2D__20231010__204554/cells.csv.gz")
locs_t2d1$cell_id = sub("-1","-T2D1",locs_t2d1$cell_id)
locs_nd2 = read.csv("~/output-XETG00126__0017308__ND16090_T2D__20231010__204554/cells.csv.gz")
locs_nd2$cell_id = sub("-1","-T2D2",locs_nd2$cell_id)
locs_t2d2 = read.csv("~/output-XETG00126__0017308__ND16091_ND__20231010__204554/cells.csv.gz")
locs_t2d2$cell_id = sub("-1","-ND2",locs_t2d2$cell_id)
locs_coords = rbind(locs_nd1,locs_t2d1,locs_nd2,locs_t2d2)
rownames(locs_coords) = locs_coords$cell_id
locs_coords = locs_coords[rownames(stDat),]

# Creating an overall dataframe for smoothing
df=data.frame(barcodes=umap_coords$Barcode,t2d_risk=preds[,1],umap1=umap_coords$UMAP.1,umap2=umap_coords$UMAP.2,locsx=locs_coords$x_centroid,locsy=locs_coords$y_centroid)
df$sample = sub(".*[-]","",df$barcodes)
rm(umap_nd1,umap_t2d1,umap_nd2,umap_t2d2,locs_nd1,locs_t2d1,locs_nd2,locs_t2d2)
gc()

df = na.omit(df)
saveRDS(df,file="~/df_v2.rds")


################################################################################
df = readRDS("~/dfv2.rds")

# Creating plots--------------
smoothed_nd1 = readRDS("~/smoothed_nd1.rds")
smoothed_nd2 = readRDS("~/smoothed_nd2.rds")
smoothed_t2d1 = readRDS("~/smoothed_t2d1.rds")
smoothed_t2d2 = readRDS("~/smoothed_t2d2.rds")
smoothed2 = rbind(smoothed_nd1,smoothed_t2d1,smoothed_nd2,smoothed_t2d2)
smoothed2_subset <- smoothed2[which(smoothed2$barcodes%in%df$barcodes),]
df_subset <- df[which(df$barcodes%in%smoothed2$barcodes),]

df_combined <- merge(df_subset,smoothed2_subset,by="barcodes")
df_tmp = df_combined[!is.na(df_combined$smoothed),]
df_tmp$scaled_haz = toCorrCoeff(df_tmp$smoothed) #scale it to a [-1,1] value

group = sub(".*[-]","",df_tmp$barcodes)
mapper = c("ND","ND","T2D","T2D")
names(mapper) = c("ND1","ND2","T2D1","T2D2")
group2 = mapper[group]
df_tmp$group = factor(group,levels=c("ND1","T2D1","ND2","T2D2"))
df_tmp$group2 = group2


# Sample p-values calculation
p_vals <- df_tmp %>%
  filter(group %in% c("ND1", "T2D1")) %>%
  summarize(p = wilcox.test(scaled_haz ~ group)$p.value) %>%
  mutate(label = ifelse(p < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(p, 3))),
         xpos = 1.5,  # Midpoint between ND1 (1) and T2D1 (2)
         ypos = 1)

p_vals2 <- df_tmp %>%
  filter(group %in% c("ND2", "T2D2")) %>%
  summarize(p = wilcox.test(scaled_haz ~ group)$p.value) %>%
  mutate(label = ifelse(p < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(p, 3))),
         xpos = 3.5,  # Midpoint between ND2 (3) and T2D2 (4)
         ypos = 1)

p_df <- bind_rows(p_vals, p_vals2)

# Your plot with ggrepel
boxplot1 <- ggplot(df_tmp, aes(x = group, y = scaled_haz, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )+
  ylab("T2D Association") +
  xlab("Sample") +
  scale_fill_manual(values = c("ND1" = "blue", "T2D1" = "red", "ND2" = "blue", "T2D2" = "red")) +
  geom_text_repel(data = p_df, aes(x = xpos, y = ypos, label = label), inherit.aes = FALSE)

ggsave("~/boxplot_bysamples_scaled.pdf",boxplot1)


# disease p-values calculation
p_vals3 <- df_tmp %>%
  filter(group2 %in% c("ND", "T2D")) %>%
  summarize(p = wilcox.test(scaled_haz ~ group2)$p.value) %>%
  mutate(label = ifelse(p < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(p, 3))),
         xpos = 1.5,  # Midpoint between ND(1) and T2D(2)
         ypos = max(df_tmp$scaled_haz[df_tmp$group2 %in% c("ND", "T2D")]) * 1.05)


boxplot2 <- ggplot(df_tmp, aes(x = group2, y = scaled_haz, fill = group2)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )+
  ylab("T2D Association") +
  xlab("") +
  scale_fill_manual(values = c("ND" = "blue", "T2D" = "red")) +
  geom_text_repel(data = p_vals3, aes(x = xpos, y = ypos, label = label), inherit.aes = FALSE)

ggsave("~/boxplot_bydisease_scaled.pdf",boxplot2)



# heatmap for samples---------

df_tmp_ND1 = df_tmp[df_tmp$sample=="ND1",]
pdf("~/ND1_sample_heatmap.pdf",width=6,height=5)

ggplot(df_tmp_ND1, aes(x = locsx, y = locsy, colour = scaled_haz)) + 
  geom_point(alpha = 0.5, size = 0.1) + 
  scale_colour_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = median(df_tmp_ND1$scaled_haz),
    name = "T2D\nAssociation"  # Legend title
  ) +
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle(expression(ND[1]))+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
dev.off()

df_tmp_T2D1 = df_tmp[df_tmp$sample=="T2D1",]
pdf("~/T2D1_sample_heatmap.pdf",width=6,height=5)

ggplot(df_tmp_T2D1, aes(x = locsx, y = locsy, colour = scaled_haz)) + 
  geom_point(alpha = 0.5, size = 0.1) + 
  scale_colour_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = median(df_tmp_T2D1$scaled_haz),
    name = "T2D\nAssociation"  # Legend title
  ) +
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle(expression(T2D[1]))+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
dev.off()


df_tmp_ND2 = df_tmp[df_tmp$sample=="ND2",]
pdf("~/ND2_sample_heatmap.pdf",width=6,height=5)

ggplot(df_tmp_ND2, aes(x = locsx, y = locsy, colour = scaled_haz)) + 
  geom_point(alpha = 0.5, size = 0.1) + 
  scale_colour_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = median(df_tmp_ND2$scaled_haz),
    name = "T2D\nAssociation"  # Legend title
  ) +
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle(expression(ND[2]))+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
dev.off()


df_tmp_T2D2 = df_tmp[df_tmp$sample=="T2D2",]
pdf("~/T2D2_sample_heatmap.pdf",width=6,height=5)

ggplot(df_tmp_T2D2, aes(x = locsx, y = locsy, colour = scaled_haz)) + 
  geom_point(alpha = 0.5, size = 0.1) + 
  scale_colour_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = median(df_tmp_T2D2$scaled_haz),
    name = "T2D\nAssociation"  # Legend title
  ) +
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle(expression(T2D[2]))+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
dev.off()


