require(factoextra) # PCA
require(gridExtra) # Visu
require(cluster) # PAM algo

# Load datasets
load("processed_data/meta.rda")
load("processed_data/non_cible.rda")

# Out directory
outdir = "inputs_sansCTRL2"

# Keep HFpEF samples
meta = meta[which(meta$pheno == "HFpEF"),]
non_cible = non_cible[rownames(meta),]

# Only identified lipids
id = readxl::read_xlsx("data/BECAME1_LipidomicData_lipidesidentifiés.xlsx")
id = as.data.frame(id)
rownames(id) = id$Compound.Name
non_cible = non_cible[,id$Compound.Name]

# PCA for non-cible data
res_pca <- prcomp(non_cible, scale = F, center = T)
p1=fviz_eig(res_pca)
p2=fviz_pca_ind(res_pca, axes = c(1,2),
                col.ind = as.factor(meta$Diabete), label='none',
                legend.title = "Diabete")
p3=fviz_pca_ind(res_pca, axes = c(1,3),
                col.ind = as.factor(meta$Diabete), label = "none",
                legend.title = "Diabete")
p4=fviz_pca_ind(res_pca, axes = c(1,4),
                col.ind = as.factor(meta$Diabete), label = "none",
                legend.title = "Diabete")
p5=fviz_pca_var(res_pca,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE, geom.var = c("point"))
p6=fviz_pca_var(res_pca, axes = c(1,3),
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE, geom.var = c("point"))
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2, top="PCA HFpEF Non-Cible data (scale=F)")

############################# H-CLUST
# # H-CLUST ON RAW DATA
# dist = dist(non_cible, method = 'euclidean')
# hclustBrut = hclust(dist)
# plot(hclustBrut)
# 
# for(k in 2:5){
#   cl = cutree(hclustBrut, k=k)
#   print(table(cl))
#   cl = data.frame("Patient"=names(cl), "Cluster"=unname(cl))
#   write.table(cl, file=paste0(outdir, "/HFpEF_NCBrut_HClust", k, ".clst"), row.names = F, col.names = T, sep='\t', quote=F)
# }

# H-CLUST ON PCA DATA
summ <- summary(res_pca)
cum = cumsum(summ$importance[2,])
# number of PC to explain >80% variance
cum80 = length(cum)-length(cum[which(cum>=0.80)])+1

distPCA = dist(res_pca$x[,c(1:cum80)], method = 'euclidean')
hclustPCA = hclust(distPCA)
plot(hclustPCA)

for(k in 2:5){
  cl = cutree(hclustPCA, k=k)
  print(table(cl))
  cl = data.frame("Patient"=names(cl), "Cluster"=unname(cl))
  write.table(cl, file=paste0(outdir, "/HFpEF_NCPCA_HClust", k, ".clst"), row.names = F, col.names = T, sep='\t', quote=F)
}


############################# PAM and KMEANS
# PAM and KMEANS ON RAW DATA
for(k in 2:5){
  pam_clusters = pam(non_cible,
                     diss = FALSE,
                     nstart = 1000,
                     metric = "euclidean",
                     k = k)

  outPam = paste(outdir, "/HFpEF_NCBrut_PAM", k, ".clst", sep="")
  pamCl = data.frame("Patient"=names(pam_clusters$clustering), "Cluster"=pam_clusters$clustering)
  write.table(pamCl, file=outPam, sep="\t", col.names=T, row.names=F, quote=F)
  print(table(pamCl$Cluster))

  kmeans_clusters = kmeans(non_cible, k, nstart = 1000)
  outKmeans = paste(outdir, "/HFpEF_NCBrut_Kmeans", k, ".clst", sep="")
  kmeansCl = data.frame("Patient"=names(kmeans_clusters$cluster), "Cluster"=kmeans_clusters$cluster)
  write.table(kmeansCl, file=outKmeans, sep="\t", col.names=T, row.names=F, quote=F)
  print(table(kmeansCl$Cluster))
}

# PAM and KMEANS ON PCA DATA
for(k in 2:5){
  pam_clusters = pam(res_pca$x[,c(1:cum80)],
                     diss = FALSE,
                     nstart = 1000,
                     metric = "euclidean",
                     k = k)
  outPam = paste(outdir, "/HFpEF_NCPCA_PAM", k, ".clst", sep="")
  pamCl = data.frame("Patient"=names(pam_clusters$clustering), "Cluster"=pam_clusters$clustering)
  write.table(pamCl, file=outPam, sep="\t", col.names=T, row.names=F, quote=F)
  print(table(pamCl$Cluster))

  kmeans_clusters = kmeans(res_pca$x[,c(1:cum80)], k, nstart = 1000)
  outKmeans = paste(outdir, "/HFpEF_NCPCA_Kmeans", k, ".clst", sep="")
  kmeansCl = data.frame("Patient"=names(kmeans_clusters$cluster), "Cluster"=kmeans_clusters$cluster)
  write.table(kmeansCl, file=outKmeans, sep="\t", col.names=T, row.names=F, quote=F)
  print(table(kmeansCl$Cluster))
}



