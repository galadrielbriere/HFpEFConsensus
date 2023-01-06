# Clinical metadata
meta = read.table("data/clinicalMetadata_revuGalaNoComma.csv", header=T, na.strings = c("", " ", "NA"), sep=";")
# Use ID_JTL as id for samples
meta$ID_JTL = paste0("JTL", meta$ID_JTL)
rownames(meta) = meta$ID_JTL
# Format echLipidomique ID
meta$echLipidomique = paste0("X",meta$echLipidomique)

# Cible data
cible = read.table("data/became_cible_ctrl_hfpef.csv", sep="\t", header=T)
cible$ID_JTL = paste0("JTL", cible$ID_JTL)
# extract metadata for cible samples
metacible = cible[,c(1:29)]
cible = cible[,-c(1:29)]
table(metacible$ID_JTL %in% rownames(meta))
# Use ID_JTL as id
rownames(metacible) = metacible$ID_JTL
rownames(cible) = metacible$ID_JTL

# Non-cible data
non_cible = read.table("data/BECAME1_FBF_corrige_non-cible.csv", header=T, sep="\t")
rownames(non_cible) = non_cible$Compound.Name
non_cible = non_cible[,-c(1:3)]
table(colnames(non_cible) %in% meta$echLipidomique)
# 25 samples to remove from non_cible data
non_cible = non_cible[,which(colnames(non_cible) %in% meta$echLipidomique)]

# Convert echLipidomique id to ID_JTL
id_jtl = meta$ID_JTL
names(id_jtl) = meta$echLipidomique
colnames(non_cible) = id_jtl[colnames(non_cible)]

# Transpose to have samples as rows
non_cible = as.data.frame(t(non_cible))

# Log2 
non_cible = log2(non_cible)

# Save processed data
save(meta, file="processed_data/meta.rda")
save(cible, file="processed_data/cible.rda")
save(metacible, file="processed_data/metacible.rda")
save(non_cible, file="processed_data/non_cible.rda")
